
subroutine dforce( NGdof, position, force, LComputeDerivatives, LComputeAxis)



  use constants, only : zero, half, one, pi, pi2

  use numerical, only : logtolerance

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdforce, Nvol, Ntor, Lrad, Igeometry, &
                        epsilon, &
                        Lconstraint, Lcheck, dRZ, &
                        Lextrap, &
                        mupftol, &
                        LHmatrix

  use cputiming, only : Tdforce

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        Mvol, NAdof, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        ImagneticOK, &
                        Energy, ForceErr, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        mn, im, in, &
                        dpflux, dtflux, sweight, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        BBe, IIo, BBo, IIe, & ! these are just used for screen diagnostics;
                        LGdof, dBdX, &
                        Ate, Aze, Ato, Azo, & ! only required for broadcasting
                        diotadxup, dItGpdxtp, & ! only required for broadcasting
                        lBBintegral, &
                        dFFdRZ,HdFFdRZ, dBBdmp, dmupfdx, &
                        hessian, dessian, Lhessianallocated, &
                        hessian2D,dessian2D, Lhessian2Dallocated, &
                        hessian3D,dessian3D,Lhessian3Dallocated, denergydrr, denergydrz,denergydzr,denergydzz, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        LocalConstraint, xoffset, &
                        solution, IPdtdPf, &
                        IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                        ext ! For outputing Lcheck = 6 test



  LOCALS

  INTEGER, parameter   :: NB = 3 ! optimal workspace block size for LAPACK:DSYSVX;

  INTEGER, intent(in)  :: NGdof               ! dimensions;
  REAL,    intent(in)  :: position(0:NGdof)
  REAL,    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives !

  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, tdoc_ntz, idoc, idof, tdof, jdof, ivol, imn, ll, ihybrd1, lwa, Ndofgl, llmodnp
  INTEGER              :: maxfev, ml, muhybr, mode, nprint, nfev, ldfjac, lr, Nbc, NN, cpu_id, ideriv
  REAL                 :: epsfcn, factor
  REAL                 :: Fdof(1:Mvol-1), Xdof(1:Mvol-1)
  INTEGER              :: ipiv(1:Mvol)
  REAL, allocatable    :: fjac(:, :), r(:), Fvec(:), dpfluxout(:)

  INTEGER              :: status(MPI_STATUS_SIZE), request_recv, request_send, cpu_send
  INTEGER              :: id
  INTEGER              :: iflag, idgesv, Lwork
  INTEGER              :: idofr,idofz,tdofr,tdofz

  CHARACTER            :: packorunpack
  EXTERNAL             :: dfp100, dfp200

  LOGICAL              :: LComputeAxis, dfp100_logical



  BEGIN(dforce)

  if( LocalConstraint ) then
    SALLOCATE( dmupfdx, (1:Mvol,    1:1,1:2,1:LGdof,0:1), zero )
  else
    SALLOCATE( dmupfdx, (1:Mvol, 1:Mvol-1,1:2,1:LGdof,1), zero ) ! TODO change the format to put vvol in last index position...
  endif
  
  Lhessianallocated = .true.
  Lhessian2Dallocated = .false.
  Lhessian3Dallocated = .false.

  if( LcomputeDerivatives ) then
    SALLOCATE( dFFdRZ, (1:LGdof,0:1,1:LGdof,0:1,1:Mvol), zero )
    SALLOCATE( dBBdmp, (1:LGdof,1:Mvol,0:1,1:2), zero )
    SALLOCATE( hessian, (1:NGdof,1:NGdof), zero )
    SALLOCATE( dessian, (1:NGdof,1:LGdof), zero )
  endif




  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;


  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                          iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack, LcomputeDerivatives, LComputeAxis ) )



  if( LcomputeDerivatives ) then

   dBBdmp(1:LGdof,1:Mvol,0:1,1:2) = zero
  endif





  Xdof(1:Mvol-1) = zero;

  if( LocalConstraint ) then

    SALLOCATE( Fvec, (1:Mvol-1), zero)

    Ndofgl = 0; Fvec(1:Mvol-1) = 0; dfp100_logical = .FALSE.;
    Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset

    dBdX%L = LComputeDerivatives
    WCALL(dforce, dfp100, (Ndofgl, Xdof, Fvec, dfp100_logical) )

    DALLOCATE( Fvec )

  else

    IPDtdPf = zero
    Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

      Ndofgl = Mvol-1


    SALLOCATE( Fvec, (1:Ndofgl), zero )

    dfp100_logical = .FALSE.

    WCALL(dforce, dfp100, (Ndofgl, Xdof(1:Mvol-1), Fvec(1:Ndofgl), dfp100_logical))

    SALLOCATE(dpfluxout, (1:Ndofgl), zero )
    if ( myid .eq. 0 ) then

        dpfluxout = Fvec
        call DGESV( Ndofgl, 1, IPdtdPf, Ndofgl, ipiv, dpfluxout, Ndofgl, idgesv )

        dpflux(2:Mvol) = dpflux(2:Mvol) - dpfluxout(1:Mvol-1)
    endif

    RlBCAST(dpfluxout(1:Ndofgl), Ndofgl, 0)
    RlBCAST(dpflux(1:Mvol)   , Mvol, 0)
    do vvol = 2, Mvol

      WCALL(dforce, IsMyVolume, (vvol))

      if( IsMyVolumeValue .EQ. 0 ) then
          cycle
      else if( IsMyVolumeValue .EQ. -1) then
          FATAL(dforce, .true., Unassociated volume)
      endif

      NN = NAdof(vvol)

      SALLOCATE( solution, (1:NN, 0:2), zero)

      packorunpack = 'P'
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! packing;
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,2), 2 ) ) ! packing;

      solution(1:NN, 0) = solution(1:NN, 0) - dpfluxout(vvol-1) * solution(1:NN, 2)


      packorunpack = 'U'
      WCALL( dforce, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! unpacking;

      DALLOCATE( solution )

    enddo ! end of do vvol = 1, Mvol

    DALLOCATE(Fvec)
    DALLOCATE(dpfluxout)


  endif !matches if( LocalConstraint )


  do vvol = 1, Mvol
    call WhichCpuID(vvol, cpu_id)

    LlBCAST( ImagneticOK(vvol)         , 1, cpu_id)

    do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
      do ii = 1, mn
        RlBCAST( Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
        RlBCAST( Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
      enddo
    enddo


    if( NOTstellsym ) then
      do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
        do ii = 1, mn
              RlBCAST( Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
              RlBCAST( Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, cpu_id)
        enddo
      enddo
    endif
  enddo



  WCALL(dforce, dfp200, ( LcomputeDerivatives, vvol) )



  do vvol = 1, Mvol

    LREGION( vvol )
    WCALL( dforce, brcast, ( vvol ) )

  enddo






  lBBintegral(1:Nvol) = lBBintegral(1:Nvol) * half

  Energy = sum( lBBintegral(1:Nvol) ) ! should also compute beta;




  ;   force(0:NGdof) = zero

  do vvol = 1, Mvol-1

    LREGION(vvol)

    tdoc = (vvol-1) * LGdof

    if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

      ;  idoc = 0           ! degree-of-constraint counter; set;

      if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! to be made redundant;
        FATAL( dforce, 2.gt.Mvol, psifactor needs attention )
        ;force(tdoc+idoc+1:tdoc+idoc+mn) = position(1:mn) - ( iRbc(1:mn,2) / psifactor(1:mn,2) )
      else
        ;force(tdoc+idoc+1:tdoc+idoc+mn    ) = ( Bemn(1:mn    ,vvol+1,0) - Bemn(1:mn    ,vvol+0,1) ) * BBweight(1:mn) ! pressure imbalance;
      endif

      ;  BBe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;

      ;  idoc = idoc + mn   ! degree-of-constraint counter; increment;

      if( Igeometry.ge.3 ) then ! add spectral constraints;

        force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = (                           Iomn(2:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                            + (                         + Somn(2:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                            - ( Somn(2:mn    ,vvol+1,0)                           ) * sweight(vvol+1)


        IIo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;

        idoc = idoc + mn-1

      endif ! end of if( Igeometry.ge.3 ) ;

      if( NOTstellsym ) then

        force(tdoc+idoc+1:tdoc+idoc+mn-1  ) = ( Bomn(2:mn    ,vvol+1,0) - Bomn(2:mn    ,vvol+0,1) ) * BBweight(2:mn) ! pressure imbalance;

        BBo(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn-1) ) ) / (mn-1), logtolerance ) ! screen diagnostics;

        idoc = idoc + mn-1 ! degree-of-constraint counter; increment;

        if( Igeometry.ge.3 ) then ! add spectral constraints;

          force(tdoc+idoc+1:tdoc+idoc+mn    ) = (                           Iemn(1:mn    ,vvol+0  ) ) * epsilon         & ! spectral constraints;
                                              + (                         + Semn(1:mn    ,vvol+0,1) ) * sweight(vvol+0) & ! poloidal length constraint;
                                              - ( Semn(1:mn    ,vvol+1,0)                           ) * sweight(vvol+1)


          IIe(vvol) = max( sum( abs( force(tdoc+idoc+1:tdoc+idoc+mn  ) ) ) / (mn  ), logtolerance ) ! screen diagnostics;

          idoc = idoc + mn   ! degree-of-constraint counter; increment;

        endif ! end of if( Igeometry.ge.3 ) ;

      endif ! end of if( NOTstellsym ) ;


    else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) );

      ;                       ; BBe(vvol) = 9.9E+09
      ;                       ; IIo(vvol) = 9.9E+09
      if ( NOTstellsym ) then ; BBo(vvol) = 9.9E+09
      ;                      ; IIe(vvol) = 9.9E+09
      endif

      ; force(tdoc+1:tdoc+LGdof) = 9.9E+09

    endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

  enddo ! end of do vvol;



  if( NGdof.ne.0 ) then ; ForceErr = sqrt( sum( force(1:NGdof)*force(1:NGdof) ) / NGdof ) ! this includes spectral constraints;
  else                  ; ForceErr = zero
  endif




4000 format("dforce : ",f10.2," : ",6x,3x,"; ",:,"|f|=",es12.5," ; ",:,"time=",f10.2,"s ;",:," log",a5,"=",28f6.2  ," ...")
4001 format("dforce : ", 10x ," : ",6x,3x,"; ",:,"    ",  12x ,"   ",:,"     ", 10x ,"  ;",:," log",a5,"=",28f6.2  ," ...")




if( Lhessian3Dallocated .and. Igeometry.ge.3) then ! construct Hessian3D;

  if (LcomputeDerivatives ) then

      hessian2D(1:NGdof,1:NGdof) = zero

      do vvol = 1, Mvol-1 ! loop over interior surfaces;

        if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

          idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry;
          idofr=0
          idofz=0

          do ii = 1, mn ! loop over degrees-of-freedom;

            do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z;


              do issym = 0, 1 ! stellarator symmetry;

                if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;

                if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
                if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};

                idof = idof + 1 ! labels degree-of-freedom;


                   if (irz.eq.0) then
                       idofr = idofr + 1 ! labels degree-of-freedom;
                   else
                       idofz = idofz + 1
                   endif
                    if ( vvol.gt.1 ) then
                       tdofr = (vvol-2) * LGdof + idofr ! labels degree-of-freedom in internal interface geometry   ;
                       tdofz = (vvol-2) * LGdof + idofz
                       tdoc = (vvol-1) * LGdof          ! labels force-balance constraint across internal interfaces; !always fix
                       idoc = 0                         ! local  force-balance constraint across internal interface ;
                       if (irz .eq.0) then
                           hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdofr) = -denergydrr(idoc+1:idoc+LGdof ,vvol+0,1,idof,0)
                       else
                           hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdofz+mn)= -denergydzr(idoc+1:idoc+LGdof,vvol+0,1,idof,0)
                       endif
                     endif ! end of if( vvol.gt.1 ) ;

                        ;tdofr = (vvol-1) * LGdof + idofr
                        ;tdofz = (vvol-1) * LGdof + idofz
                        ;tdoc = (vvol-1) * LGdof ! shorthand;
                        ;idoc = 0 ! diagonal elements;
                       if (irz .eq. 0) then
                          ;hessian2D(tdoc+1:tdoc+LGdof,tdofr) = denergydrr(idoc+1:idoc+LGdof,vvol+1,0,idof,0) - denergydrr(idoc+1:idoc+LGdof,vvol+0,1,idof,1)! &
                        else
                          ;hessian2D(tdoc+1:tdoc+LGdof,tdofz+mn) = denergydzr(idoc+1:idoc+LGdof,vvol+1,0,idof,0) - denergydzr(idoc+1:idoc+LGdof,vvol+0,1,idof,1)

                      endif


                      if ( vvol.lt.Mvol-1 ) then
                         tdofr = (vvol+0) * LGdof + idofr
                         tdofz = (vvol+0) * LGdof + idofz
                         tdoc = (vvol-1) * LGdof ! shorthand;
                         idoc = 0
                        if (irz.eq.0) then
                             hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdofr) = denergydrr(idoc+1:idoc+LGdof ,vvol+1,0,idof,1)
                        else
                             hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdofz+mn) = denergydzr(idoc+1:idoc+LGdof ,vvol+1,0,idof,1)
                        endif ! end of if( vvol.lt.Mvol-1 ) ;
                     endif







              enddo ! matches do issym ;

            enddo ! matches do irz ;
          enddo ! matches do ii ;

        else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

          FATAL( dforce, .true., need to provide suitable values for hessian2D in case of field failure )

        endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

      enddo ! end of do vvol;
   endif

endif ! end of if( LcomputeDerivatives ) ;





if( LHmatrix .and. Lhessian2Dallocated .and. Igeometry.eq.2) then ! construct Hessian2D;

  if (LcomputeDerivatives ) then

      hessian2D(1:NGdof,1:NGdof) = zero

      do vvol = 1, Mvol-1 ! loop over interior surfaces;

        if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

          idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry;

          do ii = 1, mn ! loop over degrees-of-freedom;

            do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z;

              if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;

              do issym = 0, 1 ! stellarator symmetry;

                if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;

                if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
                if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};

                idof = idof + 1 ! labels degree-of-freedom;


                if( LocalConstraint) then
                  if( vvol.gt.1 ) then
                    tdof = (vvol-2) * LGdof + idof ! labels degree-of-freedom in internal interface geometry   ;
                    tdoc = (vvol-1) * LGdof        ! labels force-balance constraint across internal interfaces;
                    idoc = 0                       ! local  force-balance constraint across internal interface ;
                    hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  - HdFFdRZ(idoc+1:idoc+LGdof,1,idof,0,vvol+0)

                  endif ! end of if( vvol.gt.1 ) ;

                  ;tdof = (vvol-1) * LGdof + idof
                  ;tdoc = (vvol-1) * LGdof ! shorthand;
                  ;idoc = 0

                  if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                      ;hessian2D(tdoc+idof                  ,tdof) = one ! diagonal elements;
                  else
                      ;hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = HdFFdRZ(idoc+1:idoc+LGdof,0,idof,0,vvol+1) - HdFFdRZ(idoc+1:idoc+LGdof,1,idof,1,vvol+0)
                  endif ! end of if( Lextrap.eq.1 .and. vvol.eq.1 )


                  if( vvol.lt.Mvol-1 ) then
                    tdof = (vvol+0) * LGdof + idof
                    tdoc = (vvol-1) * LGdof ! shorthand;
                    idoc = 0
                    if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                      if    ( im(idof).le.0                     ) then ; hessian2D(tdoc+idof,tdof) = - one
                      else                                             ; hessian2D(tdoc+idof,tdof) = - one
                      endif
                    else
                      hessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = HdFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
                    endif
                  endif ! end of if( vvol.lt.Mvol-1 ) ;


                  if( vvol.eq.Mvol-1 ) then
                    tdoc = (vvol-1) * LGdof ! shorthand ;
                    idoc = 0
                    dessian2D(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = HdFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)

                  endif ! end of if( vvol.lt.Mvol-1 ) ;

                else ! Global constraint


                  FATAL( dforce, .true., incorrect choice of Lconstraint in SPEC)

                endif ! matches if( LocalConstraint );

              enddo ! matches do issym ;

            enddo ! matches do irz ;
          enddo ! matches do ii ;

        else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

          FATAL( dforce, .true., need to provide suitable values for hessian2D in case of field failure )

        endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

      enddo ! end of do vvol;
   endif

endif ! end of if( LcomputeDerivatives ) ;







  if( LcomputeDerivatives .and. Lhessianallocated) then ! construct Hessian;


    hessian(1:NGdof,1:NGdof) = zero

    do vvol = 1, Mvol-1 ! loop over interior surfaces;

      if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

        idof = 0 ! labels degree-of-freedom = Fourier harmonic of surface geometry;



        do ii = 1, mn ! loop over degrees-of-freedom;



          do irz = 0, 1 ! Fourier harmonic of R, Fourier harmonic of Z;

            if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;

            do issym = 0, 1 ! stellarator symmetry;



              if( issym.eq.1 .and. YESstellsym ) cycle ! no dependence on the non-stellarator symmetric harmonics;

              if( ii.eq.1 .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
              if( ii.eq.1 .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};

              idof = idof + 1 ! labels degree-of-freedom;



              if( LocalConstraint ) then
                if( vvol.gt.1 ) then
                  tdof = (vvol-2) * LGdof + idof ! labels degree-of-freedom in internal interface geometry   ;
                  tdoc = (vvol-1) * LGdof        ! labels force-balance constraint across internal interfaces;
                  idoc = 0                       ! local  force-balance constraint across internal interface ;
                  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =                                           - dFFdRZ(idoc+1:idoc+LGdof,1,idof,0,vvol+0)
                  if( Lconstraint.eq.1 ) then ! this is a little clumsy; could include  or something . . . ;
                    hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                     &
                                                                - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol,1,1,idof,0) &
                                                                - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol,1,2,idof,0)
                  endif ! end of if( Lconstraint.eq.1 ) ;
                endif ! end of if( vvol.gt.1 ) ;

                ;tdof = (vvol-1) * LGdof + idof
                ;tdoc = (vvol-1) * LGdof ! shorthand;
                ;idoc = 0

                if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                    ;hessian(tdoc+idof                  ,tdof) = one ! diagonal elements;
                else
                    ;hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,0,vvol+1) - dFFdRZ(idoc+1:idoc+LGdof,1,idof,1,vvol+0)
                    if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                        hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                                    + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,0) &
                                                                    + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,0) &
                                                                    - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,1) * dmupfdx(vvol+0,1,1,idof,1) &
                                                                    - dBBdmp(idoc+1:idoc+LGdof,vvol+0,1,2) * dmupfdx(vvol+0,1,2,idof,1)
                    endif ! end of if( Lconstraint.eq.1 );
                endif ! end of if( Lextrap.eq.1 .and. vvol.eq.1 )


                if( vvol.lt.Mvol-1 ) then
                  tdof = (vvol+0) * LGdof + idof
                  tdoc = (vvol-1) * LGdof ! shorthand;
                  idoc = 0
                  if( Lextrap.eq.1 .and. vvol.eq.1 ) then
                    if    ( im(idof).le.0                     ) then ; hessian(tdoc+idof,tdof) = - one
                    else                                             ; hessian(tdoc+idof,tdof) = - one
                    endif
                  else
                    hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
                    if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                      hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof) =  hessian(tdoc+idoc+1:tdoc+idoc+LGdof,tdof)                       &
                                                                  + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,1) &
                                                                  + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,1)
                    endif ! end of if( Lconstraint.eq.1 ) then;
                  endif
                endif ! end of if( vvol.lt.Mvol-1 ) ;


                if( vvol.eq.Mvol-1 ) then
                  tdoc = (vvol-1) * LGdof ! shorthand ;
                  idoc = 0
                  dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) = dFFdRZ(idoc+1:idoc+LGdof,0,idof,1,vvol+1)
                  if( Lconstraint.eq.1 ) then ! this is a little clumsy;
                    dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof) =  dessian(tdoc+idoc+1:tdoc+idoc+LGdof,idof)                       &
                                                                + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,1) * dmupfdx(vvol+1,1,1,idof,1) &
                                                                + dBBdmp(idoc+1:idoc+LGdof,vvol+1,0,2) * dmupfdx(vvol+1,1,2,idof,1)
                  endif ! end of if( Lconstraint.eq.1 ) then;

                endif ! end of if( vvol.lt.Mvol-1 ) ;

              else ! Global constraint


                do ivol = 1, Mvol-1
                  tdoc = (ivol-1) * LGdof ! shorthand ;
                  tdof = (vvol-1) * LGdof + idof

                  if( ivol.eq.vvol-1 ) then
                    hessian(tdoc+1:tdoc+LGdof,tdof) =  dFFdRZ(1:LGdof,0,idof,1,ivol+1)
                  elseif( ivol.eq.vvol ) then 
                    hessian(tdoc+1:tdoc+LGdof,tdof) =  dFFdRZ(1:LGdof,0,idof,0,ivol+1) - dFFdRZ(1:LGdof,1,idof,1,ivol)
                  elseif( ivol.eq.vvol+1 ) then
                    hessian(tdoc+1:tdoc+LGdof,tdof) =                                  - dFFdRZ(1:LGdof,1,idof,0,ivol)
                  endif




                  hessian(tdoc+1:tdoc+LGdof,tdof) = hessian(tdoc+1:tdoc+LGdof,tdof)                              &
                                                    + dBBdmp(1:LGdof,ivol+1,0,1) * dmupfdx(ivol+1,vvol,1,idof,1) &
                                                    + dBBdmp(1:LGdof,ivol+1,0,2) * dmupfdx(ivol+1,vvol,2,idof,1) &
                                                    - dBBdmp(1:LGdof,ivol+0,1,1) * dmupfdx(ivol+0,vvol,1,idof,1) &
                                                    - dBBdmp(1:LGdof,ivol+0,1,2) * dmupfdx(ivol+0,vvol,2,idof,1)

                enddo
              endif ! matches if( LocalConstraint );
            enddo ! matches do issym ;
          enddo ! matches do irz ;
        enddo ! matches do ii ;

      else ! matches if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;

        FATAL( dforce, .true., need to provide suitable values for hessian in case of field failure )

      endif ! end of if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) ;
    enddo ! end of do vvol;


  endif ! end of if( LcomputeDerivatives ) ;





  RETURN(dforce)

end subroutine dforce

