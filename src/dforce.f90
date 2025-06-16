
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



  

  if( LocalConstraint ) then
    SALLOCATE( dmupfdx, (1:Mvol,    1:1,1:2,1:LGdof,0:1), zero )
  else
    SALLOCATE( dmupfdx, (1:Mvol, 1:Mvol-1,1:2,1:LGdof,1), zero ) ! TODO change the format to put vvol in last index position...
  endif

  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;


  WCALL( dforce, packxi,( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                          iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack, LcomputeDerivatives, LComputeAxis ) )



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




end subroutine dforce

