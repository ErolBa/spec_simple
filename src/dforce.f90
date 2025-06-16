
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



  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  real(8)      :: cput, cpui, cpuo=0

  INTEGER, parameter   :: NB = 3 ! optimal workspace block size for LAPACK:DSYSVX;

  INTEGER, intent(in)  :: NGdof               ! dimensions;
  real(8),    intent(in)  :: position(0:NGdof)
  real(8),    intent(out) :: force(0:NGdof)      ! force;
  LOGICAL, intent(in)  :: LComputeDerivatives !

  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, tdoc, tdoc_ntz, idoc, idof, tdof, jdof, ivol, imn, ll, ihybrd1, lwa, Ndofgl, llmodnp
  INTEGER              :: maxfev, ml, muhybr, mode, nprint, nfev, ldfjac, lr, Nbc, NN, cpu_id, ideriv
  real(8)                 :: epsfcn, factor
  real(8)                 :: Fdof(1:Mvol-1), Xdof(1:Mvol-1)
  INTEGER              :: ipiv(1:Mvol)
  real(8), allocatable    :: fjac(:, :), r(:), Fvec(:), dpfluxout(:)

  INTEGER              :: status(MPI_STATUS_SIZE), request_recv, request_send, cpu_send
  INTEGER              :: id
  INTEGER              :: iflag, idgesv, Lwork
  INTEGER              :: idofr,idofz,tdofr,tdofz

  CHARACTER            :: packorunpack
  EXTERNAL             :: dfp100, dfp200

  LOGICAL              :: LComputeAxis, dfp100_logical



  

  if( LocalConstraint ) then
if( allocated( dmupfdx ) ) deallocate( dmupfdx )
allocate( dmupfdx(1:Mvol,1:1,1:2,1:LGdof,0:1), stat=astat )
dmupfdx(1:Mvol,1:1,1:2,1:LGdof,0:1) = zero
  else
if( allocated( dmupfdx ) ) deallocate( dmupfdx )
allocate( dmupfdx(1:Mvol,1:Mvol-1,1:2,1:LGdof,1), stat=astat )
dmupfdx(1:Mvol,1:Mvol-1,1:2,1:LGdof,1) = zero
  endif

  packorunpack = 'U' ! unpack geometrical degrees-of-freedom;


  call packxi( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                          iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), packorunpack, LcomputeDerivatives, LComputeAxis )



  Xdof(1:Mvol-1) = zero;

  if( LocalConstraint ) then

if( allocated( Fvec ) ) deallocate( Fvec )
allocate( Fvec(1:Mvol-1), stat=astat )
Fvec(1:Mvol-1) = zero

    Ndofgl = 0; Fvec(1:Mvol-1) = 0; dfp100_logical = .FALSE.;
    Xdof(1:Mvol-1) = dpflux(2:Mvol) + xoffset

    dBdX%L = LComputeDerivatives
    call dfp100(Ndofgl, Xdof, Fvec, dfp100_logical )

    deallocate(Fvec,stat=astat)

  else

    IPDtdPf = zero
    Xdof(1:Mvol-1)   = dpflux(2:Mvol) + xoffset

      Ndofgl = Mvol-1


if( allocated( Fvec ) ) deallocate( Fvec )
allocate( Fvec(1:Ndofgl), stat=astat )
Fvec(1:Ndofgl) = zero

    dfp100_logical = .FALSE.

    call dfp100(Ndofgl, Xdof(1:Mvol-1), Fvec(1:Ndofgl), dfp100_logical)

if( allocated( dpfluxout ) ) deallocate( dpfluxout )
allocate( dpfluxout(1:Ndofgl), stat=astat )
dpfluxout(1:Ndofgl) = zero
    if ( myid .eq. 0 ) then

        dpfluxout = Fvec
        call DGESV( Ndofgl, 1, IPdtdPf, Ndofgl, ipiv, dpfluxout, Ndofgl, idgesv )

        dpflux(2:Mvol) = dpflux(2:Mvol) - dpfluxout(1:Mvol-1)
    endif

call MPI_BCAST(dpfluxout(1:Ndofgl), Ndofgl, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
call MPI_BCAST(dpflux(1:Mvol), Mvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
    do vvol = 2, Mvol

      call IsMyVolume(vvol)

      if( IsMyVolumeValue .EQ. 0 ) then
          cycle
      else if( IsMyVolumeValue .EQ. -1) then
if( .true. ) then
     write(6,'("dforce :      fatal : myid=",i3," ; .true. ; Unassociated volume;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "dforce : .true. : Unassociated volume ;"
   endif
      endif

      NN = NAdof(vvol)

if( allocated( solution ) ) deallocate( solution )
allocate( solution(1:NN,0:2), stat=astat )
solution(1:NN,0:2) = zero

      packorunpack = 'P'
      call packab( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ! packing;
      call packab( packorunpack, vvol, NN, solution(1:NN,2), 2 ) ! packing;

      solution(1:NN, 0) = solution(1:NN, 0) - dpfluxout(vvol-1) * solution(1:NN, 2)


      packorunpack = 'U'
      call packab( packorunpack, vvol, NN, solution(1:NN,0), 0  ) ! unpacking;

      deallocate(solution ,stat=astat)

    enddo ! end of do vvol = 1, Mvol

    deallocate(Fvec,stat=astat)
    deallocate(dpfluxout,stat=astat)


  endif !matches if( LocalConstraint )


  do vvol = 1, Mvol
    call WhichCpuID(vvol, cpu_id)

    call MPI_BCAST(ImagneticOK(vvol),1,MPI_LOGICAL,cpu_id,MPI_COMM_SPEC,ierr)

    do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
      do ii = 1, mn
call MPI_BCAST(Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
      enddo
    enddo


    if( NOTstellsym ) then
      do ideriv=0,2
      if( (.not.LcomputeDerivatives) .and. (ideriv.ne.0) ) cycle
        do ii = 1, mn
call MPI_BCAST(Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
        enddo
      enddo
    endif
  enddo



  call dfp200( LcomputeDerivatives, vvol )



  do vvol = 1, Mvol

    if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false.
  else                                   ; Lcoordinatesingularity = .true.
  endif
    call brcast( vvol  )

  enddo






  lBBintegral(1:Nvol) = lBBintegral(1:Nvol) * half

  Energy = sum( lBBintegral(1:Nvol) ) ! should also compute beta;




  ;   force(0:NGdof) = zero

  do vvol = 1, Mvol-1

    if( Igeometry.eq.1 .or. vvol.gt.1 ) then ; Lcoordinatesingularity = .false.
  else                                   ; Lcoordinatesingularity = .true.
  endif

    tdoc = (vvol-1) * LGdof

    if( ImagneticOK(vvol) .and. ImagneticOK(vvol+1) ) then ! the magnetic fields in the volumes adjacent to this interface are valid;

      ;  idoc = 0           ! degree-of-constraint counter; set;

      if( Lextrap.eq.1 .and. vvol.eq.1 ) then ! to be made redundant;
if( 2.gt.Mvol ) then
     write(6,'("dforce :      fatal : myid=",i3," ; 2.gt.Mvol ; psifactor needs attention;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "dforce : 2.gt.Mvol : psifactor needs attention ;"
   endif
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

