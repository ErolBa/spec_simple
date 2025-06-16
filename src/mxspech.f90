
program spec_main
  implicit none
  call xspech
end program spec_main

subroutine xspech

  use numerical
  use allglobal, only: set_mpi_comm, myid, ncpu, cpus, version, MPI_COMM_SPEC, &
                       wrtend, read_inputlists_from_file, check_inputs, broadcast_inputs, skip_write, &
                       ext
  use inputlist, only: initialize_inputs, Wxspech
  use fileunits, only: ounit
  use sphdf5,    only: init_outfile, &
                       mirror_input_to_outfile, &
                       init_convergence_output, &
                       hdfint, finish_outfile, write_grid
  use cputiming, only: Txspech

  LOCALS

  CHARACTER            :: ldate*8, ltime*10, arg*100

  call MPI_INIT( ierr )

  BEGIN(xspech)

  call set_mpi_comm(MPI_COMM_WORLD)

  cpus = GETTIME
  cpuo = cpus

  skip_write = .false.

  cput = GETTIME
  if( myid.eq.0 ) then

    write(ounit,'("xspech : ", 10x ," : version = "F5.2)') version
    call date_and_time( ldate, ltime )
    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,1000) cput-cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small

    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : parallelism : ncpu=",i3," ; nthreads=",i3," ;")') cput-cpus, ncpu, nthreads

    call read_command_args()

    call initialize_inputs()

    write(ounit,'("xspech : ", 10x ," : ")')
    write(ounit,'("xspech : ",f10.2," : begin execution ; calling global:readin ;")') cput-cpus

    call read_inputlists_from_file()

    call check_inputs()

  endif ! myid.eq.0

  call broadcast_inputs()

  call preset()

  call init_outfile()

  call mirror_input_to_outfile()



  if ( myid .eq. 0 ) then ! save restart file;
    call wrtend() ! write initial restart file
  endif

  call init_convergence_output()


  call spec()

  call final_diagnostics()

  call write_grid()

  if( myid.eq.0 ) then

    call wrtend()
  endif

  call hdfint()

  call finish_outfile()

  call ending()

  call MPI_Barrier(MPI_COMM_SPEC, ierr)

  if (myid.eq.0) then
   cput = GETTIME
   write(ounit,'("xspech : ", 10x ," :")')
   write(ounit,'("xspech : ",f10.2," : myid=",i3," : time="f8.2"m = "f6.2"h = "f5.2"d ;")') cput-cpus, myid, (cput-cpus) / (/ 60, 60*60, 24*60*60 /)
  endif

  MPIFINALIZE

  stop

1000 format("xspech : ",f10.2," : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ;")

end subroutine xspech

subroutine read_command_args

  use fileunits, only: ounit
  use inputlist, only: Wreadin
  use allglobal, only: cpus, myid, ext, MPI_COMM_SPEC, write_spec_namelist

  LOCALS

  LOGICAL              :: Lspexist
  INTEGER              :: iargc, iarg, numargs, extlen, sppos

  CHARACTER(len=100)   :: arg

  if (myid.eq.0) then

    cput = GETTIME

    call getarg( 1, arg )

    write(ounit,'("rdcmdl : ", 10x ," : ")')
    select case (trim(arg))
    case ("", "-h", "--help")
        write(ounit,'("rdcmdl : ", 10x ," : file extension must be given as first command line argument ;")')
        write(ounit,'("rdcmdl : ", 10x ," : Usage: <mpiexec> xspec input_file <arguments>")')
        write(ounit,'("rdcmdl : ", 10x ," : Other options:")')
        write(ounit,'("rdcmdl : ", 10x ," :     -h / --help :  print help information.")')
        write(ounit,'("rdcmdl : ", 10x ," :     -i / --init :  generate a template input file.")')
        write(ounit,'("rdcmdl : ", 10x ," : Additional arguments:")')
        write(ounit,'("rdcmdl : ", 10x ," :     -readin : print debugging information during reading inputs")')
        call MPI_ABORT( MPI_COMM_SPEC, 0, ierr )
    case ("-i", "--init")
        write(ounit,'("rdcmdl : ", 10x ," : write a template input file in example.sp")')
        call write_spec_namelist()
        call MPI_ABORT( MPI_COMM_SPEC, 0, ierr )
    case default
        extlen = len_trim(arg)
        sppos = index(arg, ".sp", .true.) ! search for ".sp" from the back of ext
        if (sppos.eq.extlen-2) then       ! check if ext ends with ".sp"
            arg = arg(1:extlen-3)         ! if this is the case, remove ".sp" from end of ext
        endif
        ext = trim(arg)

        write(ounit,'("rdcmdl : ", 10x ," : ")')
        write(ounit,'("rdcmdl : ",f10.2," : ext = ",a100)') cput-cpus, ext
    end select



    numargs = iargc()

    if( numargs.gt.1 ) then
      iarg = 1
      do while ( iarg < numargs )
        iarg = iarg + 1 ; call getarg( iarg, arg)
        select case( arg )
        case("-readin"   ) ; Wreadin = .true.
        case("-p4pg"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
        case("-p4wd"     ) ; iarg = iarg + 1 ; call getarg( iarg, arg) ! TODO: what is this?
        case default       ; write(ounit,'("rdcmdl : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput-cpus, myid, arg
        end select
      enddo
    endif

  end if ! check for myid.eq.0

end subroutine read_command_args


subroutine spec

  use constants, only : zero, one, pi2, mu0

  use numerical, only : vsmall, logtolerance

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wxspech, &
                        Nfp, Igeometry, Nvol, Lrad, &
                        tflux, pflux, phiedge, pressure, pscale, helicity, Ladiabatic, adiabatic, gamma, &
                        Rbc, Zbs, Rbs, Zbc, &
                        Lconstraint, &
                        Lfreebound, mfreeits, gBntol, gBnbld, vcasingtol, LautoinitBn, &
                        Lfindzero, LautoinitBn, &
                        odetol, nPpts, nPtrj, &
                        LHevalues, LHevectors, LHmatrix, Lperturbed, Lcheck, &
                        Lzerovac, &
                        mu, Isurf, Ivolume

  use cputiming, only : Txspech

  use allglobal, only : wrtend, ncpu, myid, cpus, ext, &
                        Mvol, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, &
                        Ntz, &
                        LGdof, NGdof, &
                        iRbc, iZbs, iRbs, iZbc, &
                        BBe, IIo, BBo, IIe, &
                        vvolume, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        dtflux, dpflux, &
                        ImagneticOK, &
                        ForceErr, BnsErr,&
                        efmn, ofmn, cfmn, sfmn, &
                        iBns, iBnc, iVns, iVnc, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging; 09 Mar 17;
                        nfreeboundaryiterations, &
                        beltramierror, &
                        first_free_bound, &
                        dMA, dMB, dMD, dMG, MBpsi, solution, IPDt, &
                        version, &
                        MPI_COMM_SPEC, &
                        force_final, Lhessianallocated, LocalConstraint, hessian, dBBdmp, dFFdRZ, dmupfdx, &
                        dRodR, dRodZ, dZodR, dZodZ



  LOCALS

  LOGICAL              :: LComputeDerivatives, LContinueFreeboundaryIterations, exist, LupdateBn, LComputeAxis
  INTEGER              :: imn, lmn, lNfp, lim, lin, ii, ideriv, stat
  INTEGER              :: vvol, ifail, wflag, iflag, vflag
  REAL                 :: rflag, lastcpu, lRwc, lRws, lZwc, lZws, lItor, lGpol, lgBc, lgBs
  REAL,    allocatable :: position(:), gradient(:)
  CHARACTER            :: pack
  INTEGER              :: Lfindzero_old, mfreeits_old
  REAL                 :: gBnbld_old
  INTEGER              :: lnPtrj, numTrajTotal

  

  cpuo = GETTIME

  FATAL( xspech, NGdof.lt.0, counting error )

  SALLOCATE( position, (0:NGdof), zero ) ! position ; NGdof = #geometrical degrees-of-freedom was computed in preset;



  nfreeboundaryiterations = -1

9000 nfreeboundaryiterations = nfreeboundaryiterations + 1

  if( NGdof.gt.0 ) then ! pack geometry into vector; 14 Jan 13;

   pack = 'P'
   LComputeAxis = .true.
   WCALL( xspech, packxi, ( NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), &
                            iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), pack, .false., LComputeAxis ) )

  endif




  do vvol = 1, Mvol

   LREGION(vvol)
   vflag = 0
   WCALL( xspech, volume, ( vvol, vflag ) ) ! compute volume;
   
   if( Ladiabatic.eq.0 ) adiabatic(vvol) = pressure(vvol) * vvolume(vvol)**gamma ! initialize adiabatic constants using supplied pressure profile;

  enddo ! end of do vvol = 1, Mvol;

  if( Mvol.gt.Nvol ) then ; adiabatic(Mvol) = zero ; pressure(Mvol) = zero ! these are never used; 15 May 13;
  endif

  if( Wxspech .and. myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("xspech : ",f10.2," : myid=",i3," ; adiabatic constants = "999es13.5)') cput-cpus, myid, adiabatic(1:Mvol)
  endif

  pressure(1:Mvol) = adiabatic(1:Mvol) / vvolume(1:Mvol)**gamma ! this matches construction of adiabatic above;



  lastcpu = GETTIME

  
    SALLOCATE( force_final, (0:NGdof), zero )

    LComputeDerivatives = .false.
    LComputeAxis = .true.

    WCALL( xspech, dforce, ( NGdof, position(0:NGdof), force_final(0:NGdof), LComputeDerivatives, LComputeAxis) )

1000 format("xspech : ",f10.2," : #freeits=",i3," ; ":"|f|="es12.5" ; ":"time=",f10.2,"s ;":" log"a5,:"="28f6.2" ...")
1001 format("xspech : ", 10x ," :          ",3x," ; ":"    "  12x "   ":"     ", 10x ,"  ;":" log"a5,:"="28f6.2" ...")

  select case( Igeometry )                                  ! 08 Feb 16;
  case( 1   ) ; tflux(1) = dtflux(1) ; pflux(1) = dpflux(1) ! 08 Feb 16;
  case( 2:3 ) ; tflux(1) = dtflux(1) ; pflux(1) =   zero    ! 08 Feb 16;
  end select                                                ! 08 Feb 16;

  do vvol = 2, Mvol; tflux(vvol) = tflux(vvol-1) + dtflux(vvol) ! 01 Jul 14;
  ;                  pflux(vvol) = pflux(vvol-1) + dpflux(vvol) ! 01 Jul 14;
  enddo

  tflux(1:Mvol) = tflux(1:Mvol) * pi2 / phiedge ! this is the "inverse" operation defined in preset; 19 Jul 16;
  pflux(1:Mvol) = pflux(1:Mvol) * pi2 / phiedge


  WCALL( xspech, ra00aa, ('W') ) ! this writes vector potential to file;

  if( myid.eq.0 ) then ! write restart file; note that this is inside free-boundary iteration loop; 11 Aug 14;
   WCALL( xspech, wrtend ) ! write restart file; save initial input;
  endif

end subroutine spec

subroutine final_diagnostics

  use inputlist, only: nPtrj, nPpts, Igeometry, Lcheck, Nvol, odetol, &
                       Isurf, Ivolume, mu, Wmacros, Ltransform, Lsvdiota, Lconstraint
  use fileunits, only: ounit
  use constants, only: zero
  use allglobal, only: pi2, myid, ncpu, MPI_COMM_SPEC, cpus, Mvol, Ntz, mn, &
                       beltramierror, Lcoordinatesingularity, &
                       Lplasmaregion, Lvacuumregion, &
                       Btemn, Bzemn, Btomn, Bzomn, &
                       efmn, ofmn, cfmn, sfmn, &
                       IPDt, ImagneticOK, dtflux, Iquad, lmns, Nt, Nz, diotadxup, &
                       IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                       dlambdaout, diotadxup


  LOCALS

  integer              :: iocons, llmodnp, vvol, iflag, cpu_id
  real                 :: sumI
  REAL,    allocatable :: Bt00(:,:,:)
  REAL                 :: work(0:1,-1:2) 


if( Ltransform ) then

  FATAL(xspech, Lsvdiota.ne.1, Lsvdiota needs to be one for s.f.l transformation)

  do vvol=1,Mvol
    call brcast(vvol)
  enddo

  iflag = -1
  do vvol = 1, Mvol
    call IsMyVolume(vvol)
    if (IsMyVolumeValue.eq.0) then
      cycle
    elseif (IsMyVolumeValue.eq.-1) then
      FATAL( xspech, .true., Unassociated volume )
    endif

    LREGION( vvol )

    call tr00ab( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2, vvol) ) ! stores lambda in a global variable.
  enddo

  do vvol = 1, Mvol
    call WhichCpuID( vvol, cpu_id )
    RlBCAST( diotadxup(0:1,-1:2,vvol), 8, cpu_id  )
    RlBCAST( dlambdaout(1:lmns, vvol, 0:1), 2*lmns, cpu_id  )
  enddo

endif


  

  SALLOCATE( Bt00, (1:Mvol, 0:1, -1:2) , zero)

  do vvol = 1, Mvol

    LREGION(vvol)

    do iocons = 0, 1
	  if( Lcoordinatesingularity .and. iocons.eq.0 ) cycle
          if( vvol.eq.Nvol+1 .and. iocons.eq.1 ) cycle
      call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), 0, iocons)

      Btemn(1:mn, iocons, vvol) = efmn(1:mn)
      Btomn(1:mn, iocons, vvol) = ofmn(1:mn)
      Bzemn(1:mn, iocons, vvol) = cfmn(1:mn)
      Bzomn(1:mn, iocons, vvol) = sfmn(1:mn)
    enddo
  enddo

  do vvol = 1, Mvol-1
    IPDt(vvol) = pi2 * (Bt00(vvol+1, 0, 0) - Bt00(vvol, 1, 0))
  enddo

  DALLOCATE( Bt00 )

  sumI = 0
  do vvol = 1, Mvol
    Ivolume(vvol) = mu(vvol) * dtflux(vvol) * pi2 + sumI    ! factor pi2 due to normalization in preset
    sumI = Ivolume(vvol)                                    ! Sum over all volumes since this is how Ivolume is defined
  enddo

  if (myid.eq.0) then
   cput = GETTIME

   if( nPpts.gt.0 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; Poincare plot ; odetol="es8.1" ; nPpts="i7" ;":" nPtrj="24(i5",")" ...")') &
                 cput-cpus, myid, odetol, nPpts, nPtrj(1:min(Mvol,24))
   endif

   if( Lcheck.eq.1 ) then
    write(ounit,'("xspech : ", 10x ," :")')
    write(ounit,'("xspech : ",f10.2," : myid=",i3," ; calling jo00aa; computing error in field ;")') cput-cpus, myid
   endif
  endif

  do vvol = 1, Mvol

   LREGION(vvol)

   if( myid.eq.modulo(vvol-1,ncpu) .and. myid.lt.Mvol) then ! the following is in parallel; 20 Jun 14;

    if( .not.ImagneticOK(vvol) ) then ; cput = GETTIME ; write(ounit,1002) cput-cpus ; write(ounit,1002) cput-cpus, myid, vvol, ImagneticOK(vvol) ; cycle
    endif


    if( Lcheck.eq.1 ) then
     call jo00aa( vvol, Ntz, Iquad(vvol), mn )
    endif

   endif ! myid.eq.modulo(vvol-1,ncpu)
  enddo ! end of do vvol = 1, Mvol; ! end of parallel diagnostics loop; 03 Apr 13;


1002 format("xspech : ",f10.2," :":" myid=",i3," ; vvol=",i3," ; IBeltrami="L2" ; construction of Beltrami field failed ;")



  do vvol = 1, Mvol ; llmodnp = modulo(vvol-1,ncpu)


   RlBCAST( Btemn(1:mn,0:1,vvol), mn*2, llmodnp ) ! this is computed in lbpol; 07 Dec 16;
   RlBCAST( Bzemn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Btomn(1:mn,0:1,vvol), mn*2, llmodnp )
   RlBCAST( Bzomn(1:mn,0:1,vvol), mn*2, llmodnp )

   RlBCAST( beltramierror(vvol,1:9), 9, llmodnp ) ! this is computed in jo00aa; 21 Aug 18;

  enddo ! end of do vvol = 1, Mvol; 01 Jul 14;

end subroutine final_diagnostics

subroutine ending

  use constants, only : zero

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wxspech, Ltiming

  use cputiming

  use allglobal, only : myid, cpus, mn, MPI_COMM_SPEC, ext


  LOCALS

  REAL      :: Ttotal, dcpu, ecpu
  CHARACTER :: date*8, time*10



  cpui = GETTIME ; cpuo = cpui ! see macro expansion for begin; 11 Aug 14;

  cput = GETTIME


  cput = GETTIME ; dcpu = cput-cpus

  if( Ltiming .and. myid.eq.0 ) then

   Ttotal = zero

   write(ounit,'("ending : ",f10.2," : time spent in wrtend =",f10.2," ;")') dcpu, Twrtend ; Ttotal = Ttotal + Twrtend
   write(ounit,'("ending : ",f10.2," : time spent in readin =",f10.2," ;")') dcpu, Treadin ; Ttotal = Ttotal + Treadin

   ecpu = Ttotal-dcpu ! error in actual cpu time and calculated cpu time;  7 Mar 13;

   write(ounit,'("ending : ",f10.2," : Ttotal =",f10.2," s = "f8.2" m = "f6.2" h ; Timing Error = ",f10.2,"s = ",f10.2,"%")') &
dcpu, Ttotal / (/ 1, 60, 3600 /), ecpu, 100*ecpu/dcpu

  endif ! end of if( Ltiming .and. myid.eq.0 ) then; 01 Jul 14;

  if( myid.eq.0 ) then
   call date_and_time(date,time)
   write(ounit,'("ending : ", 10x ," : ")')
   write(ounit,1000) dcpu, myid, dcpu / (/ 1, 60, 60*60, 24*60*60 /), date(1:4), date(5:6), date(7:8), time(1:2), time(3:4), time(5:6), ext
   write(ounit,'("ending : ", 10x ," : ")')
  endif ! end of if( myid.eq.0 ) ; 14 Jan 15;

1000 format("ending : ",f10.2," : myid=",i3," ; completion ; time=",f10.2,"s = "f8.2"m = "f6.2"h = "f5.2"d ; date= "&
  a4"/"a2"/"a2" ; time= "a2":"a2":"a2" ; ext = "a60)

end subroutine ending

