
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
    use cputiming, only: Txspech

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    character :: ldate*8, ltime*10, arg*100

    call MPI_INIT(ierr)

    call set_mpi_comm(MPI_COMM_WORLD)

    cpus = MPI_WTIME()
    cpuo = cpus

    skip_write = .false.

    cput = MPI_WTIME()
    if (myid == 0) then

        write (ounit, '("xspech : ", 10x ," : version = "F5.2)') version
        call date_and_time(ldate, ltime)
        write (ounit, '("xspech : ", 10x ," : ")')
        write (ounit, 1000) cput - cpus, ldate(1:4), ldate(5:6), ldate(7:8), ltime(1:2), ltime(3:4), ltime(5:6), machprec, vsmall, small

        write (ounit, '("xspech : ", 10x ," : ")')
        write (ounit, '("xspech : ",f10.2," : parallelism : ncpu=",i3," ; nthreads=",i3," ;")') cput - cpus, ncpu, 1

        call read_command_args()

        call initialize_inputs()

        write (ounit, '("xspech : ", 10x ," : ")')
        write (ounit, '("xspech : ",f10.2," : begin execution ; calling global:readin ;")') cput - cpus

        call read_inputlists_from_file()

        call check_inputs()

    end if

    call broadcast_inputs()

    call preset()

    if (myid == 0) then
        call wrtend()
    end if

    call spec()

    call final_diagnostics()

    if (myid == 0) then

        call wrtend()
    end if

    call MPI_Barrier(MPI_COMM_SPEC, ierr)

    if (myid == 0) then
        cput = MPI_WTIME()
        write (ounit, '("xspech : ", 10x ," :")')
        write (ounit, '("xspech : ",f10.2," : myid=",i3," : time="f8.2"m = "f6.2"h = "f5.2"d ;")') cput - cpus, myid, (cput - cpus)/(/60, 60*60, 24*60*60/)
    end if

    call MPI_FINALIZE(ierr)

    stop

1000 format("xspech : ", f10.2, " : date="a4"/"a2"/"a2" , "a2":"a2":"a2" ; machine precision="es9.2" ; vsmall="es9.2" ; small="es9.2" ;")

end subroutine xspech

subroutine read_command_args

    use fileunits, only: ounit
    use inputlist, only: Wreadin
    use allglobal, only: cpus, myid, ext, MPI_COMM_SPEC, write_spec_namelist

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical :: Lspexist
    integer :: iargc, iarg, numargs, extlen, sppos

    character(len=100) :: arg

    if (myid == 0) then

        cput = MPI_WTIME()

        call getarg(1, arg)

        write (ounit, '("rdcmdl : ", 10x ," : ")')
        select case (trim(arg))
        case ("", "-h", "--help")
            write (ounit, '("rdcmdl : ", 10x ," : file extension must be given as first command line argument ;")')
            write (ounit, '("rdcmdl : ", 10x ," : Usage: <mpiexec> xspec input_file <arguments>")')
            write (ounit, '("rdcmdl : ", 10x ," : Other options:")')
            write (ounit, '("rdcmdl : ", 10x ," :     -h / --help :  print help information.")')
            write (ounit, '("rdcmdl : ", 10x ," :     -i / --init :  generate a template input file.")')
            write (ounit, '("rdcmdl : ", 10x ," : Additional arguments:")')
            write (ounit, '("rdcmdl : ", 10x ," :     -readin : print debugging information during reading inputs")')
            call MPI_ABORT(MPI_COMM_SPEC, 0, ierr)
        case ("-i", "--init")
            write (ounit, '("rdcmdl : ", 10x ," : write a template input file in example.sp")')
            call write_spec_namelist()
            call MPI_ABORT(MPI_COMM_SPEC, 0, ierr)
        case default
            extlen = len_trim(arg)
            sppos = index(arg, ".sp", .true.)
            if (sppos == extlen - 2) then
                arg = arg(1:extlen - 3)
            end if
            ext = trim(arg)

            write (ounit, '("rdcmdl : ", 10x ," : ")')
            write (ounit, '("rdcmdl : ",f10.2," : ext = ",a100)') cput - cpus, ext
        end select

        numargs = iargc()

        if (numargs > 1) then
            iarg = 1
            do while (iarg < numargs)
                iarg = iarg + 1; call getarg(iarg, arg)
                select case (arg)
                case ("-readin"); Wreadin = .true.
                case ("-p4pg"); iarg = iarg + 1; call getarg(iarg, arg)
                case ("-p4wd"); iarg = iarg + 1; call getarg(iarg, arg)
                case default; write (ounit, '("rdcmdl : ",f10.2," : myid=",i3," : argument not recognized ; arg = ",a100)') cput - cpus, myid, arg
                end select
            end do
        end if

    end if

end subroutine read_command_args

subroutine spec

    use constants, only: zero, one, pi2, mu0

    use numerical, only: vsmall, logtolerance

    use fileunits, only: ounit, lunit

    use inputlist, only: Wmacros, Wxspech, &
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

    use cputiming, only: Txspech

    use allglobal, only: wrtend, ncpu, myid, cpus, ext, &
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
                         ForceErr, BnsErr, &
                         efmn, ofmn, cfmn, sfmn, &
                         iBns, iBnc, iVns, iVnc, &
                         Ate, Aze, Ato, Azo, &
                         nfreeboundaryiterations, &
                         beltramierror, &
                         first_free_bound, &
                         dMA, dMB, dMD, dMG, MBpsi, solution, IPDt, &
                         version, &
                         MPI_COMM_SPEC, &
                         force_final, LocalConstraint, dBBdmp, dFFdRZ, dmupfdx, &
                         dRodR, dRodZ, dZodR, dZodZ

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical :: LComputeDerivatives, LContinueFreeboundaryIterations, exist, LupdateBn, LComputeAxis
    integer :: imn, lmn, lNfp, lim, lin, ii, ideriv, stat
    integer :: vvol, ifail, wflag, iflag, vflag
    real(8) :: rflag, lastcpu, lRwc, lRws, lZwc, lZws, lItor, lGpol, lgBc, lgBs
    real(8), allocatable :: position(:), gradient(:)
    character :: pack
    integer :: Lfindzero_old, mfreeits_old
    real(8) :: gBnbld_old
    integer :: lnPtrj, numTrajTotal

    cpuo = MPI_WTIME()

    if (NGdof < 0) then
        write (6, '("xspech :      fatal : myid=",i3," ; NGdof.lt.0 ; counting error;")') myid
        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
        stop "xspech : NGdof.lt.0 : counting error ;"
    end if

    if (allocated(position)) deallocate (position)
    allocate (position(0:NGdof), stat=astat)
    position(0:NGdof) = zero

    nfreeboundaryiterations = -1

9000 nfreeboundaryiterations = nfreeboundaryiterations + 1

    if (NGdof > 0) then

        pack = 'P'
        LComputeAxis = .true.
        call packxi(NGdof, position(0:NGdof), Mvol, mn, iRbc(1:mn, 0:Mvol), iZbs(1:mn, 0:Mvol), &
                    iRbs(1:mn, 0:Mvol), iZbc(1:mn, 0:Mvol), pack, .false., LComputeAxis)

    end if

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if
        vflag = 0
        call volume(vvol, vflag)

        if (Ladiabatic == 0) adiabatic(vvol) = pressure(vvol)*vvolume(vvol)**gamma

    end do

    if (Mvol > Nvol) then; adiabatic(Mvol) = zero; pressure(Mvol) = zero
    end if

    if (Wxspech .and. myid == 0) then
        cput = MPI_WTIME()
        write (ounit, '("xspech : ",f10.2," : myid=",i3," ; adiabatic constants = "999es13.5)') cput - cpus, myid, adiabatic(1:Mvol)
    end if

    pressure(1:Mvol) = adiabatic(1:Mvol)/vvolume(1:Mvol)**gamma

    lastcpu = MPI_WTIME()

    if (allocated(force_final)) deallocate (force_final)
    allocate (force_final(0:NGdof), stat=astat)
    force_final(0:NGdof) = zero

    LComputeDerivatives = .false.
    LComputeAxis = .true.

    call dforce(NGdof, position(0:NGdof), force_final(0:NGdof), LComputeDerivatives, LComputeAxis)

1000 format("xspech : ", f10.2, " : #freeits=", i3, " ; ":"|f|="es12.5" ; ":"time=", f10.2, "s ;":" log"a5, :"="28f6.2" ...")
1001 format("xspech : ", 10x, " :          ", 3x, " ; ":"    "12x "   ":"     ", 10x, "  ;":" log"a5, :"="28f6.2" ...")

    select case (Igeometry)
    case (1); tflux(1) = dtflux(1); pflux(1) = dpflux(1)
    case (2:3); tflux(1) = dtflux(1); pflux(1) = zero
    end select

    do vvol = 2, Mvol; tflux(vvol) = tflux(vvol - 1) + dtflux(vvol)
        ; pflux(vvol) = pflux(vvol - 1) + dpflux(vvol)
    end do

    tflux(1:Mvol) = tflux(1:Mvol)*pi2/phiedge
    pflux(1:Mvol) = pflux(1:Mvol)*pi2/phiedge

    call ra00aa('W')

    if (myid == 0) then
        call wrtend()
    end if

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

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer :: iocons, llmodnp, vvol, iflag, cpu_id
    real :: sumI
    real(8), allocatable :: Bt00(:, :, :)
    real(8) :: work(0:1, -1:2)

    if (Ltransform) then

        if (Lsvdiota /= 1) then
            write (6, '("xspech :      fatal : myid=",i3," ; Lsvdiota.ne.1 ; Lsvdiota needs to be one for s.f.l transformation;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "xspech : Lsvdiota.ne.1 : Lsvdiota needs to be one for s.f.l transformation ;"
        end if

        do vvol = 1, Mvol
            call brcast(vvol)
        end do

        iflag = -1
        do vvol = 1, Mvol
            call IsMyVolume(vvol)
            if (IsMyVolumeValue == 0) then
                cycle
            elseif (IsMyVolumeValue == -1) then
                if (.true.) then
                    write (6, '("xspech :      fatal : myid=",i3," ; .true. ; Unassociated volume;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "xspech : .true. : Unassociated volume ;"
                end if
            end if

            if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
            else; Lcoordinatesingularity = .true.
            end if

            call tr00ab(vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, vvol))
        end do

        do vvol = 1, Mvol
            call WhichCpuID(vvol, cpu_id)
            call MPI_BCAST(diotadxup(0:1, -1:2, vvol), 8, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
            call MPI_BCAST(dlambdaout(1:lmns, vvol, 0:1), 2*lmns, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr)
        end do

    end if

    if (allocated(Bt00)) deallocate (Bt00)
    allocate (Bt00(1:Mvol, 0:1, -1:2), stat=astat)
    Bt00(1:Mvol, 0:1, -1:2) = zero

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if

        do iocons = 0, 1
            if (Lcoordinatesingularity .and. iocons == 0) cycle
            if (vvol == Nvol + 1 .and. iocons == 1) cycle
            call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), 0, iocons)

            Btemn(1:mn, iocons, vvol) = efmn(1:mn)
            Btomn(1:mn, iocons, vvol) = ofmn(1:mn)
            Bzemn(1:mn, iocons, vvol) = cfmn(1:mn)
            Bzomn(1:mn, iocons, vvol) = sfmn(1:mn)
        end do
    end do

    do vvol = 1, Mvol - 1
        IPDt(vvol) = pi2*(Bt00(vvol + 1, 0, 0) - Bt00(vvol, 1, 0))
    end do

    deallocate (Bt00, stat=astat)

    sumI = 0
    do vvol = 1, Mvol
        Ivolume(vvol) = mu(vvol)*dtflux(vvol)*pi2 + sumI
        sumI = Ivolume(vvol)
    end do

    if (myid == 0) then
        cput = MPI_WTIME()

        if (nPpts > 0) then
            write (ounit, '("xspech : ", 10x ," :")')
            write (ounit, '("xspech : ",f10.2," : myid=",i3," ; Poincare plot ; odetol="es8.1" ; nPpts="i7" ;":" nPtrj="24(i5",")" ...")') &
                cput - cpus, myid, odetol, nPpts, nPtrj(1:min(Mvol, 24))
        end if

        if (Lcheck == 1) then
            write (ounit, '("xspech : ", 10x ," :")')
            write (ounit, '("xspech : ",f10.2," : myid=",i3," ; calling jo00aa; computing error in field ;")') cput - cpus, myid
        end if
    end if

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if

        if (myid == modulo(vvol - 1, ncpu) .and. myid < Mvol) then

            if (.not. ImagneticOK(vvol)) then; cput = MPI_WTIME(); write (ounit, 1002) cput - cpus; write (ounit, 1002) cput - cpus, myid, vvol, ImagneticOK(vvol); cycle
            end if

            if (Lcheck == 1) then
                call jo00aa(vvol, Ntz, Iquad(vvol), mn)
            end if

        end if
    end do

1002 format("xspech : ", f10.2, " :":" myid=", i3, " ; vvol=", i3, " ; IBeltrami="L2" ; construction of Beltrami field failed ;")

    do vvol = 1, Mvol; llmodnp = modulo(vvol - 1, ncpu)

        call MPI_BCAST(Btemn(1:mn, 0:1, vvol), mn*2, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Bzemn(1:mn, 0:1, vvol), mn*2, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Btomn(1:mn, 0:1, vvol), mn*2, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Bzomn(1:mn, 0:1, vvol), mn*2, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(beltramierror(vvol, 1:9), 9, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

    end do

end subroutine final_diagnostics
