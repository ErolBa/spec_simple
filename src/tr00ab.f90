
subroutine tr00ab(lvol, mn, NN, Nt, Nz, iflag, ldiota)

    use constants, only: zero, third, half, one, two, pi2, goldenmean

    use numerical, only: vsmall, small, machprec, sqrtmachprec

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wtr00ab, Nvol, Lrad, Mpol, Ntor, &
                         Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol

    use cputiming, only: Ttr00ab

    use allglobal, only: ncpu, cpus, myid, MPI_COMM_SPEC, &
                         pi2nfp, &
                         Mvol, im, in, mns, ims, ins, &
                         YESstellsym, NOTstellsym, &
                         glambda, &
                         Ntz, hNt, hNz, &
                         iotakkii, iotaksub, iotakadd, iotaksgn, &
                         Ate, Aze, Ato, Azo, TT, RTT, &
                         Lcoordinatesingularity, Lvacuumregion, regumm, dlambdaout

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, mn, NN, Nt, Nz, iflag

    real(8), intent(inout) :: ldiota(0:1, -1:2)

    integer :: innout, ll, ii, jj, kk, jb, kb, mj, nj, ideriv, jderiv, id, MM, ielement, nelements, Lcurvature, idof, icon, mi, ni, imupf

    real(8) :: lcpu, mfactor, lss, Dteta, Dzeta, rfac, tol, rnorm, omega, diotaerror

    real(8) :: lAte(0:mn, -1:2), lAze(0:mn, -1:2), lAto(0:mn, -1:2), lAzo(0:mn, -1:2)

    integer :: IA, if04aaf, idgesvx, ipiv(1:NN), iwork4(1:NN)
    real(8), allocatable :: dmatrix(:, :, :), omatrix(:, :), FAA(:, :)
    real(8) :: drhs(1:NN, -1:2), dlambda(1:NN, -1:2)
    real(8) :: Rdgesvx(1:NN), Cdgesvx(1:NN), work4(1:4*NN), rcond, ferr(1), berr(1), ferr2(1:2), berr2(1:2)
    character :: equed

    integer :: maxitn, reqdits, extralength, lrwork, integerwork(1:2*Nt*Nz + 2 + 1), if11def, if11zaf, if11xaf
    integer :: IAA, if04atf, if04arf
    integer :: Ndof, label(-3:Nt + 2, -3:Nz + 2), isym

    integer :: idgelsd, Lwork, Liwork, Irank, nlvl
    real(8) :: sval(1:NN)
    real(8), allocatable :: work(:)

    real(8) :: Bsupt(1:Nt*Nz, -1:2), Bsupz(1:Nt*Nz, -1:2), tdot(1:Nt*Nz)
    real(8) :: Bsubs(1:Nt*Nz, -1:2), Bsubt(1:Nt*Nz, -1:2), Bsubz(1:Nt*Nz, -1:2)

    real(8) :: dotteta, dotzeta

    real(8), allocatable :: rmatrix(:, :, :), rrhs(:, :), rlambda(:, :), wks1(:), wks2(:), AA(:, :)

    integer :: inz(-1:2), lnz
    integer, allocatable :: irow(:, :), jcol(:, :), istr(:), iwork(:)
    real(8), allocatable :: smatrix(:, :), srhs(:, :), slambda(:, :), swork(:)
    character :: duplicate*1, zeros*1, method*8, precon*1, trans*1, check*1

    do innout = 0, 1

        if (Lcoordinatesingularity .and. innout == 0) cycle
        if (Lvacuumregion .and. innout == 1) cycle

        lAte(0:mn, -1:2) = zero
        lAze(0:mn, -1:2) = zero
        lAto(0:mn, -1:2) = zero
        lAzo(0:mn, -1:2) = zero

        do ideriv = -1, 2; id = ideriv

            if (iflag == 1 .and. ideriv /= 0) cycle
            if (iflag == 2 .and. ideriv < 0) cycle
            if (iflag == -1 .and. ideriv > 0) cycle

            do ii = 1, mn

                mi = im(ii)

                if (Lcoordinatesingularity) then
                    do ll = mi, Lrad(lvol), 2

                        ; lAte(ii, id) = lAte(ii, id) + Ate(lvol, id, ii)%s(ll)*RTT(ll, mi, innout, 1)*half
                        ; lAze(ii, id) = lAze(ii, id) - Aze(lvol, id, ii)%s(ll)*RTT(ll, mi, innout, 1)*half
                        if (NOTstellsym) then
                            lAto(ii, id) = lAto(ii, id) + Ato(lvol, id, ii)%s(ll)*RTT(ll, mi, innout, 1)*half
                            lAzo(ii, id) = lAzo(ii, id) - Azo(lvol, id, ii)%s(ll)*RTT(ll, mi, innout, 1)*half
                        end if

                    end do
                else
                    do ll = 0, Lrad(lvol)

                        ; lAte(ii, id) = lAte(ii, id) + Ate(lvol, id, ii)%s(ll)*TT(ll, innout, 1)
                        ; lAze(ii, id) = lAze(ii, id) - Aze(lvol, id, ii)%s(ll)*TT(ll, innout, 1)
                        if (NOTstellsym) then
                            lAto(ii, id) = lAto(ii, id) + Ato(lvol, id, ii)%s(ll)*TT(ll, innout, 1)
                            lAzo(ii, id) = lAzo(ii, id) - Azo(lvol, id, ii)%s(ll)*TT(ll, innout, 1)
                        end if

                    end do
                end if
            end do

            if (Lsparse > 0) then
                call invfft(mn, im(1:mn), in(1:mn), lAte(1:mn, id), lAto(1:mn, id), lAze(1:mn, id), lAzo(1:mn, id), &
                            Nt, Nz, Bsupz(1:Ntz, id), Bsupt(1:Ntz, id))
            end if

        end do

        if (Lsparse == 0 .or. Lsparse == 3) then
            if (allocated(dmatrix)) deallocate (dmatrix)
            allocate (dmatrix(1:NN, 1:NN, -1:2), stat=astat)
            dmatrix(1:NN, 1:NN, -1:2) = zero
            if (allocated(omatrix)) deallocate (omatrix)
            allocate (omatrix(1:NN, 1:NN), stat=astat)
            omatrix(1:NN, 1:NN) = zero
            if (allocated(FAA)) deallocate (FAA)
            allocate (FAA(1:NN, 1:NN), stat=astat)
            FAA(1:NN, 1:NN) = zero
        end if

        if (Lsparse == 0 .or. Lsparse == 3) then

            drhs(1:NN, -1:2) = zero

            dmatrix(1:NN, 1:NN, -1:2) = zero

            do ideriv = -1, 2

                if (iflag == 1 .and. ideriv /= 0) cycle
                if (iflag == 2 .and. ideriv < 0) cycle
                if (iflag == -1 .and. ideriv > 0) cycle
                do kk = 1, mn

                    ii = iotakkii(kk)

                    ; drhs(ii, ideriv) = +lAze(kk, ideriv)
                    if (NOTstellsym .and. kk > 0) then
                        drhs(ii + mns - 1, ideriv) = +lAzo(kk, ideriv)
                    end if

                    ; dmatrix(ii, 1, ideriv) = lAte(kk, ideriv)
                    if (NOTstellsym .and. kk > 0) then
                        dmatrix(ii + mns - 1, 1, ideriv) = lAto(kk, ideriv)
                    end if

                    do jj = 2, mns; mj = ims(jj); nj = ins(jj)

                        ii = iotakadd(kk, jj)

                        if (ii < 1) cycle

                        dmatrix(ii, jj, ideriv) = dmatrix(ii, jj, ideriv) + (-mj*lAze(kk, ideriv) + nj*lAte(kk, ideriv))*half
                        if (NOTstellsym) then
                            if (ii + mns - 1 < 1 .or. ii + mns - 1 > NN .or. jj + mns - 1 < 1 .or. jj + mns - 1 > NN) then
                                write (6, '("tr00ab :      fatal : myid=",i3," ; ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN ; illegal subscript;")') myid
                                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                                stop "tr00ab : ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN : illegal subscript ;"
                            end if
                            dmatrix(ii + mns - 1, jj, ideriv) = dmatrix(ii + mns - 1, jj, ideriv) + (-mj*lAzo(kk, ideriv) + nj*lAto(kk, ideriv))*half
                            dmatrix(ii, jj + mns - 1, ideriv) = dmatrix(ii, jj + mns - 1, ideriv) - (+mj*lAzo(kk, ideriv) - nj*lAto(kk, ideriv))*half
                            dmatrix(ii + mns - 1, jj + mns - 1, ideriv) = dmatrix(ii + mns - 1, jj + mns - 1, ideriv) + (+mj*lAze(kk, ideriv) - nj*lAte(kk, ideriv))*half
                        end if

                    end do

                    do jj = 2, mns; mj = ims(jj); nj = ins(jj)

                        ii = iotaksub(kk, jj)

                        if (ii < 1) cycle

                        if (ii > NN .or. jj > NN) then
                            write (6, '("tr00ab :      fatal : myid=",i3," ; ii.gt.NN .or. jj.gt.NN ; illegal subscript;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "tr00ab : ii.gt.NN .or. jj.gt.NN : illegal subscript ;"
                        end if
                        dmatrix(ii, jj, ideriv) = dmatrix(ii, jj, ideriv) + (-mj*lAze(kk, ideriv) + nj*lAte(kk, ideriv))*half
                        if (NOTstellsym) then
                            if (ii + mns - 1 < 1 .or. ii + mns - 1 > NN .or. jj + mns - 1 < 1 .or. jj + mns - 1 > NN) then
                                write (6, '("tr00ab :      fatal : myid=",i3," ; ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN ; illegal subscript;")') myid
                                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                                stop "tr00ab : ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN : illegal subscript ;"
                            end if
                            dmatrix(ii + mns - 1, jj, ideriv) = dmatrix(ii + mns - 1, jj, ideriv) + (-mj*lAzo(kk, ideriv) + nj*lAto(kk, ideriv))*half*iotaksgn(kk, jj)
                            dmatrix(ii, jj + mns - 1, ideriv) = dmatrix(ii, jj + mns - 1, ideriv) + (+mj*lAzo(kk, ideriv) - nj*lAto(kk, ideriv))*half
                            dmatrix(ii + mns - 1, jj + mns - 1, ideriv) = dmatrix(ii + mns - 1, jj + mns - 1, ideriv) - (+mj*lAze(kk, ideriv) - nj*lAte(kk, ideriv))*half*iotaksgn(kk, jj)
                        end if

                    end do

                end do
            end do

            call DCOPY(NN*NN, dmatrix(1, 1, 0), 1, omatrix(1, 1), 1)

            do jderiv = 0, 1

                if (iflag == 1 .and. jderiv /= 0) cycle

                select case (jderiv)
                case (0); 
                case (1); 
                    if (iflag == 2) then; call DGEMV('N', NN, NN, -one, dmatrix(1, 1, 1), NN, dlambda(1, 0), 1, one, drhs(1, 1), 1)
                        ; ; call DGEMV('N', NN, NN, -one, dmatrix(1, 1, 2), NN, dlambda(1, 0), 1, one, drhs(1, 2), 1)
                    end if
                    if (iflag == -1) then; call DGEMV('N', NN, NN, -one, dmatrix(1, 1, -1), NN, dlambda(1, 0), 1, one, drhs(1, -1), 1)
                    end if
                case default
                    if (.true.) then
                        write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; invalid jderiv;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "tr00ab : .true. : invalid jderiv ;"
                    end if
                end select

                lcpu = MPI_WTIME()

                select case (Lsvdiota)

                case (0)

                    if04aaf = 1

                    select case (jderiv)

                    case (0)

                        MM = 1

                        call dgesvx('N', 'N', NN, MM, dmatrix(1:NN, 1:NN, 0), NN, FAA(1:NN, 1:NN), NN, ipiv(1:NN), &
                                    equed, Rdgesvx(1:NN), Cdgesvx(1:NN), drhs(1:NN, 0:0), NN, dlambda(1:NN, 0:0), &
                                    NN, rcond, ferr, berr, work4(1:4*NN), iwork4(1:NN), idgesvx)

                        ; ldiota(innout, 0) = dlambda(1, 0)
                        ; dlambdaout(1:NN, lvol, innout) = dlambda(1:NN, 0)

                    case (1)

                        MM = 2
                        if (iflag == -1) then; drhs(1:NN, 1) = drhs(1:NN, -1)
                            ; ; drhs(1:NN, 2) = zero
                        end if

                        call DCOPY(NN*NN, omatrix(1, 1), 1, dmatrix(1, 1, 0), 1)

                        call dgesvx('N', 'N', NN, MM, dmatrix(1:NN, 1:NN, 0), NN, FAA(1:NN, 1:NN), NN, ipiv(1:NN), &
                                    equed, Rdgesvx(1:NN), Cdgesvx(1:NN), drhs(1:NN, 1:MM), NN, dlambda(1:NN, 1:MM), &
                                    NN, rcond, ferr2(1:MM), berr2(1:MM), work4(1:4*NN), iwork4(1:NN), idgesvx)

                        if (iflag == 2) ldiota(innout, 1:2) = dlambda(1, 1:2)
                        if (iflag == -1) ldiota(innout, -1) = dlambda(1, 1)

                    case default

                        if (.true.) then
                            write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; invalid jderiv;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "tr00ab : .true. : invalid jderiv ;"
                        end if

                    end select

                    cput = MPI_WTIME()

                    select case (idgesvx)
                    case (0); if (Wtr00ab) write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput - lcpu, "solved Fourier ; ", dlambda(1, 0)
                    case (1:); write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput - lcpu, "singular ;       "
                    case (:-1); write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput - lcpu, "input error ;    "
                    case default; 
                        if (.true.) then
                            write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; illegal ifail returned by dgesvx;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "tr00ab : .true. : illegal ifail returned by dgesvx ;"
                        end if
                    end select

                    if (idgesvx /= 0) then
                        write (6, '("tr00ab :      fatal : myid=",i3," ; idgesvx.ne.0 ; failed to construct straight-fieldline angle using dgesvx;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "tr00ab : idgesvx.ne.0 : failed to construct straight-fieldline angle using dgesvx ;"
                    end if

                case (1)

                    nlvl = max(0, int(log(real(NN)/26)/log(2.0d0)) + 1)
                    Lwork = (63 + 8*nlvl)*NN + 676
                    Liwork = max(1, 11*NN + 3*nlvl*NN)

                    if (allocated(work)) deallocate (work)
                    allocate (work(1:Lwork), stat=astat)
                    work(1:Lwork) = zero
                    if (allocated(iwork)) then
                        deallocate (iwork, stat=astat)
                    end if
                    if (allocated(iwork)) deallocate (iwork)
                    allocate (iwork(1:Liwork), stat=astat)
                    iwork(1:Liwork) = zero

                    select case (jderiv)

                    case (0)

                        dlambda(1:NN, 0) = drhs(1:NN, 0)

                        call dgelsd(NN, NN, 1, dmatrix(1:NN, 1:NN, 0), NN, dlambda(1:NN, 0), NN, sval(1:NN), rcond, Irank, &
                                    work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd)

                        ldiota(innout, 0) = dlambda(1, 0)
                        dlambdaout(1:NN, lvol, innout) = dlambda(1:NN, 0)

                    case (1)

                        if (iflag == 2) then
                            do imupf = 1, 2
                                dmatrix(1:NN, 1:NN, 0) = omatrix(1:NN, 1:NN); dlambda(1:NN, imupf) = drhs(1:NN, imupf)

                                call dgelsd(NN, NN, 1, dmatrix(1:NN, 1:NN, 0), NN, dlambda(1:NN, imupf), NN, sval(1:NN), rcond, Irank, &
                                            work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd)

                                ldiota(innout, imupf) = dlambda(1, imupf)
                            end do
                        elseif (iflag == -1) then
                            do imupf = -1, -1
                                dmatrix(1:NN, 1:NN, 0) = omatrix(1:NN, 1:NN); dlambda(1:NN, imupf) = drhs(1:NN, imupf)

                                call dgelsd(NN, NN, 1, dmatrix(1:NN, 1:NN, 0), NN, dlambda(1:NN, imupf), NN, sval(1:NN), rcond, Irank, &
                                            work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd)

                                ldiota(innout, imupf) = dlambda(1, imupf)
                            end do
                        else
                            if (.true.) then
                                write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; invalid iflag;")') myid
                                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                                stop "tr00ab : .true. : invalid iflag ;"
                            end if
                        end if

                    case default

                        if (.true.) then
                            write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; invalid jderiv;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "tr00ab : .true. : invalid jderiv ;"
                        end if

                    end select

                    deallocate (work, stat=astat)

                    cput = MPI_WTIME()

                    select case (idgelsd)
                    case (0); if (Wtr00ab) write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput - lcpu, "solved Fourier ; ", dlambda(1, 0)
                    case (:-1); write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput - lcpu, "input error ;    "
                    case (1:); write (ounit, 1030) cput - cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput - lcpu, "QR failed ;      "
                    case default
                        if (.true.) then
                            write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; illegal ifail returned by f04arf;")') myid
                            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                            stop "tr00ab : .true. : $illegal ifail returned by f04arf ;"
                        end if
                    end select

                    if (idgelsd /= 0) then
                        write (6, '("tr00ab :      fatal : myid=",i3," ; idgelsd.ne.0 ; failed to construct straight-fieldline angle using dgelsd;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "tr00ab : idgelsd.ne.0 : failed to construct straight-fieldline angle using dgelsd ;"
                    end if

                    dmatrix(1:NN, 1:NN, 0) = omatrix(1:NN, 1:NN)

                case default

                    if (.true.) then
                        write (6, '("tr00ab :      fatal : myid=",i3," ; .true. ; illegal Lsvdiota;")') myid
                        call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                        stop "tr00ab : .true. : illegal Lsvdiota ;"
                    end if

                end select

1030            format("tr00ab : ", f10.2, " ; myid=", i3, " ; lvol=", i3, " ; innout="i2" ; jderiv="i2" ; "a7"="i2" ; time="f10.4" ; "a17, :" [d]iota="es17.09" ;")

            end do

        end if

        if (Lsparse == 0 .or. Lsparse == 3) then
            deallocate (dmatrix, stat=astat)
            deallocate (omatrix, stat=astat)
            deallocate (FAA, stat=astat)
        end if

    end do

end subroutine tr00ab

