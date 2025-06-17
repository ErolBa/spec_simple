
subroutine curent(lvol, mn, Nt, Nz, iflag, ldItGp)

    use constants, only: zero, one, two, pi2

    use numerical, only:

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wcurent, Lrad

    use cputiming, only: Tcurent

    use allglobal, only: ncpu, cpus, myid, MPI_COMM_SPEC, &
                         Mvol, im, in, mne, ime, ine, &
                         YESstellsym, NOTstellsym, &
                         sg, guvij, &
                         Ntz, ijreal, ijimag, jireal, jiimag, &
                         efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                         Ate, Aze, Ato, Azo, TT, &
                         build_vector_potential

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, mn, Nt, Nz, iflag
    real(8), intent(out) :: ldItGp(0:1, -1:2)

    integer :: innout, ideriv, ii, ll, Lcurvature, ifail
    real(8) :: lss
    real(8) :: Bsupt(1:Nt*Nz, -1:2), Bsupz(1:Nt*Nz, -1:2)

    innout = 0; lss = two*innout - one

    do ideriv = -1, 2

        if (iflag == 1 .and. ideriv /= 0) cycle
        if (iflag == 2 .and. ideriv < 0) cycle
        if (iflag == -1 .and. ideriv > 0) cycle

        call build_vector_potential(lvol, innout, ideriv, 1)

        call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                    Nt, Nz, Bsupz(1:Ntz, ideriv), Bsupt(1:Ntz, ideriv))

    end do

    select case (iflag)
    case (-1); Lcurvature = 3
    case (1); Lcurvature = 1
    case (2); Lcurvature = 1
    end select

    call coords(lvol, lss, Lcurvature, Ntz, mn)

    do ideriv = -1, 2

        if (iflag == 1 .and. ideriv /= 0) cycle
        if (iflag == 2 .and. ideriv < 0) cycle
        if (iflag == -1 .and. ideriv > 0) cycle

        ijreal(1:Ntz) = (-Bsupt(1:Ntz, ideriv)*guvij(1:Ntz, 2, 2, 0) + Bsupz(1:Ntz, ideriv)*guvij(1:Ntz, 2, 3, 0))/sg(1:Ntz, 0)
        ijimag(1:Ntz) = (-Bsupt(1:Ntz, ideriv)*guvij(1:Ntz, 2, 3, 0) + Bsupz(1:Ntz, ideriv)*guvij(1:Ntz, 3, 3, 0))/sg(1:Ntz, 0)
        if (ideriv == -1) then
            ijreal(1:Ntz) = ijreal(1:Ntz) + (-Bsupt(1:Ntz, 0)*guvij(1:Ntz, 2, 2, 1) + Bsupz(1:Ntz, 0)*guvij(1:Ntz, 2, 3, 1))/sg(1:Ntz, 0)
            ijimag(1:Ntz) = ijimag(1:Ntz) + (-Bsupt(1:Ntz, 0)*guvij(1:Ntz, 2, 3, 1) + Bsupz(1:Ntz, 0)*guvij(1:Ntz, 3, 3, 1))/sg(1:Ntz, 0)
        end if

        ifail = 0
        call tfft(Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                  mne, ime(1:mne), ine(1:mne), efmn(1:mne), ofmn(1:mne), cfmn(1:mne), sfmn(1:mne), ifail)

        ldItGp(0, ideriv) = efmn(1)*pi2
        ldItGp(1, ideriv) = cfmn(1)*pi2

    end do

end subroutine curent

