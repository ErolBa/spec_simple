
subroutine metrix(lquad, lvol)

    use constants, only: zero, one

    use numerical, only: small

    use fileunits, only: ounit

    use inputlist, only: Wmetrix

    use cputiming, only: Tmetrix

    use allglobal, only: myid, ncpu, cpus, &
                         dBdX, &
                         mn, im, in, mne, ime, ine, &
                         Nt, Nz, Ntz, efmn, ofmn, cfmn, sfmn, & ! 10 Dec 15;
                         ijreal, & ! workspace;
                         sg, guvij, & ! calculated in coords;
                         gvuij, & ! this is workspace: nowhere used outside of this routine;
                         goomne, goomno, &
                         gssmne, gssmno, &
                         gstmne, gstmno, &
                         gszmne, gszmno, &
                         gttmne, gttmno, &
                         gtzmne, gtzmno, &
                         gzzmne, gzzmno, &
                         guvijsave

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, lquad

    integer :: Lcurvature, ifail, ideriv, jquad

    do jquad = 1, lquad

        gvuij(1:Ntz, 0, 0) = one

        gvuij(1:Ntz, 1, 1) = guvijsave(1:Ntz, 1, 1, jquad)
        gvuij(1:Ntz, 1, 2) = guvijsave(1:Ntz, 1, 2, jquad)
        gvuij(1:Ntz, 1, 3) = guvijsave(1:Ntz, 1, 3, jquad)
        gvuij(1:Ntz, 2, 2) = guvijsave(1:Ntz, 2, 2, jquad)
        gvuij(1:Ntz, 2, 3) = guvijsave(1:Ntz, 2, 3, jquad)
        gvuij(1:Ntz, 3, 3) = guvijsave(1:Ntz, 3, 3, jquad)

        ifail = 0
        call tfft(Nt, Nz, gvuij(1:Ntz, 0, 0), ijreal(1:Ntz), &
                  mne, ime(1:mne), ine(1:mne), goomne(1:mne, jquad), goomno(1:mne, jquad), cfmn(1:mne), sfmn(1:mne), ifail)
        goomne(0, jquad) = zero; goomno(0, jquad) = zero

        ifail = 0
        call tfft(Nt, Nz, gvuij(1:Ntz, 1, 1), gvuij(1:Ntz, 1, 2), &
                  mne, ime(1:mne), ine(1:mne), gssmne(1:mne, jquad), gssmno(1:mne, jquad), gstmne(1:mne, jquad), gstmno(1:mne, jquad), ifail)
        gssmne(0, jquad) = zero; gssmno(0, jquad) = zero
        gstmne(0, jquad) = zero; gstmno(0, jquad) = zero

        ifail = 0
        call tfft(Nt, Nz, gvuij(1:Ntz, 1, 3), gvuij(1:Ntz, 2, 2), &
                  mne, ime(1:mne), ine(1:mne), gszmne(1:mne, jquad), gszmno(1:mne, jquad), gttmne(1:mne, jquad), gttmno(1:mne, jquad), ifail)
        gszmne(0, jquad) = zero; gszmno(0, jquad) = zero
        gttmne(0, jquad) = zero; gttmno(0, jquad) = zero

        ifail = 0
        call tfft(Nt, Nz, gvuij(1:Ntz, 2, 3), gvuij(1:Ntz, 3, 3), &
                  mne, ime(1:mne), ine(1:mne), gtzmne(1:mne, jquad), gtzmno(1:mne, jquad), gzzmne(1:mne, jquad), gzzmno(1:mne, jquad), ifail)
        gtzmne(0, jquad) = zero; gtzmno(0, jquad) = zero
        gzzmne(0, jquad) = zero; gzzmno(0, jquad) = zero

    end do

end subroutine metrix

subroutine compute_guvijsave(lquad, vvol, ideriv, Lcurvature)

    use allglobal, only: gaussianabscissae, Ntz, mn, guvij, guvijsave, &
                         sg

    implicit none

    integer, intent(in) :: vvol, lquad, ideriv, Lcurvature
    integer :: jquad, ii, jj
    real(8) :: lss

    do jquad = 1, lquad
        lss = gaussianabscissae(jquad, vvol)
        call coords(vvol, lss, Lcurvature, Ntz, mn)
        guvijsave(1:Ntz, 1:3, 1:3, jquad) = guvij(1:Ntz, 1:3, 1:3, ideriv)
        do ii = 1, 3
            do jj = 1, 3
                guvijsave(1:Ntz, jj, ii, jquad) = guvijsave(1:Ntz, jj, ii, jquad)/sg(1:Ntz, 0)
            end do
        end do
    end do

end subroutine compute_guvijsave
