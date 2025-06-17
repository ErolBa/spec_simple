
subroutine lbpol(lvol, Bt00, ideriv, iocons)

    use constants, only: mu0, pi, pi2, two, one, half, zero

    use allglobal, only: Ate, Aze, Ato, Azo, TT, &
                         YESstellsym, NOTstellsym, &
                         im, in, mne, ime, ine, Mvol, mn, &
                         sg, guvij, &
                         Ntz, Lcoordinatesingularity, &
                         efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                         Nt, Nz, &
                         regumm, &
                         cpus, myid, dBdX, &
                         build_vector_potential

    use inputlist, only: Lrad, Wlbpol, Igeometry, Lcheck

    use fileunits, only: ounit

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer :: Lcurvature, ideriv, ii, ll, ifail, lvol, mi, ni, iocons
    real(8) :: lss, Bt00(1:Mvol, 0:1, -1:2)
    real(8) :: lAte(1:mn), lAze(1:mn), lAto(1:mn), lAzo(1:mn)
    real(8) :: dAt(1:Ntz), dAz(1:Ntz), Bt(1:Ntz), Bz(1:Ntz), dAt0(1:Ntz), dAz0(1:Ntz)
    real(8) :: dBtzero
    real(8) :: mfactor
    logical :: LGeometricDerivative

    lss = two*iocons - one

    Lcurvature = 1
    call coords(lvol, lss, Lcurvature, Ntz, mn)

    call build_vector_potential(lvol, iocons, ideriv, 1)

    call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz), dAz(1:Ntz))

    Bt(1:Ntz) = (-dAz(1:Ntz)*guvij(1:Ntz, 2, 2, 0) + dAt(1:Ntz)*guvij(1:Ntz, 2, 3, 0))/sg(1:Ntz, 0)
    Bz(1:Ntz) = (-dAz(1:Ntz)*guvij(1:Ntz, 2, 3, 0) + dAt(1:Ntz)*guvij(1:Ntz, 3, 3, 0))/sg(1:Ntz, 0)

    if (ideriv == -1) then

        Lcurvature = 3
        call coords(lvol, lss, Lcurvature, Ntz, mn)

        call build_vector_potential(lvol, iocons, 0, 1)

        call invfft(mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt0(1:Ntz), dAz0(1:Ntz))

        Bt(1:Ntz) = Bt(1:Ntz) + (-dAz0(1:Ntz)*guvij(1:Ntz, 2, 2, 1) + dAt0(1:Ntz)*guvij(1:Ntz, 2, 3, 1))/sg(1:Ntz, 0)
        Bz(1:Ntz) = Bz(1:Ntz) + (-dAz0(1:Ntz)*guvij(1:Ntz, 2, 3, 1) + dAt0(1:Ntz)*guvij(1:Ntz, 3, 3, 1))/sg(1:Ntz, 0)

    end if

    ifail = 0
    call tfft(Nt, Nz, Bt(1:Ntz), Bz(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail)

    Bt00(lvol, iocons, ideriv) = efmn(1)

5555 continue

end subroutine lbpol

