
subroutine gi00ab(Mpol, Ntor, Nfp, mn, im, in)

    implicit none

    integer, intent(in) :: Mpol, Ntor, Nfp, mn
    integer, intent(out) :: im(mn), in(mn)

    integer :: imn, mm, nn

    imn = 0

    ; mm = 0
    ; do nn = 0, Ntor
        ; imn = imn + 1; im(imn) = mm; in(imn) = nn*Nfp
        ; end do
    ; 
    do mm = 1, Mpol
        do nn = -Ntor, Ntor
            imn = imn + 1; im(imn) = mm; in(imn) = nn*Nfp
        end do
    end do

    return

end subroutine gi00ab

subroutine getimn(Mpol, Ntor, Nfp, mi, ni, idx)
    implicit none
    integer, intent(in) :: Mpol, Ntor, Nfp, mi, ni
    integer, intent(out) :: idx

    if (mi > Mpol .or. mi < 0 .or. ni > Ntor*Nfp .or. ni < -Ntor*Nfp) then
        idx = 0
    elseif (mi == 0) then
        idx = 1 + ni/Nfp
    else
        idx = 1 + Ntor + (2*Ntor + 1)*(mi - 1) + (ni/Nfp + Ntor + 1)
    end if

end subroutine getimn

subroutine tfft(Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn, ifail)

    use constants, only: half, zero, pi2

    use fileunits, only: ounit

    use inputlist, only: Nfp
    use allglobal, only: pi2nfp

    use fftw_interface
#ifdef OPENMP
    use OMP_LIB
#endif
    implicit none

    intrinsic aimag

    integer :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, imn, ifail, mm, nn
    real(8) :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)

    logical :: Lcheck = .false.
    integer :: jj, kk, ithread
    real(8) :: arg, ca, sa
    complex(c_double_complex) :: z1, z2, z3

    ithread = 1

    do jj = 1, Nz; cplxin(:, jj, ithread) = cmplx(ijreal((jj - 1)*Nt + 1:jj*Nt), ijimag((jj - 1)*Nt + 1:jj*Nt), KIND=c_double_complex)
    end do

    call fftw_execute_dft(planf, cplxin(:, :, ithread), cplxout(:, :, ithread)) !Forward transform
    Ntz = Nt*Nz
    cplxout(:, :, ithread) = cplxout(:, :, ithread)/Ntz
    cplxout(1, 1, ithread) = half*cplxout(1, 1, ithread)

    do imn = 1, mn
        mm = im(imn); nn = in(imn)/Nfp

        z1 = cplxout(1 + mod(Nt - mm, Nt), 1 + mod(Nz + nn, Nz), ithread)
        z2 = cplxout(1 + mm, 1 + mod(Nz - nn, Nz), ithread)

        z3 = z1 + z2
        efmn(imn) = real(z3); cfmn(imn) = aimag(z3)

        z3 = z1 - z2
        ofmn(imn) = aimag(z3); sfmn(imn) = -real(z3)
    end do

    if (.not. Lcheck) return

    ijreal(1:Ntz) = zero; ijimag(1:Ntz) = zero

    do jj = 0, Nt - 1

        do kk = 0, Nz - 1

            do imn = 1, mn; arg = im(imn)*jj*pi2/Nt - in(imn)*kk*pi2nfp/Nz; ca = cos(arg); sa = sin(arg)

                ijreal(1 + jj + kk*Nt) = ijreal(1 + jj + kk*Nt) + efmn(imn)*ca + ofmn(imn)*sa
                ijimag(1 + jj + kk*Nt) = ijimag(1 + jj + kk*Nt) + cfmn(imn)*ca + sfmn(imn)*sa

            end do
        end do
    end do

    return

end subroutine tfft

subroutine invfft(mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag)

    use constants, only: zero, two, half
    use inputlist, only: Nfp
    use fftw_interface
#ifdef OPENMP
    use OMP_LIB
#endif

    implicit none

    integer, intent(in) :: mn, im(mn), in(mn)
    real(8), intent(in) :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
    integer, intent(in) :: Nt, Nz
    real(8), intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;

    integer :: imn, jj, mm, nn, ithread

    ithread = 1

    cplxin(:, :, ithread) = zero

    do imn = 1, mn; mm = im(imn); nn = in(imn)/Nfp
        cplxin(1 + mod(Nt - mm, Nt), 1 + mod(Nz + nn, Nz), ithread) = &
            half*cmplx(efmn(imn) - sfmn(imn), cfmn(imn) + ofmn(imn), KIND=c_double_complex)
        cplxin(1 + mm, 1 + mod(Nz - nn, Nz), ithread) = &
            half*cmplx(efmn(imn) + sfmn(imn), cfmn(imn) - ofmn(imn), KIND=c_double_complex)
    end do
    cplxin(1, 1, ithread) = two*cplxin(1, 1, ithread)

    call fftw_execute_dft(planb, cplxin(:, :, ithread), cplxout(:, :, ithread)) !Inverse transform

    do jj = 1, Nz
        ijreal((jj - 1)*Nt + 1:jj*Nt) = real(cplxout(:, jj, ithread))
        ijimag((jj - 1)*Nt + 1:jj*Nt) = aimag(cplxout(:, jj, ithread))
    end do

    return

end subroutine invfft

subroutine gauleg(n, weight, abscis, ifail)

    use constants, only: zero, one, two, pi

    implicit none

    intrinsic abs, cos, epsilon

    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: weight, abscis
    integer, intent(out) :: ifail

    integer, parameter :: maxiter = 16
    integer :: m, j, i, irefl, iter
    real(8) :: z1, z, pp, p3, p2, p1
    real(8), parameter :: eps = epsilon(z)

    if (n < 1) then; ifail = 2; return
    end if

    m = (n + 1)/2 !Roots are symmetric in interval, so we only need half
    do i = 1, m !Loop over desired roots
        irefl = n + 1 - i
        if (i /= irefl) then
            z = cos(pi*(i - 0.25)/(n + 0.5)) ! Approximate ith root
        else !For an odd number of abscissae, the center must be at zero by symmetry.
            z = 0.0
        end if

        do iter = 1, maxiter
            p1 = one; p2 = zero ! Initialize recurrence relation

            do j = 1, n !Recurrence relation to get P(x)
                p3 = p2; p2 = p1
                p1 = ((two*j - one)*z*p2 - (j - one)*p3)/j
            end do !j

            pp = n*(z*p1 - p2)/(z*z - one) !Derivative of P(x)
            z1 = z; z = z1 - p1/pp !Newton iteration
            if (abs(z - z1) <= eps) exit !Convergence test
        end do !iter
        if (iter > maxiter) then
            ifail = 1; return
        end if

        abscis(i) = -z; abscis(irefl) = z
        weight(i) = two/((one - z*z)*pp*pp)
        weight(irefl) = weight(i)
    end do !i

    ifail = 0
end subroutine gauleg

#ifdef DELETETHIS

real(8) function pythag(a, b)
    implicit none
    real(8) :: a, b
    real(8) :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
        pythag = absa*sqrt(1.+(absb/absa)**2)
    else
        if (absb == 0.) then
            pythag = 0.
        else
            pythag = absb*sqrt(1.+(absa/absb)**2)
        end if
    end if
    return
end function pythag

#endif

