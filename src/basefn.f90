
subroutine get_cheby(lss, lrad, cheby)

    use constants, only: zero, one, two

    implicit none

    real(8), intent(in) :: lss
    integer, intent(in) :: lrad
    real(8), intent(inout) :: cheby(0:lrad, 0:1)

    integer :: ll

    cheby = zero

    ; cheby(0, 0:1) = (/one, zero/)
    ; cheby(1, 0:1) = (/lss, one/)
    do ll = 2, lrad; cheby(ll, 0:1) = (/two*lss*cheby(ll - 1, 0) - cheby(ll - 2, 0), two*cheby(ll - 1, 0) + two*lss*cheby(ll - 1, 1) - cheby(ll - 2, 1)/)
    end do

    do ll = 1, lrad
        cheby(ll, 0) = cheby(ll, 0) - (-1)**ll
    end do

    do ll = 0, lrad
        cheby(ll, 0:1) = cheby(ll, 0:1)/real(ll + 1)
    end do

    return
end subroutine get_cheby

subroutine get_cheby_d2(lss, lrad, cheby)

    use constants, only: zero, one, two

    implicit none

    real(8), intent(in) :: lss
    integer, intent(in) :: lrad
    real(8), intent(inout) :: cheby(0:lrad, 0:2)

    integer :: ll

    cheby = zero

    ; ; cheby(0, 0:2) = (/one, zero, zero/)
    ; ; cheby(1, 0:2) = (/lss, one, zero/)
    do ll = 2, lrad
        cheby(ll, 0:2) = (/two*lss*cheby(ll - 1, 0) - cheby(ll - 2, 0), &
                           two*cheby(ll - 1, 0) + two*lss*cheby(ll - 1, 1) - cheby(ll - 2, 1), &
                           two*cheby(ll - 1, 1) + two*cheby(ll - 1, 1) + two*lss*cheby(ll - 1, 2) - cheby(ll - 2, 2)/)
    end do

    do ll = 1, lrad
        cheby(ll, 0) = cheby(ll, 0) - (-1)**ll
    end do

    do ll = 0, lrad
        cheby(ll, 0:2) = cheby(ll, 0:2)/real(ll + 1)
    end do

    return
end subroutine get_cheby_d2

subroutine get_zernike(r, lrad, mpol, zernike)

    use constants, only: zero, one, two

    implicit none

    real(8), intent(in) :: r
    integer, intent(in) :: lrad, mpol
    real(8), intent(inout) :: zernike(0:lrad, 0:mpol, 0:1)

    real(8) :: rm, rm1
    real(8) :: factor1, factor2, factor3, factor4
    integer :: m, n

    rm = one
    rm1 = zero
    zernike(:, :, :) = zero
    do m = 0, mpol
        if (lrad >= m) then
            zernike(m, m, 0:1) = (/rm, real(m)*rm1/)
        end if

        if (lrad >= m + 2) then
            zernike(m + 2, m, 0) = real(m + 2)*rm*r**2 - real(m + 1)*rm
            zernike(m + 2, m, 1) = real((m + 2)**2)*rm*r - real((m + 1)*m)*rm1
        end if

        do n = m + 4, lrad, 2
            factor1 = real(n)/real(n**2 - m**2)
            factor2 = real(4*(n - 1))
            factor3 = real((n - 2 + m)**2)/real(n - 2) + real((n - m)**2)/real(n)
            factor4 = real((n - 2)**2 - m**2)/real(n - 2)

            zernike(n, m, 0) = factor1*((factor2*r**2 - factor3)*zernike(n - 2, m, 0) - factor4*zernike(n - 4, m, 0))
            zernike(n, m, 1) = factor1*(two*factor2*r*zernike(n - 2, m, 0) + (factor2*r**2 - factor3)*zernike(n - 2, m, 1) - factor4*zernike(n - 4, m, 1))
        end do

        rm1 = rm
        rm = rm*r

    end do

    do n = 2, lrad, 2
        zernike(n, 0, 0) = zernike(n, 0, 0) - (-1)**(n/2)
    end do

    if (mpol >= 1) then
        do n = 3, lrad, 2
            zernike(n, 1, 0) = zernike(n, 1, 0) - (-1)**((n - 1)/2)*real((n + 1)/2)*r
            zernike(n, 1, 1) = zernike(n, 1, 1) - (-1)**((n - 1)/2)*real((n + 1)/2)
        end do
    end if

    do m = 0, mpol
        do n = m, lrad, 2
            zernike(n, m, :) = zernike(n, m, :)/real(n + 1)
        end do
    end do
end subroutine get_zernike

subroutine get_zernike_d2(r, lrad, mpol, zernike)

    use constants, only: zero, one, two

    implicit none

    real(8), intent(in) :: r
    integer, intent(in) :: lrad, mpol
    real(8), intent(inout) :: zernike(0:lrad, 0:mpol, 0:2)

    real(8) :: rm, rm1, rm2
    real(8) :: factor1, factor2, factor3, factor4
    integer :: m, n

    rm = one
    rm1 = zero
    rm2 = zero
    zernike(:, :, :) = zero
    do m = 0, mpol
        if (lrad >= m) then
            zernike(m, m, 0:2) = (/rm, real(m)*rm1, real(m*(m - 1))*rm2/)
        end if

        if (lrad >= m + 2) then
            zernike(m + 2, m, 0) = real(m + 2)*rm*r**2 - real(m + 1)*rm
            zernike(m + 2, m, 1) = real((m + 2)**2)*rm*r - real((m + 1)*m)*rm1
            zernike(m + 2, m, 2) = real((m + 2)**2*(m + 1))*rm - real((m + 1)*m*(m - 1))*rm2
        end if

        do n = m + 4, lrad, 2
            factor1 = real(n)/real(n**2 - m**2)
            factor2 = real(4*(n - 1))
            factor3 = real((n - 2 + m)**2)/real(n - 2) + real((n - m)**2)/real(n)
            factor4 = real((n - 2)**2 - m**2)/real(n - 2)

            zernike(n, m, 0) = factor1*((factor2*r**2 - factor3)*zernike(n - 2, m, 0) - factor4*zernike(n - 4, m, 0))
            zernike(n, m, 1) = factor1*(two*factor2*r*zernike(n - 2, m, 0) + (factor2*r**2 - factor3)*zernike(n - 2, m, 1) - factor4*zernike(n - 4, m, 1))
            zernike(n, m, 2) = factor1*(two*factor2*(two*r*zernike(n - 2, m, 1) + zernike(n - 2, m, 0)) &
                                        + (factor2*r**2 - factor3)*zernike(n - 2, m, 2) - factor4*zernike(n - 4, m, 2))
        end do

        rm2 = rm1
        rm1 = rm
        rm = rm*r

    end do
    do n = 2, lrad, 2
        zernike(n, 0, 0) = zernike(n, 0, 0) - (-1)**(n/2)
    end do
    if (mpol >= 1) then
        do n = 3, lrad, 2
            zernike(n, 1, 0) = zernike(n, 1, 0) - (-1)**((n - 1)/2)*real((n + 1)/2)*r
            zernike(n, 1, 1) = zernike(n, 1, 1) - (-1)**((n - 1)/2)*real((n + 1)/2)
        end do
    end if

    do m = 0, mpol
        do n = m, lrad, 2
            zernike(n, m, :) = zernike(n, m, :)/real(n + 1)
        end do
    end do
end subroutine get_zernike_d2

subroutine get_zernike_rm(r, lrad, mpol, zernike)

    use constants, only: zero, one, two

    implicit none

    real(8), intent(in) :: r
    integer, intent(in) :: lrad, mpol
    real(8), intent(inout) :: zernike(0:lrad, 0:mpol)

    real(8) :: factor1, factor2, factor3, factor4
    integer :: m, n

    zernike(:, :) = zero
    do m = 0, mpol
        if (lrad >= m) then
            zernike(m, m) = one
        end if

        if (lrad >= m + 2) then
            zernike(m + 2, m) = real(m + 2)*r**2 - real(m + 1)
        end if

        do n = m + 4, lrad, 2
            factor1 = real(n)/real(n**2 - m**2)
            factor2 = real(4*(n - 1))
            factor3 = real((n - 2 + m)**2)/real(n - 2) + real((n - m)**2)/real(n)
            factor4 = real((n - 2)**2 - m**2)/real(n - 2)

            zernike(n, m) = factor1*((factor2*r**2 - factor3)*zernike(n - 2, m) - factor4*zernike(n - 4, m))
        end do

    end do
    do n = 2, lrad, 2
        zernike(n, 0) = zernike(n, 0) - (-1)**(n/2)
    end do
    if (mpol >= 1) then
        do n = 3, lrad, 2
            zernike(n, 1) = zernike(n, 1) - (-1)**((n - 1)/2)*real((n + 1)/2)
        end do
    end if

    do m = 0, mpol
        do n = m, lrad, 2
            zernike(n, m) = zernike(n, m)/real(n + 1)
        end do
    end do
end subroutine get_zernike_rm
