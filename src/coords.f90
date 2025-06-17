
subroutine coords(lvol, lss, Lcurvature, Ntz, mn)

    use constants, only: zero, half, one, two, pi2

    use numerical, only: vsmall, small

    use fileunits, only: ounit

    use inputlist, only: Wcoords, Igeometry, Ntor, rpol, rtor

    use allglobal, only: myid, cpus, pi2nfp, &
                         Mvol, im, in, halfmm, &
                         iRbc, iZbs, iRbs, iZbc, &
                         NOTstellsym, Lcoordinatesingularity, &
                         Nt, Nz, &
                         Rij, Zij, &
                         cosi, sini, &
                         sg, guvij, &
                         dBdX, &
                         dRodR, dRodZ, dZodR, dZodZ, &
                         Remn_ext, Romn_ext, Zemn_ext, Zomn_ext, Iquad, gaussianabscissae, use_ext_mesh

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, Lcurvature, Ntz, mn
    real(8), intent(in) :: lss

    integer :: ii, jj, kk, irz, innout, issym, signlss, mi, ni, imn, radial_index
    real(8) :: Remn(1:mn, 0:2), Zomn(1:mn, 0:2), Romn(1:mn, 0:2), Zemn(1:mn, 0:2), alss, blss, sbar, sbarhim(1:mn), fj(1:mn, 0:2)

    real(8) :: Dij(1:Ntz, 0:3), dguvij(1:Ntz, 1:3, 1:3), DRxij(1:Ntz, 0:3), DZxij(1:Ntz, 0:3)

    Rij(1:Ntz, 0:3, 0:3) = zero; sg(1:Ntz, 0:3) = zero; guvij(1:Ntz, 1:3, 1:3, 0:3) = zero
    Zij(1:Ntz, 0:3, 0:3) = zero

    Remn(1:mn, 0:2) = zero
    Zomn(1:mn, 0:2) = zero
    Romn(1:mn, 0:2) = zero
    Zemn(1:mn, 0:2) = zero

    if (Lcoordinatesingularity) then

        if (use_ext_mesh) then

            if (.not. allocated(Remn_ext)) write (*, *) "WARNING: Remn_ext not allocated"

            radial_index = minloc(abs(gaussianabscissae(:, 1) - lss), DIM=1)

            Remn(1:mn, 0:2) = Remn_ext(1:mn, radial_index, 0:2)
            Romn(1:mn, 0:2) = Romn_ext(1:mn, radial_index, 0:2)
            Zemn(1:mn, 0:2) = Zemn_ext(1:mn, radial_index, 0:2)
            Zomn(1:mn, 0:2) = Zomn_ext(1:mn, radial_index, 0:2)

        else

            sbar = (lss + one)*half

            select case (Igeometry)
            case (2); fj(1:Ntor + 1, 0) = sbar
                ; ; fj(Ntor + 2:mn, 0) = sbar**(im(Ntor + 2:mn) + 1)
            case (3); fj(1:Ntor + 1, 0) = sbar**2
                ; ; fj(Ntor + 2:mn, 0) = sbar**im(Ntor + 2:mn)
            case default
                if (.true.) then
                    write (6, '("coords :      fatal : myid=",i3," ; .true. ; invalid Igeometry for Lcoordinatesingularity=T;")') myid

                    stop "coords : .true. : invalid Igeometry for Lcoordinatesingularity=T ;"
                end if
            end select

            Remn(1:mn, 0) = iRbc(1:mn, 0) + (iRbc(1:mn, 1) - iRbc(1:mn, 0))*fj(1:mn, 0)
            if (NOTstellsym) then
                Romn(1:mn, 0) = iRbs(1:mn, 0) + (iRbs(1:mn, 1) - iRbs(1:mn, 0))*fj(1:mn, 0)
            end if
            if (Igeometry == 3) then
                Zomn(1:mn, 0) = iZbs(1:mn, 0) + (iZbs(1:mn, 1) - iZbs(1:mn, 0))*fj(1:mn, 0)
                if (NOTstellsym) then
                    Zemn(1:mn, 0) = iZbc(1:mn, 0) + (iZbc(1:mn, 1) - iZbc(1:mn, 0))*fj(1:mn, 0)
                end if
            end if

        end if

    else

        alss = half*(one - lss); blss = half*(one + lss)

        Remn(1:mn, 0) = alss*iRbc(1:mn, lvol - 1) + blss*iRbc(1:mn, lvol)
        if (NOTstellsym) then
            Romn(1:mn, 0) = alss*iRbs(1:mn, lvol - 1) + blss*iRbs(1:mn, lvol)
        end if
        if (Igeometry == 3) then
            Zomn(1:mn, 0) = alss*iZbs(1:mn, lvol - 1) + blss*iZbs(1:mn, lvol)
            if (NOTstellsym) then
                Zemn(1:mn, 0) = alss*iZbc(1:mn, lvol - 1) + blss*iZbc(1:mn, lvol)
            end if
        end if

    end if

    call invfft(mn, im(1:mn), in(1:mn), Remn(1:mn, 0), Romn(1:mn, 0), Zemn(1:mn, 0), Zomn(1:mn, 0), &
                Nt, Nz, Rij(1:Ntz, 0, 0), Zij(1:Ntz, 0, 0))

    if (Lcurvature == 0) return

    if (Lcoordinatesingularity) then

        if (.not. use_ext_mesh) then

            select case (Igeometry)
            case (2); fj(1:Ntor + 1, 1) = half
                ; ; fj(Ntor + 2:mn, 1) = half*(im(Ntor + 2:mn) + one)*fj(Ntor + 2:mn, 0)/sbar
            case (3); fj(1:Ntor + 1, 1) = sbar
                ; ; fj(Ntor + 2:mn, 1) = half*im(Ntor + 2:mn)*fj(Ntor + 2:mn, 0)/sbar
            case default
                if (.true.) then
                    write (6, '("coords :      fatal : myid=",i3," ; .true. ; invalid Igeometry for Lcoordinatesingularity=T and Lcurvature.ne.0;")') myid

                    stop "coords : .true. : invalid Igeometry for Lcoordinatesingularity=T and Lcurvature.ne.0 ;"
                end if
            end select

            Remn(1:mn, 1) = (iRbc(1:mn, 1) - iRbc(1:mn, 0))*fj(1:mn, 1)
            if (NOTstellsym) then
                Romn(1:mn, 1) = (iRbs(1:mn, 1) - iRbs(1:mn, 0))*fj(1:mn, 1)
            end if
            if (Igeometry == 3) then
                Zomn(1:mn, 1) = (iZbs(1:mn, 1) - iZbs(1:mn, 0))*fj(1:mn, 1)
                if (NOTstellsym) then
                    Zemn(1:mn, 1) = (iZbc(1:mn, 1) - iZbc(1:mn, 0))*fj(1:mn, 1)
                end if
            end if
        end if

    else

        Remn(1:mn, 1) = (-iRbc(1:mn, lvol - 1) + iRbc(1:mn, lvol))*half
        if (NOTstellsym) then
            Romn(1:mn, 1) = (-iRbs(1:mn, lvol - 1) + iRbs(1:mn, lvol))*half
        end if
        if (Igeometry == 3) then
            Zomn(1:mn, 1) = (-iZbs(1:mn, lvol - 1) + iZbs(1:mn, lvol))*half
            if (NOTstellsym) then
                Zemn(1:mn, 1) = (-iZbc(1:mn, lvol - 1) + iZbc(1:mn, lvol))*half
            end if
        end if

    end if

    call invfft(mn, im(1:mn), in(1:mn), Remn(1:mn, 1), Romn(1:mn, 1), Zemn(1:mn, 1), Zomn(1:mn, 1), &
                Nt, Nz, Rij(1:Ntz, 1, 0), Zij(1:Ntz, 1, 0))

    call invfft(mn, im(1:mn), in(1:mn), im(1:mn)*Romn(1:mn, 0), -im(1:mn)*Remn(1:mn, 0), im(1:mn)*Zomn(1:mn, 0), -im(1:mn)*Zemn(1:mn, 0), &
                Nt, Nz, Rij(1:Ntz, 2, 0), Zij(1:Ntz, 2, 0))

    call invfft(mn, im(1:mn), in(1:mn), -in(1:mn)*Romn(1:mn, 0), in(1:mn)*Remn(1:mn, 0), -in(1:mn)*Zomn(1:mn, 0), in(1:mn)*Zemn(1:mn, 0), &
                Nt, Nz, Rij(1:Ntz, 3, 0), Zij(1:Ntz, 3, 0))

    do ii = 1, 3; Rij(1:Ntz, 0, ii) = Rij(1:Ntz, ii, 0)
        ; ; Zij(1:Ntz, 0, ii) = Zij(1:Ntz, ii, 0)
    end do

    guvij(1:Ntz, 0, 0, 0) = one

    select case (Igeometry)

    case (1)

        sg(1:Ntz, 0) = Rij(1:Ntz, 1, 0)*rpol*rtor

        do ii = 1, 3
            do jj = ii, 3; guvij(1:Ntz, ii, jj, 0) = Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, 0)
            end do
        end do

        guvij(1:Ntz, 2, 2, 0) = guvij(1:Ntz, 2, 2, 0) + rpol*rpol
        guvij(1:Ntz, 3, 3, 0) = guvij(1:Ntz, 3, 3, 0) + rtor*rtor

    case (2)

        sg(1:Ntz, 0) = Rij(1:Ntz, 1, 0)*Rij(1:Ntz, 0, 0)

        do ii = 1, 3
            do jj = ii, 3; guvij(1:Ntz, ii, jj, 0) = Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, 0)
            end do
        end do

        guvij(1:Ntz, 2, 2, 0) = guvij(1:Ntz, 2, 2, 0) + Rij(1:Ntz, 0, 0)*Rij(1:Ntz, 0, 0)
        guvij(1:Ntz, 3, 3, 0) = guvij(1:Ntz, 3, 3, 0) + one

    case (3)

        sg(1:Ntz, 0) = Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, 0))

        do ii = 1, 3
            do jj = ii, 3; guvij(1:Ntz, ii, jj, 0) = Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, 0) + Zij(1:Ntz, ii, 0)*Zij(1:Ntz, jj, 0)
            end do
        end do

        guvij(1:Ntz, 3, 3, 0) = guvij(1:Ntz, 3, 3, 0) + (Rij(1:Ntz, 0, 0))**2

    case default

        if (.true.) then
            write (6, '("coords :      fatal : myid=",i3," ; .true. ; selected Igeometry not supported;")') myid

            stop "coords : .true. : selected Igeometry not supported ;"
        end if

    end select

    do ii = 2, 3
        do jj = 1, ii - 1; guvij(1:Ntz, ii, jj, 0) = guvij(1:Ntz, jj, ii, 0)
        end do
    end do

    if (Lcurvature <= 1) return

    select case (Lcurvature)

    case (2)

        if (Lcoordinatesingularity) then

            if (.not. use_ext_mesh) then

                select case (Igeometry)
                case (2); fj(1:Ntor + 1, 2) = zero
                    ; ; fj(Ntor + 2:mn, 2) = half*(im(Ntor + 2:mn))*fj(Ntor + 2:mn, 1)/sbar
                case (3); fj(1:Ntor + 1, 2) = half
                    ; ; fj(Ntor + 2:mn, 2) = half*(im(Ntor + 2:mn) - one)*fj(Ntor + 2:mn, 1)/sbar
                case default; 
                    ; 
                    if (.true.) then
                        write (6, '("coords :      fatal : myid=",i3," ; .true ; invalid Igeometry for Lcoordinatesingularity=T and Lcurvature=2;")') myid

                        stop "coords : .true : invalid Igeometry for Lcoordinatesingularity=T and Lcurvature=2 ;"
                    end if
                end select; 
                Remn(1:mn, 2) = (iRbc(1:mn, 1) - iRbc(1:mn, 0))*fj(1:mn, 2)
                if (NOTstellsym) then
                    Romn(1:mn, 2) = (iRbs(1:mn, 1) - iRbs(1:mn, 0))*fj(1:mn, 2)
                end if
                if (Igeometry == 3) then
                    Zomn(1:mn, 2) = (iZbs(1:mn, 1) - iZbs(1:mn, 0))*fj(1:mn, 2)
                    if (NOTstellsym) then
                        Zemn(1:mn, 2) = (iZbc(1:mn, 1) - iZbc(1:mn, 0))*fj(1:mn, 2)
                    end if
                end if

            end if

        else

            Remn(1:mn, 2) = zero
            if (NOTstellsym) then
                Romn(1:mn, 2) = zero
            end if
            if (Igeometry == 3) then
                Zomn(1:mn, 2) = zero
                if (NOTstellsym) then
                    Zemn(1:mn, 2) = zero
                end if
            end if

        end if

        call invfft(mn, im(1:mn), in(1:mn), &
                    Remn(1:mn, 2), Romn(1:mn, 2), Zemn(1:mn, 2), Zomn(1:mn, 2), &
                    Nt, Nz, Rij(1:Ntz, 1, 1), Zij(1:Ntz, 1, 1))

        call invfft(mn, im(1:mn), in(1:mn), &
                    +im(1:mn)*Romn(1:mn, 1), -im(1:mn)*Remn(1:mn, 1), im(1:mn)*Zomn(1:mn, 1), -im(1:mn)*Zemn(1:mn, 1), &
                    Nt, Nz, Rij(1:Ntz, 1, 2), Zij(1:Ntz, 1, 2))

        call invfft(mn, im(1:mn), in(1:mn), &
                    -in(1:mn)*Romn(1:mn, 1), in(1:mn)*Remn(1:mn, 1), -in(1:mn)*Zomn(1:mn, 1), in(1:mn)*Zemn(1:mn, 1), &
                    Nt, Nz, Rij(1:Ntz, 1, 3), Zij(1:Ntz, 1, 3))

        call invfft(mn, im(1:mn), in(1:mn), &
                    -im(1:mn)*im(1:mn)*Remn(1:mn, 0), -im(1:mn)*im(1:mn)*Romn(1:mn, 0), -im(1:mn)*im(1:mn)*Zemn(1:mn, 0), -im(1:mn)*im(1:mn)*Zomn(1:mn, 0), &
                    Nt, Nz, Rij(1:Ntz, 2, 2), Zij(1:Ntz, 2, 2))

        call invfft(mn, im(1:mn), in(1:mn), &
                    +im(1:mn)*in(1:mn)*Remn(1:mn, 0), im(1:mn)*in(1:mn)*Romn(1:mn, 0), im(1:mn)*in(1:mn)*Zemn(1:mn, 0), im(1:mn)*in(1:mn)*Zomn(1:mn, 0), &
                    Nt, Nz, Rij(1:Ntz, 2, 3), Zij(1:Ntz, 2, 3))

        call invfft(mn, im(1:mn), in(1:mn), &
                    -in(1:mn)*in(1:mn)*Remn(1:mn, 0), -in(1:mn)*in(1:mn)*Romn(1:mn, 0), -in(1:mn)*in(1:mn)*Zemn(1:mn, 0), -in(1:mn)*in(1:mn)*Zomn(1:mn, 0), &
                    Nt, Nz, Rij(1:Ntz, 3, 3), Zij(1:Ntz, 3, 3))

        do ii = 2, 3
            do jj = 1, ii - 1; Rij(1:Ntz, ii, jj) = Rij(1:Ntz, jj, ii); Zij(1:Ntz, ii, jj) = Zij(1:Ntz, jj, ii)
            end do
        end do

        select case (Igeometry)

        case (1)

            do kk = 1, 3

                sg(1:Ntz, kk) = Rij(1:Ntz, 1, kk)*rpol*rtor

                do ii = 1, 3
                    do jj = ii, 3; guvij(1:Ntz, ii, jj, kk) = Rij(1:Ntz, ii, kk)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, kk)
                    end do
                end do

            end do

        case (2)

            do kk = 1, 3

                sg(1:Ntz, kk) = Rij(1:Ntz, 1, kk)*Rij(1:Ntz, 0, 0) + Rij(1:Ntz, 1, 0)*Rij(1:Ntz, 0, kk)

                do ii = 1, 3
                    do jj = ii, 3; guvij(1:Ntz, ii, jj, kk) = Rij(1:Ntz, ii, kk)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, kk)
                    end do
                end do

                guvij(1:Ntz, 2, 2, kk) = guvij(1:Ntz, 2, 2, kk) + Rij(1:Ntz, 0, kk)*Rij(1:Ntz, 0, 0) + Rij(1:Ntz, 0, 0)*Rij(1:Ntz, 0, kk)

            end do

        case (3)

            do kk = 1, 3

                sg(1:Ntz, kk) = Rij(1:Ntz, kk, 0)*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, 0)) &
                                + Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 1, kk)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, kk)*Zij(1:Ntz, 2, 0)) &
                                + Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, kk) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, kk))

                sg(1:Ntz, kk) = sg(1:Ntz, kk)

                do ii = 1, 3
                    do jj = ii, 3
                        guvij(1:Ntz, ii, jj, kk) = Rij(1:Ntz, ii, kk)*Rij(1:Ntz, jj, 0) &
                                                   + Rij(1:Ntz, ii, 0)*Rij(1:Ntz, jj, kk) &
                                                   + Zij(1:Ntz, ii, kk)*Zij(1:Ntz, jj, 0) &
                                                   + Zij(1:Ntz, ii, 0)*Zij(1:Ntz, jj, kk)
                    end do
                end do

                guvij(1:Ntz, 3, 3, kk) = guvij(1:Ntz, 3, 3, kk) + (Rij(1:Ntz, 0, kk)*Rij(1:Ntz, 0, 0) + Rij(1:Ntz, 0, 0)*Rij(1:Ntz, 0, kk))

            end do

        case default

            if (.true.) then
                write (6, '("coords :      fatal : myid=",i3," ; .true. ; selected Igeometry not supported for Lcurvature.eq.2;")') myid

                stop "coords : .true. : selected Igeometry not supported for Lcurvature.eq.2 ;"
            end if

        end select

        do ii = 2, 3
            do jj = 1, ii - 1; guvij(1:Ntz, ii, jj, 1:3) = guvij(1:Ntz, jj, ii, 1:3)
            end do
        end do

    case (3, 4, 5)

        ii = dBdX%ii; innout = dBdX%innout; irz = dBdX%irz; issym = dBdX%issym

        if (Lcoordinatesingularity) then

            if ((irz == 0 .and. issym == 0) .or. (irz == 1 .and. issym == 1)) then

                Dij(1:Ntz, 0) = fj(ii, 0)*cosi(1:Ntz, ii)
                Dij(1:Ntz, 1) = fj(ii, 1)*cosi(1:Ntz, ii)
                Dij(1:Ntz, 2) = fj(ii, 0)*sini(1:Ntz, ii)*(-im(ii))
                Dij(1:Ntz, 3) = fj(ii, 0)*sini(1:Ntz, ii)*(+in(ii))

            else

                Dij(1:Ntz, 0) = fj(ii, 0)*sini(1:Ntz, ii)
                Dij(1:Ntz, 1) = fj(ii, 1)*sini(1:Ntz, ii)
                Dij(1:Ntz, 2) = fj(ii, 0)*cosi(1:Ntz, ii)*(+im(ii))
                Dij(1:Ntz, 3) = fj(ii, 0)*cosi(1:Ntz, ii)*(-in(ii))

            end if

            if (Lcurvature == 3 .or. Lcurvature == 4) then
                if (Igeometry == 3) then
                    if (irz == 0) then
                        DRxij(1:Ntz, 0) = (one - sbar**2)*dRodR(1:Ntz, issym, ii)
                        DRxij(1:Ntz, 1) = -sbar*dRodR(1:Ntz, issym, ii)
                        DRxij(1:Ntz, 2) = zero
                        DRxij(1:Ntz, 3) = (one - sbar**2)*dRodR(1:Ntz, issym + 2, ii)

                        DZxij(1:Ntz, 0) = (one - sbar**2)*dZodR(1:Ntz, issym, ii)
                        DZxij(1:Ntz, 1) = -sbar*dZodR(1:Ntz, issym, ii)
                        DZxij(1:Ntz, 2) = zero
                        DZxij(1:Ntz, 3) = (one - sbar**2)*dZodR(1:Ntz, issym + 2, ii)
                    else
                        DRxij(1:Ntz, 0) = (one - sbar**2)*dRodZ(1:Ntz, 1 - issym, ii)
                        DRxij(1:Ntz, 1) = -sbar*dRodZ(1:Ntz, 1 - issym, ii)
                        DRxij(1:Ntz, 2) = zero
                        DRxij(1:Ntz, 3) = (one - sbar**2)*dRodZ(1:Ntz, 1 - issym + 2, ii)

                        DZxij(1:Ntz, 0) = (one - sbar**2)*dZodZ(1:Ntz, 1 - issym, ii)
                        DZxij(1:Ntz, 1) = -sbar*dZodZ(1:Ntz, 1 - issym, ii)
                        DZxij(1:Ntz, 2) = zero
                        DZxij(1:Ntz, 3) = (one - sbar**2)*dZodZ(1:Ntz, 1 - issym + 2, ii)
                    end if
                end if
            end if

        else

            if (innout == 0) signlss = -1
            if (innout == 1) signlss = +1

            if ((irz == 0 .and. issym == 0) .or. (irz == 1 .and. issym == 1)) then
                Dij(1:Ntz, 0) = (one + signlss*lss)*half*cosi(1:Ntz, ii)
                Dij(1:Ntz, 1) = (signlss)*half*cosi(1:Ntz, ii)
                Dij(1:Ntz, 2) = (one + signlss*lss)*half*sini(1:Ntz, ii)*(-im(ii))
                Dij(1:Ntz, 3) = (one + signlss*lss)*half*sini(1:Ntz, ii)*(+in(ii))
            else
                Dij(1:Ntz, 0) = (one + signlss*lss)*half*sini(1:Ntz, ii)
                Dij(1:Ntz, 1) = (signlss)*half*sini(1:Ntz, ii)
                Dij(1:Ntz, 2) = (one + signlss*lss)*half*cosi(1:Ntz, ii)*(+im(ii))
                Dij(1:Ntz, 3) = (one + signlss*lss)*half*cosi(1:Ntz, ii)*(-in(ii))
            end if

        end if

        if (Lcurvature == 5) then
            if (Igeometry == 3) then

                if (irz == 0) sg(1:Ntz, 1) = (Zij(1:Ntz, 1, 0)*Dij(1:Ntz, 2) - Dij(1:Ntz, 1)*Zij(1:Ntz, 2, 0))
                if (irz == 1) sg(1:Ntz, 1) = (Dij(1:Ntz, 1)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Dij(1:Ntz, 2))

            else
                if (Lcurvature == 5 .and. Igeometry /= 3) then
                    write (6, '("coords :      fatal : myid=",i3," ; Lcurvature.eq.5 .and. Igeometry.ne.3 ; Lcurvature.eq.5 can only be combined with Igeometry.ne.3;")') myid

                    stop "coords : Lcurvature.eq.5 .and. Igeometry.ne.3 : Lcurvature.eq.5 can only be combined with Igeometry.ne.3 ;"
                end if
            end if

        else

            select case (Igeometry)

            case (1)

                sg(1:Ntz, 1) = Dij(1:Ntz, 1)*rpol*rtor

                do ii = 1, 3
                    do jj = ii, 3
                        dguvij(1:Ntz, ii, jj) = Dij(1:Ntz, ii)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*Dij(1:Ntz, jj)
                    end do
                end do

            case (2)

                if (irz == 0) sg(1:Ntz, 1) = Dij(1:Ntz, 1)*Rij(1:Ntz, 0, 0) &
                                             + Rij(1:Ntz, 1, 0)*Dij(1:Ntz, 0)

                do ii = 1, 3
                    do jj = ii, 3
                        if (irz == 0) dguvij(1:Ntz, ii, jj) = Dij(1:Ntz, ii)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*Dij(1:Ntz, jj)
                        if (irz == 1) then
                            if (.true.) then
                                write (6, '("coords :      fatal : myid=",i3," ; .true. ; No Z-geometrical degree of freedom when Igeometry=2;")') myid

                                stop "coords : .true. : No Z-geometrical degree of freedom when Igeometry=2 ;"
                            end if
                        end if
                    end do
                end do

                dguvij(1:Ntz, 2, 2) = dguvij(1:Ntz, 2, 2) + two*Dij(1:Ntz, 0)*Rij(1:Ntz, 0, 0)

            case (3)

                if (LcoordinateSingularity) then
                    if (irz == 0) sg(1:Ntz, 1) = (Dij(1:Ntz, 0) + DRxij(1:Ntz, 0))*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, 0)) &
                                                 + Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 1, 0)*Dij(1:Ntz, 2) - (Dij(1:Ntz, 1) + DRxij(1:Ntz, 1))*Zij(1:Ntz, 2, 0)) &
                                                 + Rij(1:Ntz, 0, 0)*(DZxij(1:Ntz, 1)*Rij(1:Ntz, 2, 0))
                    if (irz == 1) sg(1:Ntz, 1) = DRxij(1:Ntz, 0)*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, 0)) &
                                                 + Rij(1:Ntz, 0, 0)*((Dij(1:Ntz, 1) + DZxij(1:Ntz, 1))*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Dij(1:Ntz, 2)) &
                                                 + Rij(1:Ntz, 0, 0)*(-DRxij(1:Ntz, 1)*Zij(1:Ntz, 2, 0))

                    do ii = 1, 3
                        do jj = ii, 3
                            if (irz == 0) dguvij(1:Ntz, ii, jj) = (Dij(1:Ntz, ii) + DRxij(1:Ntz, ii))*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*(Dij(1:Ntz, jj) + DRxij(1:Ntz, jj)) &
                                                                  + DZxij(1:Ntz, ii)*Zij(1:Ntz, jj, 0) + Zij(1:Ntz, ii, 0)*DZxij(1:Ntz, jj)
                            if (irz == 1) dguvij(1:Ntz, ii, jj) = (Dij(1:Ntz, ii) + DZxij(1:Ntz, ii))*Zij(1:Ntz, jj, 0) + Zij(1:Ntz, ii, 0)*(Dij(1:Ntz, jj) + DZxij(1:Ntz, jj)) &
                                                                  + DRxij(1:Ntz, ii)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*DRxij(1:Ntz, jj)
                        end do
                    end do

                    if (irz == 0) dguvij(1:Ntz, 3, 3) = dguvij(1:Ntz, 3, 3) + two*(Dij(1:Ntz, 0) + DRxij(1:Ntz, 0))*Rij(1:Ntz, 0, 0)
                    if (irz == 1) dguvij(1:Ntz, 3, 3) = dguvij(1:Ntz, 3, 3) + two*(DRxij(1:Ntz, 0))*Rij(1:Ntz, 0, 0)

                else

                    if (irz == 0) sg(1:Ntz, 1) = Dij(1:Ntz, 0)*(Zij(1:Ntz, 1, 0)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Zij(1:Ntz, 2, 0)) &
                                                 + Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 1, 0)*Dij(1:Ntz, 2) - Dij(1:Ntz, 1)*Zij(1:Ntz, 2, 0))
                    if (irz == 1) sg(1:Ntz, 1) = Rij(1:Ntz, 0, 0)*(Dij(1:Ntz, 1)*Rij(1:Ntz, 2, 0) - Rij(1:Ntz, 1, 0)*Dij(1:Ntz, 2))

                    do ii = 1, 3
                        do jj = ii, 3
                            if (irz == 0) dguvij(1:Ntz, ii, jj) = Dij(1:Ntz, ii)*Rij(1:Ntz, jj, 0) + Rij(1:Ntz, ii, 0)*Dij(1:Ntz, jj)
                            if (irz == 1) dguvij(1:Ntz, ii, jj) = Dij(1:Ntz, ii)*Zij(1:Ntz, jj, 0) + Zij(1:Ntz, ii, 0)*Dij(1:Ntz, jj)
                        end do
                    end do

                    if (irz == 0) dguvij(1:Ntz, 3, 3) = dguvij(1:Ntz, 3, 3) + two*Dij(1:Ntz, 0)*Rij(1:Ntz, 0, 0)
                end if

            case default

                if (.true.) then
                    write (6, '("coords :      fatal : myid=",i3," ; .true. ; supplied Igeometry is not yet supported for Lcurvature.eq.3 or Lcurvature.eq.4;")') myid

                    stop "coords : .true. : supplied Igeometry is not yet supported for Lcurvature.eq.3 or Lcurvature.eq.4 ;"
                end if

            end select

            do ii = 2, 3
                do jj = 1, ii - 1; dguvij(1:Ntz, ii, jj) = dguvij(1:Ntz, jj, ii)
                end do
            end do

            guvij(1:Ntz, 0, 0, 1) = zero

            if (Lcurvature == 3) then

                do ii = 1, 3
                do jj = 1, 3; guvij(1:Ntz, ii, jj, 1) = dguvij(1:Ntz, ii, jj) - guvij(1:Ntz, ii, jj, 0)*sg(1:Ntz, 1)/sg(1:Ntz, 0)
                end do
                end do

            else

                do ii = 1, 3
                do jj = 1, 3; guvij(1:Ntz, ii, jj, 1) = dguvij(1:Ntz, ii, jj)
                end do
                end do

            end if

        end if
    end select

end subroutine coords

