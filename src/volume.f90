
subroutine volume(lvol, vflag)

    use constants, only: zero, half, one, two, four, third, quart, pi2

    use numerical, only: vsmall, small

    use fileunits, only: ounit

    use inputlist, only: Wvolume, Igeometry, Nvol, pscale

    use allglobal, only: myid, cpus, &
                         YESstellsym, Mvol, &
                         Ntz, mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                         djkp, djkm, &
                         vvolume, dvolume, &
                         Rij, Zij, cosi, sini, &
                         dBdX, &
                         pi2nfp, pi2pi2nfp, pi2pi2nfpquart

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol

    integer :: vflag, Lcurvature

    integer :: jvol, ii, jj, kk, mi, ni, mj, nj, mk, nk, innout

    real(8) :: vol(0:1), vint(1:Ntz)

    real(8) :: Rei, Roi, Zei, Zoi, Rej, Roj, Zej, Zoj, Rek, Rok, Zek, Zok

    real(8) :: AA, BB, CC, DD, lss

    if (lvol > Nvol) then; vvolume(lvol) = one; dvolume = zero; return
    end if

    ; vol(0:1) = zero
    ; dvolume = zero

    do innout = 0, 1

        jvol = lvol - 1 + innout

        select case (Igeometry)

        case (1)

            vol(innout) = iRbc(1, jvol)

            if (dBdX%L .and. dBdX%innout == innout .and. dBdX%ii == 1) then
                if (dBdX%issym == 0) dvolume = one
            end if

        case (2)

            if (YESstellsym) then

                do ii = 1, mn; mi = im(ii); ni = in(ii)

                    do jj = 1, mn; mj = im(jj); nj = in(jj)

                        vol(innout) = vol(innout) + iRbc(ii, jvol)*iRbc(jj, jvol)*(djkp(ii, jj) + djkm(ii, jj))

                        if (dBdX%L .and. dBdX%innout == innout .and. dBdX%ii == ii) then
                            dvolume = dvolume + iRbc(jj, jvol)*(djkp(jj, ii) + djkm(jj, ii) + djkp(ii, jj) + djkm(ii, jj))
                        end if

                    end do

                end do

            else

                do ii = 1, mn; mi = im(ii); ni = in(ii)
                    do jj = 1, mn; mj = im(jj); nj = in(jj)

                        vol(innout) = vol(innout) + iRbc(ii, jvol)*iRbc(jj, jvol)*(djkp(ii, jj) + djkm(ii, jj)) &
                                      + iRbs(ii, jvol)*iRbs(jj, jvol)*(djkp(ii, jj) - djkm(ii, jj))

                        if (dBdX%L .and. dBdX%innout == innout .and. dBdX%ii == ii) then
                            if (dBdX%issym == 0) then
                                dvolume = dvolume + iRbc(jj, jvol)*(djkp(jj, ii) + djkm(jj, ii) + djkp(ii, jj) + djkm(ii, jj))
                            else
                                if (.true.) then
                                    write (6, '("volume :      fatal : myid=",i3," ; .true. ; derivatives of volume under construction;")') myid

                                    stop "volume : .true. : derivatives of volume under construction ;"
                                end if
                                dvolume = dvolume + iRbs(jj, jvol)*(djkp(jj, ii) - djkm(jj, ii) + djkp(ii, jj) - djkm(ii, jj))
                            end if
                        end if

                    end do
                end do

            end if

        case (3)

            if (lvol == 1 .and. innout == 0) then
                vol(1) = zero
                dvolume = zero
            else

                Lcurvature = 1

                lss = innout*two - one
                call coords(lvol, lss, Lcurvature, Ntz, mn)

                vint = Rij(1:Ntz, 0, 0)*(Zij(1:Ntz, 0, 0)*Rij(1:Ntz, 2, 0) - Zij(1:Ntz, 2, 0)*Rij(1:Ntz, 0, 0))
                vol(innout) = four*sum(vint)/float(Ntz)

                if (dBdX%L .and. dBdX%innout == innout) then

                    ii = dBdX%ii

                    if (dBdX%irz == 0) then

                        if (dBdX%issym == 0) then
                            vint = cosi(1:Ntz, ii)*(Zij(1:Ntz, 0, 0)*Rij(1:Ntz, 2, 0) - Zij(1:Ntz, 2, 0)*Rij(1:Ntz, 0, 0)) &
                                   + Rij(1:Ntz, 0, 0)*(-im(ii)*Zij(1:Ntz, 0, 0)*sini(1:Ntz, ii)) &
                                   + Rij(1:Ntz, 0, 0)*(-Zij(1:Ntz, 2, 0)*cosi(1:Ntz, ii))
                        else
                            vint = sini(1:Ntz, ii)*(Zij(1:Ntz, 0, 0)*Rij(1:Ntz, 2, 0) - Zij(1:Ntz, 2, 0)*Rij(1:Ntz, 0, 0)) &
                                   + Rij(1:Ntz, 0, 0)*(+im(ii)*Zij(1:Ntz, 0, 0)*cosi(1:Ntz, ii)) &
                                   + Rij(1:Ntz, 0, 0)*(-Zij(1:Ntz, 2, 0)*sini(1:Ntz, ii))
                        end if

                    else

                        if (dBdX%issym == 0) then
                            vint = Rij(1:Ntz, 0, 0)*(sini(1:Ntz, ii)*Rij(1:Ntz, 2, 0)) &
                                   + Rij(1:Ntz, 0, 0)*(-im(ii)*cosi(1:Ntz, ii)*Rij(1:Ntz, 0, 0))
                        else
                            vint = Rij(1:Ntz, 0, 0)*(cosi(1:Ntz, ii)*Rij(1:Ntz, 2, 0)) &
                                   + Rij(1:Ntz, 0, 0)*(+im(ii)*sini(1:Ntz, ii)*Rij(1:Ntz, 0, 0))
                        end if

                    end if

                    dvolume = four*sum(vint)/float(Ntz)

                else

                    dvolume = zero

                end if
            end if

        end select

    end do

    select case (Igeometry)
    case (1); vvolume(lvol) = (vol(1) - vol(0))*pi2pi2nfp; dvolume = dvolume*pi2pi2nfp
    case (2); vvolume(lvol) = (vol(1) - vol(0))*pi2pi2nfpquart; dvolume = dvolume*pi2pi2nfpquart
    case (3); vvolume(lvol) = (vol(1) - vol(0))*pi2pi2nfpquart*third; dvolume = dvolume*pi2pi2nfpquart*third
    case (4); vvolume(lvol) = one; dvolume = zero
        if (abs(pscale) > vsmall) then
            write (6, '("volume :      fatal : myid=",i3," ; abs(pscale).gt.vsmall ; need to compute volume;")') myid

            stop "volume : abs(pscale).gt.vsmall : need to compute volume ;"
        end if
    case default
        if (.true.) then
            write (6, '("volume :      fatal : myid=",i3," ; .true. ; invalid Igeometry;")') myid

            stop "volume : .true. : invalid Igeometry ;"
        end if
    end select

    if (dBdX%innout == 0) dvolume = -dvolume

    if (vflag == 0 .and. vvolume(lvol) < small) then
        write (6, '("volume :      fatal : myid=",i3," ; vflag.eq.0 .and. vvolume(lvol).lt.small ; volume cannot be zero or negative;")') myid

        stop "volume : vflag.eq.0 .and. vvolume(lvol).lt.small : volume cannot be zero or negative ;"
    end if

    if (vvolume(lvol) < small) then
        write (ounit, '("volume : ", 10x ," : myid=",i3," ; lvol=",i3," ; vvolume=",es13.5," ; volume cannot be zero or negative ;")') myid, lvol, vvolume(lvol)
        vvolume(lvol) = +9.9e+09
        vflag = 1
    else
        vflag = 0
    end if

end subroutine volume

