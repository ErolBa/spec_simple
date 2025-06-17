
subroutine bfield(zeta, st, Bst)

    use constants, only: zero, one, half, two

    use numerical, only: vsmall, small

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wbfield, Lrad, Mpol

    use cputiming, only: Tbfield

    use allglobal, only: myid, ncpu, cpus, &
                         mn, im, in, halfmm, regumm, &
                         ivol, gBzeta, Ate, Aze, Ato, Azo, &
                         NOTstellsym, &
                         Lcoordinatesingularity, Mvol, &
                         Node

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    real(8), intent(in) :: zeta, st(1:Node)
    real(8), intent(out) :: Bst(1:Node)

    integer :: lvol, ii, ll, mi, ni, ideriv
    real(8) :: teta, lss, sbar, sbarhm(0:1), arg, carg, sarg, dBu(1:3)
    real(8) :: cheby(0:Lrad(ivol), 0:1), zernike(0:Lrad(1), 0:Mpol, 0:1)

    real(8) :: TT(0:Lrad(ivol), 0:1)

    lvol = ivol; ideriv = 0

    Bst(1:Node) = (/zero, zero/)

    lss = st(1); teta = st(2)

    if (abs(lss) > one) return

    if (Lcoordinatesingularity) sbar = max((lss + one)*half, small)

    if (Lcoordinatesingularity) then
        call get_zernike(sbar, Lrad(lvol), Mpol, zernike(:, :, 0:1))
    else
        call get_cheby(lss, Lrad(lvol), cheby(0:Lrad(lvol), 0:1))
    end if

    dBu(1:3) = zero

    do ii = 1, mn; mi = im(ii); ni = in(ii); arg = mi*teta - ni*zeta; carg = cos(arg); sarg = sin(arg)

        if (Lcoordinatesingularity) then

            if (abs(sbar) < vsmall) then
                write (6, '("bfield :      fatal : myid=",i3," ; abs(sbar).lt.vsmall ; need to avoid divide-by-zero;")') myid

                stop "bfield : abs(sbar).lt.vsmall : need to avoid divide-by-zero ;"
            end if

            do ll = 0, Lrad(lvol); TT(ll, 0:1) = (/zernike(ll, mi, 0), zernike(ll, mi, 1)*half/)
            end do

        else

            do ll = 0, Lrad(lvol); TT(ll, 0:1) = (/cheby(ll, 0), cheby(ll, 1)/)
            end do

        end if

        do ll = 0, Lrad(lvol)
            ; dBu(1) = dBu(1) + (-mi*Aze(lvol, ideriv, ii)%s(ll) - ni*Ate(lvol, ideriv, ii)%s(ll))*TT(ll, 0)*sarg
            ; dBu(2) = dBu(2) + (-Aze(lvol, ideriv, ii)%s(ll))*TT(ll, 1)*carg
            ; dBu(3) = dBu(3) + (Ate(lvol, ideriv, ii)%s(ll))*TT(ll, 1)*carg
            if (NOTstellsym) then
                dBu(1) = dBu(1) + (+mi*Azo(lvol, ideriv, ii)%s(ll) + ni*Ato(lvol, ideriv, ii)%s(ll))*TT(ll, 0)*carg
                dBu(2) = dBu(2) + (-Azo(lvol, ideriv, ii)%s(ll))*TT(ll, 1)*sarg
                dBu(3) = dBu(3) + (Ato(lvol, ideriv, ii)%s(ll))*TT(ll, 1)*sarg
            end if
        end do

    end do

    gBzeta = dBu(3)

    if (abs(gBzeta) < vsmall) then

        cput = 0

        write (ounit, '("bfield : ",f10.2," : lvol=",i3," ; zeta="es23.15" ; (s,t)=("es23.15" ,"es23.15" ) ; B^z="es23.15" ;")') &
            cput - cpus, lvol, zeta, st(1:2), dBu(3)

        if (abs(dBu(3)) < vsmall) then
            write (6, '("bfield :      fatal : myid=",i3," ; abs(dBu(3)).lt.vsmall ; field is not toroidal;")') myid

            stop "bfield : abs(dBu(3)).lt.vsmall : field is not toroidal ;"
        end if
    end if

    Bst(1:2) = dBu(1:2)/gBzeta

end subroutine bfield
