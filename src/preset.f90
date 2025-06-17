
subroutine preset

    use constants, only: zero, one, mu0

    use numerical, only: sqrtmachprec, vsmall, small

    use fileunits, only: ounit

    use inputlist

    use cputiming, only: Tpreset

    use allglobal

    use fftw_interface

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer :: innout, idof, jk, ll, ii, ifail, ideriv, vvol, mi, ni, mj, nj, mk, nk, mimj, ninj, mkmj, nknj, jj, kk, lvol, mm, nn, imn
    integer :: lquad, igauleg, maxIquad, Mrad, jquad, Lcurvature, zerdof, iret, work1, work2
    real(8) :: teta, zeta, arg, lss, cszeta(0:1), error
    logical :: LComputeAxis

    logical :: Lchangeangle
    integer :: nb, ix, ij, ip, idx_mode
    real(8) :: xx

    select case (Istellsym)
    case (0); YESstellsym = .false.; NOTstellsym = .true.
    case (1); YESstellsym = .true.; NOTstellsym = .false.
    case default; 
        if (.true.) then
            write (6, '("readin :      fatal : myid=",i3," ; .true. ; illegal Istellsym;")') myid

            stop "readin : .true. : illegal Istellsym ;"
        end if
    end select

    Mvol = Nvol + Lfreebound

    if (allocated(beltramierror)) deallocate (beltramierror)
    allocate (beltramierror(1:Mvol, 1:9), stat=astat)
    beltramierror(1:Mvol, 1:9) = zero

    mn = 1 + Ntor + Mpol*(2*Ntor + 1)

    if (allocated(im)) deallocate (im)
    allocate (im(1:mn), stat=astat)
    im(1:mn) = 0
    if (allocated(in)) deallocate (in)
    allocate (in(1:mn), stat=astat)
    in(1:mn) = 0

    call gi00ab(Mpol, Ntor, Nfp, mn, im(1:mn), in(1:mn))

    if (allocated(halfmm)) deallocate (halfmm)
    allocate (halfmm(1:mn), stat=astat)
    halfmm(1:mn) = im(1:mn)*half
    if (allocated(regumm)) deallocate (regumm)
    allocate (regumm(1:mn), stat=astat)
    regumm(1:mn) = im(1:mn)*half

    if (Mregular >= 2) then

        where (im > Mregular) regumm = Mregular*half

    end if

    lMpol = 4*Mpol; lNtor = 4*Ntor

    mne = 1 + lNtor + lMpol*(2*lNtor + 1)

    if (allocated(ime)) deallocate (ime)
    allocate (ime(1:mne), stat=astat)
    ime(1:mne) = 0
    if (allocated(ine)) deallocate (ine)
    allocate (ine(1:mne), stat=astat)
    ine(1:mne) = 0

    call gi00ab(lMpol, lNtor, Nfp, mne, ime(1:mne), ine(1:mne))

    sMpol = iMpol; sNtor = iNtor

    if (iMpol <= 0) sMpol = Mpol - iMpol
    if (iNtor <= 0) sNtor = Ntor - iNtor
    if (Ntor == 0) sNtor = 0

    mns = 1 + sNtor + sMpol*(2*sNtor + 1)

    if (allocated(ims)) deallocate (ims)
    allocate (ims(1:mns), stat=astat)
    ims(1:mns) = 0
    if (allocated(ins)) deallocate (ins)
    allocate (ins(1:mns), stat=astat)
    ins(1:mns) = 0

    call gi00ab(sMpol, sNtor, Nfp, mns, ims(1:mns), ins(1:mns))

    if (Lcheck == 5) then; forcetol = 1.0e+12; nPpts = 0
    end if

    if (allocated(iRbc)) deallocate (iRbc)
    allocate (iRbc(1:mn, 0:Mvol), stat=astat)
    iRbc(1:mn, 0:Mvol) = zero
    if (allocated(iZbs)) deallocate (iZbs)
    allocate (iZbs(1:mn, 0:Mvol), stat=astat)
    iZbs(1:mn, 0:Mvol) = zero
    if (allocated(iRbs)) deallocate (iRbs)
    allocate (iRbs(1:mn, 0:Mvol), stat=astat)
    iRbs(1:mn, 0:Mvol) = zero
    if (allocated(iZbc)) deallocate (iZbc)
    allocate (iZbc(1:mn, 0:Mvol), stat=astat)
    iZbc(1:mn, 0:Mvol) = zero

    if (Lperturbed == 1) then
        if (allocated(dRbc)) deallocate (dRbc)
        allocate (dRbc(1:mn, 0:Mvol), stat=astat)
        dRbc(1:mn, 0:Mvol) = zero
        if (allocated(dZbs)) deallocate (dZbs)
        allocate (dZbs(1:mn, 0:Mvol), stat=astat)
        dZbs(1:mn, 0:Mvol) = zero
        if (allocated(dRbs)) deallocate (dRbs)
        allocate (dRbs(1:mn, 0:Mvol), stat=astat)
        dRbs(1:mn, 0:Mvol) = zero
        if (allocated(dZbc)) deallocate (dZbc)
        allocate (dZbc(1:mn, 0:Mvol), stat=astat)
        dZbc(1:mn, 0:Mvol) = zero
    end if

    if (allocated(iVns)) deallocate (iVns)
    allocate (iVns(1:mn), stat=astat)
    iVns(1:mn) = zero
    if (allocated(iBns)) deallocate (iBns)
    allocate (iBns(1:mn), stat=astat)
    iBns(1:mn) = zero
    if (allocated(iVnc)) deallocate (iVnc)
    allocate (iVnc(1:mn), stat=astat)
    iVnc(1:mn) = zero
    if (allocated(iBnc)) deallocate (iBnc)
    allocate (iBnc(1:mn), stat=astat)
    iBnc(1:mn) = zero

    if (allocated(ajk)) deallocate (ajk)
    allocate (ajk(1:mn), stat=astat)
    ajk(1:mn) = zero

    do kk = 1, mn; mk = im(kk); nk = in(kk)

        if (mk == 0) ajk(kk) = pi2

    end do

    if (myid == 0) then

        if (Igeometry == 3 .and. Rbc(0, +1) + Rbc(0, -1) > zero .and. Zbs(0, +1) - Zbs(0, -1) > zero) then; Lchangeangle = .true.
        else; Lchangeangle = .false.
        end if

        if (Lchangeangle) write (ounit, '("readin : " 10x " : CHANGING ANGLE ;")')

        do ii = 1, mn; mm = im(ii); nn = in(ii)/Nfp

            if (Lchangeangle) then; jj = -1; kk = -nn
            else; jj = +1; kk = +nn
            end if

            if (mm == 0 .and. nn == 0) then

                ; iRbc(ii, Nvol) = Rbc(nn, mm)
                ; iZbs(ii, Nvol) = zero
                if (NOTstellsym) then
                    ; iRbs(ii, Nvol) = zero
                    ; iZbc(ii, Nvol) = Zbc(nn, mm)
                else
                    ; iRbs(ii, Nvol) = zero
                    ; iZbc(ii, Nvol) = zero
                end if

            else

                ; iRbc(ii, Nvol) = Rbc(kk, mm) + Rbc(-kk, -mm)
                ; iZbs(ii, Nvol) = (Zbs(kk, mm) - Zbs(-kk, -mm))*jj
                if (NOTstellsym) then
                    ; iRbs(ii, Nvol) = (Rbs(kk, mm) - Rbs(-kk, -mm))*jj
                    ; iZbc(ii, Nvol) = Zbc(kk, mm) + Zbc(-kk, -mm)
                else
                    ; iRbs(ii, Nvol) = zero
                    ; iZbc(ii, Nvol) = zero
                end if

            end if

        end do

        select case (Linitialize)

        case (:0)

            if (Lchangeangle) then; jj = -1
            else; jj = +1
            end if

            do idx_mode = 1, num_modes
                mm = mmRZRZ(idx_mode)
                nn = nnRZRZ(idx_mode)

                do ii = 1, mn; mi = im(ii); ni = in(ii)
                    if (mm == 0 .and. mi == 0 .and. nn*Nfp == ni) then
                        iRbc(ii, 1:Nvol - 1) = allRZRZ(1, 1:Nvol - 1, idx_mode)
                        iZbs(ii, 1:Nvol - 1) = allRZRZ(2, 1:Nvol - 1, idx_mode)
                        if (NOTstellsym) then
                            iRbs(ii, 1:Nvol - 1) = allRZRZ(3, 1:Nvol - 1, idx_mode)
                            iZbc(ii, 1:Nvol - 1) = allRZRZ(4, 1:Nvol - 1, idx_mode)
                        else
                            iRbs(ii, 1:Nvol - 1) = zero
                            iZbc(ii, 1:Nvol - 1) = zero
                        end if
                    elseif (mm == mi .and. nn*Nfp == jj*ni) then
                        iRbc(ii, 1:Nvol - 1) = allRZRZ(1, 1:Nvol - 1, idx_mode)
                        iZbs(ii, 1:Nvol - 1) = jj*allRZRZ(2, 1:Nvol - 1, idx_mode)
                        if (NOTstellsym) then
                            iRbs(ii, 1:Nvol - 1) = jj*allRZRZ(3, 1:Nvol - 1, idx_mode)
                            iZbc(ii, 1:Nvol - 1) = allRZRZ(4, 1:Nvol - 1, idx_mode)
                        else
                            iRbs(ii, 1:Nvol - 1) = zero
                            iZbc(ii, 1:Nvol - 1) = zero
                        end if
                    end if
                end do

            end do

        end select

        if (Igeometry == 3) then
            if (Rac(0) > zero) then
                iRbc(1:Ntor + 1, 0) = Rac(0:Ntor)
                iZbs(1:Ntor + 1, 0) = Zas(0:Ntor)
                iRbs(1:Ntor + 1, 0) = Ras(0:Ntor)
                iZbc(1:Ntor + 1, 0) = Zac(0:Ntor)
            else
            end if
        end if

    end if

    if (Igeometry == 1 .or. Igeometry == 2) then
        ; iRbc(1:mn, 0) = zero
        if (NOTstellsym) then
            iRbs(1:mn, 0) = zero
        end if
    end if

    if (Igeometry == 3) then
        iZbs(1, 0:Mvol) = zero
    end if
    if (NOTstellsym) then
        iRbs(1, 0:Mvol) = zero
    end if

    if (Igeometry == 1 .and. Lreflect == 1) then
        iRbc(2:mn, Mvol) = iRbc(2:mn, Mvol)*half
        iRbc(2:mn, 0) = -iRbc(2:mn, Mvol)
        if (NOTstellsym) then
            iRbs(2:mn, Mvol) = iRbs(2:mn, Mvol)*half
            iRbs(2:mn, 0) = -iRbs(2:mn, Mvol)
        end if
    end if

    Rscale = iRbc(1, Mvol)

    if (myid == 0) write (ounit, '("readin : ", 10x ," : myid=",i3," ; Rscale=",es22.15," ;")') myid, Rscale

    call random_seed()

    pi2nfp = pi2/Nfp

    pi2pi2nfp = pi2*pi2nfp
    pi2pi2nfphalf = pi2*pi2nfp*half
    pi2pi2nfpquart = pi2*pi2nfp*quart

    Mrad = maxval(Lrad(1:Mvol))

    if (myid == 0) write (ounit, '("preset : ",10x," : myid=",i3," ; Mrad=",i3," : Lrad=",257(i3,",",:))') myid, Mrad, Lrad(1:Mvol)

    select case (Igeometry)
    case (1:2)
        if (YESstellsym) LGdof = mn
        if (NOTstellsym) LGdof = mn + mn - 1
    case (3)
        if (YESstellsym) LGdof = mn + mn - 1
        if (NOTstellsym) LGdof = mn + mn - 1 + mn - 1 + mn
    end select

    NGdof = (Mvol - 1)*LGdof

    do vvol = 0, Nvol

        if (ql(vvol) == 0 .and. qr(vvol) == 0) then; iota(vvol) = iota(vvol)
        else; iota(vvol) = (pl(vvol) + goldenmean*pr(vvol))/(ql(vvol) + goldenmean*qr(vvol))
        end if

        if (lq(vvol) == 0 .and. rq(vvol) == 0) then; oita(vvol) = oita(vvol)
        else; oita(vvol) = (lp(vvol) + goldenmean*rp(vvol))/(lq(vvol) + goldenmean*rq(vvol))
        end if

1002    format("preset : ", 10x, " :      ", 3x, " ; transform : ", i3, " : (", i3, " /", i3, " ) * (", i3, " /", i3, " ) = ", f18.15, " ; ", &
               "(", i3, " /", i3, " ) * (", i3, " /", i3, " ) = ", f18.15, " ; ")

    end do

    if (allocated(dtflux)) deallocate (dtflux)
    allocate (dtflux(1:Mvol), stat=astat)
    dtflux(1:Mvol) = zero
    if (allocated(dpflux)) deallocate (dpflux)
    allocate (dpflux(1:Mvol), stat=astat)
    dpflux(1:Mvol) = zero

    select case (Igeometry)
    case (1); dtflux(1) = tflux(1); dpflux(1) = pflux(1)
    case (2:3); dtflux(1) = tflux(1); dpflux(1) = zero
    end select

    dtflux(2:Mvol) = tflux(2:Mvol) - tflux(1:Mvol - 1)
    dpflux(2:Mvol) = pflux(2:Mvol) - pflux(1:Mvol - 1)

    dtflux(1:Mvol) = dtflux(1:Mvol)*phiedge/pi2
    dpflux(1:Mvol) = dpflux(1:Mvol)*phiedge/pi2

    if (Lconstraint == 3) then

        mu(1) = Ivolume(1)/(tflux(1)*phiedge)

        do vvol = 2, Mvol
            mu(vvol) = (Ivolume(vvol) - Ivolume(vvol - 1))/((tflux(vvol) - tflux(vvol - 1))*phiedge)
        end do

    end if

    if (allocated(sweight)) deallocate (sweight)
    allocate (sweight(1:Mvol), stat=astat)
    sweight(1:Mvol) = zero
    do vvol = 1, Mvol; sweight(vvol) = upsilon*(vvol*one/Nvol)**wpoloidal
    end do

    if (allocated(TT)) deallocate (TT)
    allocate (TT(0:Mrad, 0:1, 0:1), stat=astat)
    TT(0:Mrad, 0:1, 0:1) = zero
    if (allocated(RTT)) deallocate (RTT)
    allocate (RTT(0:Lrad(1), 0:Mpol, 0:1, 0:1), stat=astat)
    RTT(0:Lrad(1), 0:Mpol, 0:1, 0:1) = zero
    if (allocated(RTM)) deallocate (RTM)
    allocate (RTM(0:Lrad(1), 0:Mpol), stat=astat)
    RTM(0:Lrad(1), 0:Mpol) = zero

    call get_cheby(-one, Mrad, TT(:, 0, :))
    call get_cheby(one, Mrad, TT(:, 1, :))

    call get_zernike(zero, Lrad(1), Mpol, RTT(:, :, 0, :))
    call get_zernike(one, Lrad(1), Mpol, RTT(:, :, 1, :))
    call get_zernike_rm(zero, Lrad(1), Mpol, RTM(:, :))

    if (allocated(ImagneticOK)) deallocate (ImagneticOK)
    allocate (ImagneticOK(1:Mvol), stat=astat)
    ImagneticOK(1:Mvol) = .false.

    if (allocated(ki)) deallocate (ki)
    allocate (ki(1:mn, 0:1), stat=astat)
    ki(1:mn, 0:1) = 0
    if (allocated(kija)) deallocate (kija)
    allocate (kija(1:mn, 1:mn, 0:1), stat=astat)
    kija(1:mn, 1:mn, 0:1) = 0
    if (allocated(kijs)) deallocate (kijs)
    allocate (kijs(1:mn, 1:mn, 0:1), stat=astat)
    kijs(1:mn, 1:mn, 0:1) = 0

    do ii = 1, mn; mi = im(ii); ni = in(ii)

        call getimn(lMpol, lNtor, Nfp, mi, ni, kk)
        if (kk > 0) then
            if (mi == 0 .and. ni == 0) then; ki(ii, 0:1) = (/kk, 1/)
            else; ki(ii, 0:1) = (/kk, 2/)
            end if
        end if

        do jj = 1, mn; mj = im(jj); nj = in(jj); mimj = mi + mj; ninj = ni + nj

            call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
            if (kk > 0) then
                if (mimj == 0 .and. ninj == 0) then; kija(ii, jj, 0:1) = (/kk, 1/)
                else; kija(ii, jj, 0:1) = (/kk, 2/)
                end if
            end if
            ; ; mimj = mi - mj; ninj = ni - nj

            if (mimj > 0 .or. (mimj == 0 .and. ninj >= 0)) then
                call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
                if (kk > 0) then
                    if (mimj == 0 .and. ninj == 0) then; kijs(ii, jj, 0:1) = (/kk, 1/)
                    else; kijs(ii, jj, 0:1) = (/kk, 2/)
                    end if
                end if
            else
                call getimn(lMpol, lNtor, Nfp, -mimj, -ninj, kk)
                if (kk > 0) then
                    ; ; kijs(ii, jj, 0:1) = (/kk, -2/)
                end if
            end if

        end do

    end do

    if (Igeometry == 2) then

        if (allocated(djkp)) deallocate (djkp)
        allocate (djkp(1:mn, 1:mn), stat=astat)
        djkp(1:mn, 1:mn) = 0
        if (allocated(djkm)) deallocate (djkm)
        allocate (djkm(1:mn, 1:mn), stat=astat)
        djkm(1:mn, 1:mn) = 0

        do ii = 1, mn; mi = im(ii); ni = in(ii)
            do jj = 1, mn; mj = im(jj); nj = in(jj)
                if (mi - mj == 0 .and. ni - nj == 0) djkp(ii, jj) = 1
                if (mi + mj == 0 .and. ni + nj == 0) djkm(ii, jj) = 1
            end do
        end do

    end if

    if (allocated(iotakkii)) deallocate (iotakkii)
    allocate (iotakkii(1:mn), stat=astat)
    iotakkii(1:mn) = 0

    if (allocated(iotaksub)) deallocate (iotaksub)
    allocate (iotaksub(1:mn, 1:mns), stat=astat)
    iotaksub(1:mn, 1:mns) = 0
    if (allocated(iotaksgn)) deallocate (iotaksgn)
    allocate (iotaksgn(1:mn, 1:mns), stat=astat)
    iotaksgn(1:mn, 1:mns) = 0
    if (allocated(iotakadd)) deallocate (iotakadd)
    allocate (iotakadd(1:mn, 1:mns), stat=astat)
    iotakadd(1:mn, 1:mns) = 0

    do kk = 1, mn; mk = im(kk); nk = in(kk)

        call getimn(sMpol, sNtor, Nfp, mk, nk, ii)
        if (ii > 0) iotakkii(kk) = ii

        do jj = 1, mns; mj = ims(jj); nj = ins(jj)

            mkmj = mk - mj; nknj = nk - nj

            if (mkmj > 0 .or. (mkmj == 0 .and. nknj >= 0)) then

                call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
                if (ii > 0) then; iotaksub(kk, jj) = ii; iotaksgn(kk, jj) = 1
                end if

            else

                call getimn(sMpol, sNtor, Nfp, -mkmj, -nknj, ii)
                if (ii > 0) then; iotaksub(kk, jj) = ii; iotaksgn(kk, jj) = -1
                end if

            end if

            mkmj = mk + mj; nknj = nk + nj

            call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
            if (ii > 0) then; iotakadd(kk, jj) = ii
            end if

        end do

    end do

    if (allocated(IPDt)) deallocate (IPDt)
    allocate (IPDt(1:Mvol), stat=astat)
    IPDt(1:Mvol) = zero
    if (allocated(IPDtDpf)) deallocate (IPDtDpf)
    allocate (IPDtDpf(1:Mvol - 1, 1:Mvol - 1), stat=astat)
    IPDtDpf(1:Mvol - 1, 1:Mvol - 1) = zero

    if (allocated(cheby)) deallocate (cheby)
    allocate (cheby(0:Mrad, 0:2), stat=astat)
    cheby(0:Mrad, 0:2) = zero
    if (allocated(zernike)) deallocate (zernike)
    allocate (zernike(0:Lrad(1), 0:Mpol, 0:2), stat=astat)
    zernike(0:Lrad(1), 0:Mpol, 0:2) = zero

    if (allocated(Iquad)) deallocate (Iquad)
    allocate (Iquad(1:Mvol), stat=astat)
    Iquad(1:Mvol) = 0

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if

        if (Nquad > 0) then; Iquad(vvol) = Nquad
        else
            if (Lcoordinatesingularity) Iquad(vvol) = Mpol + 2*Lrad(vvol) - Nquad
            if (.not. Lcoordinatesingularity) Iquad(vvol) = 2*Lrad(vvol) - Nquad
        end if

    end do

    maxIquad = maxval(Iquad(1:Mvol))

    if (allocated(gaussianweight)) deallocate (gaussianweight)
    allocate (gaussianweight(1:maxIquad, 1:Mvol), stat=astat)
    gaussianweight(1:maxIquad, 1:Mvol) = zero
    if (allocated(gaussianabscissae)) deallocate (gaussianabscissae)
    allocate (gaussianabscissae(1:maxIquad, 1:Mvol), stat=astat)
    gaussianabscissae(1:maxIquad, 1:Mvol) = zero

    do vvol = 1, Mvol

        lquad = Iquad(vvol)

        call gauleg(lquad, gaussianweight(1:lquad, vvol), gaussianabscissae(1:lquad, vvol), igauleg)

    end do

    LBsequad = .false.
    LBnewton = .false.
    LBlinear = .false.

    if (LBeltrami == 1 .or. LBeltrami == 3 .or. LBeltrami == 5 .or. LBeltrami == 7) LBsequad = .true.
    if (LBeltrami == 2 .or. LBeltrami == 3 .or. LBeltrami == 6 .or. LBeltrami == 7) LBnewton = .true.
    if (LBeltrami == 4 .or. LBeltrami == 5 .or. LBeltrami == 6 .or. LBeltrami == 7) LBlinear = .true.

    if (Lconstraint == 2) then
        if (Lfreebound == 1) then
            write (6, '("preset :      fatal : myid=",i3," ; Lfreebound.eq.1 ; The combination of helicity constraint and free boundary is under construction;")') myid

            stop "preset : Lfreebound.eq.1 : The combination of helicity constraint and free boundary is under construction ;"
        end if
        if (Igeometry == 3 .and. myid == 0) then
            write (ounit, *) 'WARNING: The Hessian matrix needs further review for Igeometry = 3'
            write (ounit, *) '         However, it can still serve the purpose of Lfindzero = 2'
        end if
    end if

    if (myid == 0) then
        cput = 0
        write (ounit, '("preset : ",f10.2," : LBsequad="L2" , LBnewton="L2" , LBlinear="L2" ;")') cput - cpus, LBsequad, LBnewton, LBlinear
    end if

    if (allocated(BBweight)) deallocate (BBweight)
    allocate (BBweight(1:mn), stat=astat)
    BBweight(1:mn) = opsilon*exp(-escale*(im(1:mn)**2 + (in(1:mn)/Nfp)**2))

    if (myid == 0 .and. escale > small) then
        do ii = 1, mn; write (ounit, '("preset : " 10x " : myid="i3" ; ("i3","i3") : BBweight="es13.5" ;")') myid, im(ii), in(ii)/Nfp, BBweight(ii)
        end do
    end if

    if (allocated(mmpp)) deallocate (mmpp)
    allocate (mmpp(1:mn), stat=astat)
    mmpp(1:mn) = zero

    do ii = 1, mn; mi = im(ii)

        if (mi == 0) then; mmpp(ii) = zero
        else; mmpp(ii) = mi**pcondense
        end if

    end do

    if (allocated(NAdof)) deallocate (NAdof)
    allocate (NAdof(1:Mvol), stat=astat)
    NAdof(1:Mvol) = 0
    if (allocated(Nfielddof)) deallocate (Nfielddof)
    allocate (Nfielddof(1:Mvol), stat=astat)
    Nfielddof(1:Mvol) = 0
    if (allocated(NdMASmax)) deallocate (NdMASmax)
    allocate (NdMASmax(1:Mvol), stat=astat)
    NdMASmax(1:Mvol) = 0
    if (allocated(NdMAS)) deallocate (NdMAS)
    allocate (NdMAS(1:Mvol), stat=astat)
    NdMAS(1:Mvol) = 0

    if (allocated(Ate)) deallocate (Ate)
    allocate (Ate(1:Mvol, -2:2, 1:mn), stat=astat)
    if (allocated(Aze)) deallocate (Aze)
    allocate (Aze(1:Mvol, -2:2, 1:mn), stat=astat)
    if (allocated(Ato)) deallocate (Ato)
    allocate (Ato(1:Mvol, -2:2, 1:mn), stat=astat)
    if (allocated(Azo)) deallocate (Azo)
    allocate (Azo(1:Mvol, -2:2, 1:mn), stat=astat)

    if (allocated(Fso)) deallocate (Fso)
    allocate (Fso(1:Mvol, 1:mn), stat=astat)
    Fso(1:Mvol, 1:mn) = 0
    if (allocated(Fse)) deallocate (Fse)
    allocate (Fse(1:Mvol, 1:mn), stat=astat)
    Fse(1:Mvol, 1:mn) = 0

    if (allocated(Lma)) deallocate (Lma)
    allocate (Lma(1:Mvol, 1:mn), stat=astat)
    Lma(1:Mvol, 1:mn) = 0
    if (allocated(Lmb)) deallocate (Lmb)
    allocate (Lmb(1:Mvol, 1:mn), stat=astat)
    Lmb(1:Mvol, 1:mn) = 0
    if (allocated(Lmc)) deallocate (Lmc)
    allocate (Lmc(1:Mvol, 1:mn), stat=astat)
    Lmc(1:Mvol, 1:mn) = 0
    if (allocated(Lmd)) deallocate (Lmd)
    allocate (Lmd(1:Mvol, 1:mn), stat=astat)
    Lmd(1:Mvol, 1:mn) = 0
    if (allocated(Lme)) deallocate (Lme)
    allocate (Lme(1:Mvol, 1:mn), stat=astat)
    Lme(1:Mvol, 1:mn) = 0
    if (allocated(Lmf)) deallocate (Lmf)
    allocate (Lmf(1:Mvol, 1:mn), stat=astat)
    Lmf(1:Mvol, 1:mn) = 0
    if (allocated(Lmg)) deallocate (Lmg)
    allocate (Lmg(1:Mvol, 1:mn), stat=astat)
    Lmg(1:Mvol, 1:mn) = 0
    if (allocated(Lmh)) deallocate (Lmh)
    allocate (Lmh(1:Mvol, 1:mn), stat=astat)
    Lmh(1:Mvol, 1:mn) = 0

    if (allocated(Lmavalue)) deallocate (Lmavalue)
    allocate (Lmavalue(1:Mvol, 1:mn), stat=astat)
    Lmavalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmbvalue)) deallocate (Lmbvalue)
    allocate (Lmbvalue(1:Mvol, 1:mn), stat=astat)
    Lmbvalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmcvalue)) deallocate (Lmcvalue)
    allocate (Lmcvalue(1:Mvol, 1:mn), stat=astat)
    Lmcvalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmdvalue)) deallocate (Lmdvalue)
    allocate (Lmdvalue(1:Mvol, 1:mn), stat=astat)
    Lmdvalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmevalue)) deallocate (Lmevalue)
    allocate (Lmevalue(1:Mvol, 1:mn), stat=astat)
    Lmevalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmfvalue)) deallocate (Lmfvalue)
    allocate (Lmfvalue(1:Mvol, 1:mn), stat=astat)
    Lmfvalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmgvalue)) deallocate (Lmgvalue)
    allocate (Lmgvalue(1:Mvol, 1:mn), stat=astat)
    Lmgvalue(1:Mvol, 1:mn) = zero
    if (allocated(Lmhvalue)) deallocate (Lmhvalue)
    allocate (Lmhvalue(1:Mvol, 1:mn), stat=astat)
    Lmhvalue(1:Mvol, 1:mn) = zero

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if

        if (Lcoordinatesingularity) then
            zerdof = 0
            do ii = 2, Mpol
                do jj = ii, Lrad(vvol), 2
                    zerdof = zerdof + 2*ntor + 1
                    if (NOTstellsym) zerdof = zerdof + 2*ntor + 1
                end do
            end do
            zerdof = zerdof*2

            do jj = 0, Lrad(vvol), 2
                zerdof = zerdof + ntor + 1
                if (jj >= 2) zerdof = zerdof + ntor + 1

                if (NOTstellsym) then
                    zerdof = zerdof + ntor
                    if (jj >= 2) zerdof = zerdof + ntor
                end if
            end do

            if (Mpol >= 1) then
                do jj = 1, Lrad(vvol), 2
                    zerdof = zerdof + 2*ntor + 1
                    if (jj >= 2) zerdof = zerdof + 2*ntor + 1

                    if (NOTstellsym) then
                        zerdof = zerdof + 2*ntor + 1
                        if (jj >= 2) zerdof = zerdof + 2*ntor + 1
                    end if
                end do
            end if

            Nfielddof(vvol) = zerdof
            if (YESstellsym) NAdof(vvol) = zerdof + mn + Ntor + 1 + mn - 1 + 1 + 0
            if (NOTstellsym) NAdof(vvol) = zerdof + mn + mn - 1 + Ntor + 1 + Ntor + mn - 1 + mn - 1 + 1 + 0

            NAdof(vvol) = NAdof(vvol) - (ntor + 1)
            if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - ntor

            if (Mpol >= 1) then
                NAdof(vvol) = NAdof(vvol) - (2*ntor + 1)
                if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - (2*ntor + 1)
            end if

        else
            if (YESstellsym) NAdof(vvol) = 2*(mn)*(Lrad(vvol)) + mn - 1 + 1 + 1
            if (NOTstellsym) NAdof(vvol) = 2*(mn + mn - 1)*(Lrad(vvol)) + mn - 1 + mn - 1 + 1 + 1

            if (YESstellsym) Nfielddof(vvol) = 2*(mn)*(Lrad(vvol))
            if (NOTstellsym) Nfielddof(vvol) = 2*(mn + mn - 1)*(Lrad(vvol))

        end if

        do ii = 1, mn

            do ideriv = -2, 2

                if (allocated(Ate(vvol, ideriv, ii)%s)) deallocate (Ate(vvol, ideriv, ii)%s)
                allocate (Ate(vvol, ideriv, ii)%s(0:Lrad(vvol)), stat=astat)
                Ate(vvol, ideriv, ii)%s(0:Lrad(vvol)) = zero
                if (allocated(Aze(vvol, ideriv, ii)%s)) deallocate (Aze(vvol, ideriv, ii)%s)
                allocate (Aze(vvol, ideriv, ii)%s(0:Lrad(vvol)), stat=astat)
                Aze(vvol, ideriv, ii)%s(0:Lrad(vvol)) = zero
                if (allocated(Ato(vvol, ideriv, ii)%s)) deallocate (Ato(vvol, ideriv, ii)%s)
                allocate (Ato(vvol, ideriv, ii)%s(0:Lrad(vvol)), stat=astat)
                Ato(vvol, ideriv, ii)%s(0:Lrad(vvol)) = zero
                if (allocated(Azo(vvol, ideriv, ii)%s)) deallocate (Azo(vvol, ideriv, ii)%s)
                allocate (Azo(vvol, ideriv, ii)%s(0:Lrad(vvol)), stat=astat)
                Azo(vvol, ideriv, ii)%s(0:Lrad(vvol)) = zero

            end do

            ; ideriv = 0

            if (allocated(Ate(vvol, ideriv, ii)%i)) deallocate (Ate(vvol, ideriv, ii)%i)
            allocate (Ate(vvol, ideriv, ii)%i(0:Lrad(vvol)), stat=astat)
            Ate(vvol, ideriv, ii)%i(0:Lrad(vvol)) = 0
            if (allocated(Aze(vvol, ideriv, ii)%i)) deallocate (Aze(vvol, ideriv, ii)%i)
            allocate (Aze(vvol, ideriv, ii)%i(0:Lrad(vvol)), stat=astat)
            Aze(vvol, ideriv, ii)%i(0:Lrad(vvol)) = 0
            if (allocated(Ato(vvol, ideriv, ii)%i)) deallocate (Ato(vvol, ideriv, ii)%i)
            allocate (Ato(vvol, ideriv, ii)%i(0:Lrad(vvol)), stat=astat)
            Ato(vvol, ideriv, ii)%i(0:Lrad(vvol)) = 0
            if (allocated(Azo(vvol, ideriv, ii)%i)) deallocate (Azo(vvol, ideriv, ii)%i)
            allocate (Azo(vvol, ideriv, ii)%i(0:Lrad(vvol)), stat=astat)
            Azo(vvol, ideriv, ii)%i(0:Lrad(vvol)) = 0

        end do

        select case (Linitgues)
        case (0); 
        case (1); Ate(vvol, 0, 1)%s(0:1) = dtflux(vvol)*half
            ; ; Aze(vvol, 0, 1)%s(0:1) = dpflux(vvol)*half
            if (Lcoordinatesingularity) then
                ; ; Ate(vvol, 0, 1)%s(2) = dtflux(vvol)*half*half
            end if
        case (2); 
        case (3); 
            do ii = 1, mn

                do ideriv = -2, 2

                    call random_number(Ate(vvol, ideriv, ii)%s)
                    call random_number(Aze(vvol, ideriv, ii)%s)
                    Ate(vvol, ideriv, ii)%s = Ate(vvol, ideriv, ii)%s*maxrndgues
                    Aze(vvol, ideriv, ii)%s = Aze(vvol, ideriv, ii)%s*maxrndgues
                    if (.not. YESstellsym) then
                        call random_number(Ato(vvol, ideriv, ii)%s)
                        call random_number(Azo(vvol, ideriv, ii)%s)
                        Ato(vvol, ideriv, ii)%s = Ato(vvol, ideriv, ii)%s*maxrndgues
                        Azo(vvol, ideriv, ii)%s = Azo(vvol, ideriv, ii)%s*maxrndgues
                    end if

                end do

            end do

        end select

        idof = 0

        if (Lcoordinatesingularity) then

            do ii = 1, mn; mi = im(ii); ni = in(ii)

                do ll = 0, Lrad(vvol)
                    if (ll >= mi .and. mod(mi + ll, 2) == 0) then
                    if (.not. ((ll == 0 .and. mi == 0) .or. (ll == 1 .and. mi == 1))) then
                        ; idof = idof + 1; Ate(vvol, 0, ii)%i(ll) = idof
                    end if
                    ; ; idof = idof + 1; Aze(vvol, 0, ii)%i(ll) = idof
                    if (NOTstellsym .and. ii > 1) then
                        if (.not. ((ll == 0 .and. mi == 0) .or. (ll == 1 .and. mi == 1))) then
                            ; idof = idof + 1; Ato(vvol, 0, ii)%i(ll) = idof
                        end if
                        ; ; idof = idof + 1; Azo(vvol, 0, ii)%i(ll) = idof
                    end if
                    end if
                end do

            end do

            do ii = 1, mn; mi = im(ii); ni = in(ii)
                if (mi /= 0 .and. mi /= 1) then; idof = idof + 1; Lma(vvol, ii) = idof
                end if
                if (mi == 0) then; idof = idof + 1; Lmb(vvol, ii) = idof
                end if
                if (ii > 1) then; idof = idof + 1; Lme(vvol, ii) = idof
                end if
                if (ii == 1) then; idof = idof + 1; Lmg(vvol, ii) = idof
                end if
                if (NOTstellsym) then
                    if (mi /= 0 .and. mi /= 1) then; idof = idof + 1; Lmc(vvol, ii) = idof
                    end if
                    if (ii > 1) then; idof = idof + 1; Lmf(vvol, ii) = idof
                    end if
                    if (ii > 1 .and. mi == 0) then; idof = idof + 1; Lmd(vvol, ii) = idof
                    end if
                end if

            end do

            if (idof /= NAdof(vvol)) then
                write (6, '("preset :      fatal : myid=",i3," ; idof.ne.NAdof(vvol) ; need to count Beltrami degrees-of-freedom more carefully  for coordinate singularity;")') myid

                stop "preset : idof.ne.NAdof(vvol) : need to count Beltrami degrees-of-freedom more carefully  for coordinate singularity ;"
            end if

            if ((idof + 1)**2 >= huge(idof)) then
                write (6, '("preset :      fatal : myid=",i3," ; (idof+1)**2.ge.HUGE(idof)) ; NAdof too big, should be smaller than maximum of int32 type;")') myid

                stop "preset : (idof+1)**2.ge.HUGE(idof)) : NAdof too big, should be smaller than maximum of int32 type ;"
            end if

        else

            do ii = 1, mn
                do ll = 1, Lrad(vvol); idof = idof + 1; Ate(vvol, 0, ii)%i(ll) = idof
                    ; ; idof = idof + 1; Aze(vvol, 0, ii)%i(ll) = idof
                    if (ii > 1 .and. NOTstellsym) then; idof = idof + 1; Ato(vvol, 0, ii)%i(ll) = idof
                        ; ; idof = idof + 1; Azo(vvol, 0, ii)%i(ll) = idof
                    end if
                end do
            end do

            do ii = 1, mn
                if (ii > 1) then; idof = idof + 1; Lme(vvol, ii) = idof
                end if
                if (ii > 1 .and. NOTstellsym) then; idof = idof + 1; Lmf(vvol, ii) = idof
                end if
                if (ii == 1) then; idof = idof + 1; Lmg(vvol, ii) = idof
                    ; ; idof = idof + 1; Lmh(vvol, ii) = idof
                end if
            end do

            if (idof /= NAdof(vvol)) then
                write (6, '("preset :      fatal : myid=",i3," ; idof.ne.NAdof(vvol) ; need to count degrees-of-freedom more carefully for new matrix;")') myid

                stop "preset : idof.ne.NAdof(vvol) : need to count degrees-of-freedom more carefully for new matrix ;"
            end if

            if ((idof + 1)**2 >= huge(idof)) then
                write (6, '("preset :      fatal : myid=",i3," ; (idof+1)**2.ge.HUGE(idof)) ; NAdof too big, should be smaller than maximum of int32 type;")') myid

                stop "preset : (idof+1)**2.ge.HUGE(idof)) : NAdof too big, should be smaller than maximum of int32 type ;"
            end if

        end if

        if (idof /= NAdof(vvol)) then
            write (6, '("preset :      fatal : myid=",i3," ; idof.ne.NAdof(vvol) ; impossible logic;")') myid

            stop "preset : idof.ne.NAdof(vvol) : impossible logic ;"
        end if

        do ii = 1, mn
            do jj = 0, Lrad(vvol)
                if (Ate(vvol, 0, ii)%i(jj) == 0) Ate(vvol, 0, ii)%s(jj) = zero
                if (Aze(vvol, 0, ii)%i(jj) == 0) Aze(vvol, 0, ii)%s(jj) = zero
                if (.not. YESstellsym) then
                    if (Ato(vvol, 0, ii)%i(jj) == 0) Azo(vvol, 0, ii)%s(jj) = zero
                    if (Azo(vvol, 0, ii)%i(jj) == 0) Azo(vvol, 0, ii)%s(jj) = zero
                end if
            end do
        end do

    end do

    if (Linitgues == 2) then; call ra00aa('R')
    end if

    if (myid == 0) then
        cput = 0
        write (ounit, '("preset : ", 10x ," : ")')
        write (ounit, '("preset : ",f10.2," : Nquad="i4" ; mn="i5" ; NGdof="i6" ; NAdof="16(i6",")" ...")') cput - cpus, Nquad, mn, NGdof, NAdof(1:min(Mvol, 16))
    end if

    Nt = max(Ndiscrete*4*Mpol, 1); Nz = max(Ndiscrete*4*Ntor, 1); Ntz = Nt*Nz; soNtz = one/sqrt(one*Ntz)

    ; ; hNt = Nt/2
    if (Nz > 1) then; hNz = Nz/2
    else; hNz = 0
    end if

    if (myid == 0) then
        cput = 0
        write (ounit, '("preset : ", 10x ," : ")')
        write (ounit, '("preset : ",f10.2," : Nt="i6" ; Nz="i6" ; Ntz="i9" ;")') cput - cpus, Nt, Nz, Ntz
    end if

    if (allocated(iRij)) deallocate (iRij)
    allocate (iRij(1:Ntz, 0:Mvol), stat=astat)
    iRij(1:Ntz, 0:Mvol) = zero
    if (allocated(iZij)) deallocate (iZij)
    allocate (iZij(1:Ntz, 0:Mvol), stat=astat)
    iZij(1:Ntz, 0:Mvol) = zero
    if (allocated(dRij)) deallocate (dRij)
    allocate (dRij(1:Ntz, 1:Mvol), stat=astat)
    dRij(1:Ntz, 1:Mvol) = zero
    if (allocated(dZij)) deallocate (dZij)
    allocate (dZij(1:Ntz, 1:Mvol), stat=astat)
    dZij(1:Ntz, 1:Mvol) = zero
    if (allocated(tRij)) deallocate (tRij)
    allocate (tRij(1:Ntz, 0:Mvol), stat=astat)
    tRij(1:Ntz, 0:Mvol) = zero
    if (allocated(tZij)) deallocate (tZij)
    allocate (tZij(1:Ntz, 0:Mvol), stat=astat)
    tZij(1:Ntz, 0:Mvol) = zero

    if (allocated(Rij)) deallocate (Rij)
    allocate (Rij(1:Ntz, 0:3, 0:3), stat=astat)
    Rij(1:Ntz, 0:3, 0:3) = zero
    if (allocated(Zij)) deallocate (Zij)
    allocate (Zij(1:Ntz, 0:3, 0:3), stat=astat)
    Zij(1:Ntz, 0:3, 0:3) = zero
    if (allocated(sg)) deallocate (sg)
    allocate (sg(1:Ntz, 0:3), stat=astat)
    sg(1:Ntz, 0:3) = zero
    if (allocated(guvij)) deallocate (guvij)
    allocate (guvij(1:Ntz, 0:3, 0:3, -1:3), stat=astat)
    guvij(1:Ntz, 0:3, 0:3, -1:3) = zero
    if (allocated(gvuij)) deallocate (gvuij)
    allocate (gvuij(1:Ntz, 0:3, 0:3), stat=astat)
    gvuij(1:Ntz, 0:3, 0:3) = zero

    if ((Lfindzero == 2) .or. (Lcheck == 5 .or. LHevalues .or. LHevectors .or. LHmatrix .or. Lperturbed == 1)) then
        if (allocated(dRadR)) deallocate (dRadR)
        allocate (dRadR(1:mn, 0:1, 0:1, 1:mn), stat=astat)
        dRadR(1:mn, 0:1, 0:1, 1:mn) = zero
        if (allocated(dRadZ)) deallocate (dRadZ)
        allocate (dRadZ(1:mn, 0:1, 0:1, 1:mn), stat=astat)
        dRadZ(1:mn, 0:1, 0:1, 1:mn) = zero
        if (allocated(dZadR)) deallocate (dZadR)
        allocate (dZadR(1:mn, 0:1, 0:1, 1:mn), stat=astat)
        dZadR(1:mn, 0:1, 0:1, 1:mn) = zero
        if (allocated(dZadZ)) deallocate (dZadZ)
        allocate (dZadZ(1:mn, 0:1, 0:1, 1:mn), stat=astat)
        dZadZ(1:mn, 0:1, 0:1, 1:mn) = zero

        if (allocated(dRodR)) deallocate (dRodR)
        allocate (dRodR(1:Ntz, 0:3, 1:mn), stat=astat)
        dRodR(1:Ntz, 0:3, 1:mn) = zero
        if (allocated(dRodZ)) deallocate (dRodZ)
        allocate (dRodZ(1:Ntz, 0:3, 1:mn), stat=astat)
        dRodZ(1:Ntz, 0:3, 1:mn) = zero
        if (allocated(dZodR)) deallocate (dZodR)
        allocate (dZodR(1:Ntz, 0:3, 1:mn), stat=astat)
        dZodR(1:Ntz, 0:3, 1:mn) = zero
        if (allocated(dZodZ)) deallocate (dZodZ)
        allocate (dZodZ(1:Ntz, 0:3, 1:mn), stat=astat)
        dZodZ(1:Ntz, 0:3, 1:mn) = zero
    end if

    if (allocated(goomne)) deallocate (goomne)
    allocate (goomne(0:mne, maxIquad), stat=astat)
    goomne(0:mne, maxIquad) = zero
    if (allocated(goomno)) deallocate (goomno)
    allocate (goomno(0:mne, maxIquad), stat=astat)
    goomno(0:mne, maxIquad) = zero
    if (allocated(gssmne)) deallocate (gssmne)
    allocate (gssmne(0:mne, maxIquad), stat=astat)
    gssmne(0:mne, maxIquad) = zero
    if (allocated(gssmno)) deallocate (gssmno)
    allocate (gssmno(0:mne, maxIquad), stat=astat)
    gssmno(0:mne, maxIquad) = zero
    if (allocated(gstmne)) deallocate (gstmne)
    allocate (gstmne(0:mne, maxIquad), stat=astat)
    gstmne(0:mne, maxIquad) = zero
    if (allocated(gstmno)) deallocate (gstmno)
    allocate (gstmno(0:mne, maxIquad), stat=astat)
    gstmno(0:mne, maxIquad) = zero
    if (allocated(gszmne)) deallocate (gszmne)
    allocate (gszmne(0:mne, maxIquad), stat=astat)
    gszmne(0:mne, maxIquad) = zero
    if (allocated(gszmno)) deallocate (gszmno)
    allocate (gszmno(0:mne, maxIquad), stat=astat)
    gszmno(0:mne, maxIquad) = zero
    if (allocated(gttmne)) deallocate (gttmne)
    allocate (gttmne(0:mne, maxIquad), stat=astat)
    gttmne(0:mne, maxIquad) = zero
    if (allocated(gttmno)) deallocate (gttmno)
    allocate (gttmno(0:mne, maxIquad), stat=astat)
    gttmno(0:mne, maxIquad) = zero
    if (allocated(gtzmne)) deallocate (gtzmne)
    allocate (gtzmne(0:mne, maxIquad), stat=astat)
    gtzmne(0:mne, maxIquad) = zero
    if (allocated(gtzmno)) deallocate (gtzmno)
    allocate (gtzmno(0:mne, maxIquad), stat=astat)
    gtzmno(0:mne, maxIquad) = zero
    if (allocated(gzzmne)) deallocate (gzzmne)
    allocate (gzzmne(0:mne, maxIquad), stat=astat)
    gzzmne(0:mne, maxIquad) = zero
    if (allocated(gzzmno)) deallocate (gzzmno)
    allocate (gzzmno(0:mne, maxIquad), stat=astat)
    gzzmno(0:mne, maxIquad) = zero

    if (allocated(ijreal)) deallocate (ijreal)
    allocate (ijreal(1:Ntz), stat=astat)
    ijreal(1:Ntz) = zero
    if (allocated(ijimag)) deallocate (ijimag)
    allocate (ijimag(1:Ntz), stat=astat)
    ijimag(1:Ntz) = zero
    if (allocated(jireal)) deallocate (jireal)
    allocate (jireal(1:Ntz), stat=astat)
    jireal(1:Ntz) = zero
    if (allocated(jiimag)) deallocate (jiimag)
    allocate (jiimag(1:Ntz), stat=astat)
    jiimag(1:Ntz) = zero

    if (allocated(jkreal)) deallocate (jkreal)
    allocate (jkreal(1:Ntz), stat=astat)
    jkreal(1:Ntz) = zero
    if (allocated(jkimag)) deallocate (jkimag)
    allocate (jkimag(1:Ntz), stat=astat)
    jkimag(1:Ntz) = zero
    if (allocated(kjreal)) deallocate (kjreal)
    allocate (kjreal(1:Ntz), stat=astat)
    kjreal(1:Ntz) = zero
    if (allocated(kjimag)) deallocate (kjimag)
    allocate (kjimag(1:Ntz), stat=astat)
    kjimag(1:Ntz) = zero

    if (allocated(cplxin)) deallocate (cplxin)
    allocate (cplxin(1:Nt, 1:Nz, 1), stat=astat)
    cplxin(1:Nt, 1:Nz, 1) = zero
    if (allocated(cplxout)) deallocate (cplxout)
    allocate (cplxout(1:Nt, 1:Nz, 1), stat=astat)
    cplxout(1:Nt, 1:Nz, 1) = zero

    planf = fftw_plan_dft_2d(Nz, Nt, cplxin(:, :, 1), cplxout(:, :, 1), FFTW_FORWARD, FFTW_MEASURE + FFTW_DESTROY_INPUT)
    planb = fftw_plan_dft_2d(Nz, Nt, cplxin(:, :, 1), cplxout(:, :, 1), FFTW_BACKWARD, FFTW_MEASURE + FFTW_DESTROY_INPUT)

    if (allocated(efmn)) deallocate (efmn)
    allocate (efmn(1:mne), stat=astat)
    efmn(1:mne) = zero
    if (allocated(ofmn)) deallocate (ofmn)
    allocate (ofmn(1:mne), stat=astat)
    ofmn(1:mne) = zero
    if (allocated(cfmn)) deallocate (cfmn)
    allocate (cfmn(1:mne), stat=astat)
    cfmn(1:mne) = zero
    if (allocated(sfmn)) deallocate (sfmn)
    allocate (sfmn(1:mne), stat=astat)
    sfmn(1:mne) = zero
    if (allocated(evmn)) deallocate (evmn)
    allocate (evmn(1:mne), stat=astat)
    evmn(1:mne) = zero
    if (allocated(odmn)) deallocate (odmn)
    allocate (odmn(1:mne), stat=astat)
    odmn(1:mne) = zero
    if (allocated(comn)) deallocate (comn)
    allocate (comn(1:mne), stat=astat)
    comn(1:mne) = zero
    if (allocated(simn)) deallocate (simn)
    allocate (simn(1:mne), stat=astat)
    simn(1:mne) = zero

    if (allocated(gteta)) deallocate (gteta)
    allocate (gteta(1:Ntz), stat=astat)
    gteta(1:Ntz) = zero
    if (allocated(gzeta)) deallocate (gzeta)
    allocate (gzeta(1:Ntz), stat=astat)
    gzeta(1:Ntz) = zero

    if (allocated(cosi)) deallocate (cosi)
    allocate (cosi(1:Ntz, 1:mn), stat=astat)
    cosi(1:Ntz, 1:mn) = zero
    if (allocated(sini)) deallocate (sini)
    allocate (sini(1:Ntz, 1:mn), stat=astat)
    sini(1:Ntz, 1:mn) = zero

    if (Nz == 0) then
        write (6, '("preset :      fatal : myid=",i3," ; Nz.eq.0 ; illegal division;")') myid

        stop "preset : Nz.eq.0 : illegal division ;"
    end if
    if (Nt == 0) then
        write (6, '("preset :      fatal : myid=",i3," ; Nt.eq.0 ; illegal division;")') myid

        stop "preset : Nt.eq.0 : illegal division ;"
    end if

    do ii = 1, mn; mi = im(ii); ni = in(ii)

        do kk = 0, Nz - 1; zeta = kk*pi2nfp/Nz
            do jj = 0, Nt - 1; teta = jj*pi2/Nt; jk = 1 + jj + kk*Nt; arg = mi*teta - ni*zeta
                gteta(jk) = teta
                gzeta(jk) = zeta
                cosi(jk, ii) = cos(arg)
                sini(jk, ii) = sin(arg)
            end do
        end do

    end do

    if (Igeometry == 3 .and. iRbc(1, 0) < small) then

        write (*, *) "Finding initial axis..."
        select case (Linitialize)
        case (:-1); vvol = Nvol + Linitialize
        case (0); vvol = 1
        case (1); vvol = Nvol
        case (2); vvol = Mvol
        end select

        call rzaxis(Mvol, mn, iRbc(1:mn, 0:Mvol), iZbs(1:mn, 0:Mvol), iRbs(1:mn, 0:Mvol), iZbc(1:mn, 0:Mvol), vvol, .false.)

    end if

    if (allocated(psifactor)) deallocate (psifactor)
    allocate (psifactor(1:mn, 1:Mvol), stat=astat)
    psifactor(1:mn, 1:Mvol) = zero
    if (allocated(inifactor)) deallocate (inifactor)
    allocate (inifactor(1:mn, 1:Mvol), stat=astat)
    inifactor(1:mn, 1:Mvol) = zero

    psifactor(1:mn, 1:Mvol) = one
    inifactor(1:mn, 1:Mvol) = one

    select case (Igeometry)

    case (1)

        psifactor(1:mn, 1:Nvol) = one

    case (2)

        do vvol = 1, Nvol
            do ii = 1, mn
                if (im(ii) == 0) then; psifactor(ii, vvol) = tflux(vvol)**(+half)
                else; psifactor(ii, vvol) = tflux(vvol)**(halfmm(ii) - half)
                end if
            end do
        end do

    case (3)

        do vvol = 1, Nvol
            do ii = 1, mn
                if (im(ii) == 0) then; psifactor(ii, vvol) = Rscale*tflux(vvol)**zero
                    ; inifactor(ii, vvol) = Rscale*tflux(vvol)**half
                else; psifactor(ii, vvol) = Rscale*tflux(vvol)**halfmm(ii)
                    ; inifactor(ii, vvol) = Rscale*tflux(vvol)**halfmm(ii)
                end if
            end do
        end do

    case default

        if (.true.) then
            write (6, '("readin :      fatal : myid=",i3," ; .true. ; invalid Igeometry for construction of psifactor;")') myid

            stop "readin : .true. : invalid Igeometry for construction of psifactor ;"
        end if

    end select

    if (Linitialize /= 0) then

        select case (Igeometry)

        case (1)

            do vvol = 1, Nvol
                ; iRbc(1:mn, vvol) = iRbc(1:mn, Mvol)*tflux(vvol)/tflux(Mvol)
                if (NOTstellsym) then
                    iRbs(2:mn, vvol) = iRbs(2:mn, Mvol)*tflux(vvol)/tflux(Mvol)
                end if
            end do

        case (2)

            if (Linitialize /= 1) then
                write (6, '("preset :      fatal : myid=",i3," ; Linitialize.ne.1 ; geometrical initialization under construction for cylindrical;")') myid

                stop "preset : Linitialize.ne.1 : geometrical initialization under construction for cylindrical ;"
            end if

            do vvol = 1, Nvol - 1
                ; iRbc(1:mn, vvol) = iRbc(1:mn, Nvol)*psifactor(1:mn, vvol)
                if (NOTstellsym) then
                    iRbs(2:mn, vvol) = iRbs(2:mn, Nvol)*psifactor(2:mn, vvol)
                end if
            end do

        case (3)

            if (Linitialize < 0) then
                write (6, '("preset :      fatal : myid=",i3," ; Linitialize.lt.0 ; geometrical initialization under construction for toroidal;")') myid

                stop "preset : Linitialize.lt.0 : geometrical initialization under construction for toroidal ;"
            end if

            lvol = Nvol - 1 + Linitialize

            do vvol = 1, lvol - 1
                ; iRbc(1:mn, vvol) = iRbc(1:mn, 0) + (iRbc(1:mn, lvol) - iRbc(1:mn, 0))*(inifactor(1:mn, vvol)/Rscale)/tflux(lvol)**halfmm(1:mn)
                ; iZbs(2:mn, vvol) = iZbs(2:mn, 0) + (iZbs(2:mn, lvol) - iZbs(2:mn, 0))*(inifactor(2:mn, vvol)/Rscale)/tflux(lvol)**halfmm(2:mn)
                if (NOTstellsym) then
                    iRbs(2:mn, vvol) = iRbs(2:mn, 0) + (iRbs(2:mn, lvol) - iRbs(2:mn, 0))*(inifactor(2:mn, vvol)/Rscale)/tflux(lvol)**halfmm(2:mn)
                    iZbc(1:mn, vvol) = iZbc(1:mn, 0) + (iZbc(1:mn, lvol) - iZbc(1:mn, 0))*(inifactor(1:mn, vvol)/Rscale)/tflux(lvol)**halfmm(1:mn)
                end if
            end do

        end select

    end if

    if (allocated(Bsupumn)) deallocate (Bsupumn)
    allocate (Bsupumn(1:Nvol, 0:1, 1:mn), stat=astat)
    Bsupumn(1:Nvol, 0:1, 1:mn) = zero
    if (allocated(Bsupvmn)) deallocate (Bsupvmn)
    allocate (Bsupvmn(1:Nvol, 0:1, 1:mn), stat=astat)
    Bsupvmn(1:Nvol, 0:1, 1:mn) = zero

    if (allocated(diotadxup)) deallocate (diotadxup)
    allocate (diotadxup(0:1, -1:2, 1:Mvol), stat=astat)
    diotadxup(0:1, -1:2, 1:Mvol) = zero
    if (allocated(dItGpdxtp)) deallocate (dItGpdxtp)
    allocate (dItGpdxtp(0:1, -1:2, 1:Mvol), stat=astat)
    dItGpdxtp(0:1, -1:2, 1:Mvol) = zero

    if (allocated(glambda)) deallocate (glambda)
    allocate (glambda(1:Ntz + 1, 0:2, 0:1, 1:Mvol), stat=astat)
    glambda(1:Ntz + 1, 0:2, 0:1, 1:Mvol) = zero

    if (allocated(Bemn)) deallocate (Bemn)
    allocate (Bemn(1:mn, 1:Mvol, 0:1), stat=astat)
    Bemn(1:mn, 1:Mvol, 0:1) = zero
    if (allocated(Bomn)) deallocate (Bomn)
    allocate (Bomn(1:mn, 1:Mvol, 0:1), stat=astat)
    Bomn(1:mn, 1:Mvol, 0:1) = zero
    if (allocated(Iomn)) deallocate (Iomn)
    allocate (Iomn(1:mn, 1:Mvol), stat=astat)
    Iomn(1:mn, 1:Mvol) = zero
    if (allocated(Iemn)) deallocate (Iemn)
    allocate (Iemn(1:mn, 1:Mvol), stat=astat)
    Iemn(1:mn, 1:Mvol) = zero
    if (allocated(Somn)) deallocate (Somn)
    allocate (Somn(1:mn, 1:Mvol, 0:1), stat=astat)
    Somn(1:mn, 1:Mvol, 0:1) = zero
    if (allocated(Semn)) deallocate (Semn)
    allocate (Semn(1:mn, 1:Mvol, 0:1), stat=astat)
    Semn(1:mn, 1:Mvol, 0:1) = zero
    if (allocated(Pomn)) deallocate (Pomn)
    allocate (Pomn(1:mn, 1:Mvol, 0:2), stat=astat)
    Pomn(1:mn, 1:Mvol, 0:2) = zero
    if (allocated(Pemn)) deallocate (Pemn)
    allocate (Pemn(1:mn, 1:Mvol, 0:2), stat=astat)
    Pemn(1:mn, 1:Mvol, 0:2) = zero

    if (allocated(BBe)) deallocate (BBe)
    allocate (BBe(1:Mvol - 1), stat=astat)
    BBe(1:Mvol - 1) = zero
    if (allocated(IIo)) deallocate (IIo)
    allocate (IIo(1:Mvol - 1), stat=astat)
    IIo(1:Mvol - 1) = zero
    if (allocated(BBo)) deallocate (BBo)
    allocate (BBo(1:Mvol - 1), stat=astat)
    BBo(1:Mvol - 1) = zero
    if (allocated(IIe)) deallocate (IIe)
    allocate (IIe(1:Mvol - 1), stat=astat)
    IIe(1:Mvol - 1) = zero

    if (allocated(Btemn)) deallocate (Btemn)
    allocate (Btemn(1:mn, 0:1, 1:Mvol), stat=astat)
    Btemn(1:mn, 0:1, 1:Mvol) = zero
    if (allocated(Bzemn)) deallocate (Bzemn)
    allocate (Bzemn(1:mn, 0:1, 1:Mvol), stat=astat)
    Bzemn(1:mn, 0:1, 1:Mvol) = zero
    if (allocated(Btomn)) deallocate (Btomn)
    allocate (Btomn(1:mn, 0:1, 1:Mvol), stat=astat)
    Btomn(1:mn, 0:1, 1:Mvol) = zero
    if (allocated(Bzomn)) deallocate (Bzomn)
    allocate (Bzomn(1:mn, 0:1, 1:Mvol), stat=astat)
    Bzomn(1:mn, 0:1, 1:Mvol) = zero

    if (allocated(Bloweremn)) deallocate (Bloweremn)
    allocate (Bloweremn(1:mn, 3), stat=astat)
    Bloweremn(1:mn, 3) = zero
    if (allocated(Bloweromn)) deallocate (Bloweromn)
    allocate (Bloweromn(1:mn, 3), stat=astat)
    Bloweromn(1:mn, 3) = zero

    if (allocated(vvolume)) deallocate (vvolume)
    allocate (vvolume(1:Mvol), stat=astat)
    vvolume(1:Mvol) = zero
    if (allocated(lBBintegral)) deallocate (lBBintegral)
    allocate (lBBintegral(1:Mvol), stat=astat)
    lBBintegral(1:Mvol) = zero
    if (allocated(lABintegral)) deallocate (lABintegral)
    allocate (lABintegral(1:Mvol), stat=astat)
    lABintegral(1:Mvol) = zero

    if (YESstellsym) lmns = 1 + (mns - 1)
    if (NOTstellsym) lmns = 1 + (mns - 1) + (mns - 1)

    if (allocated(dlambdaout)) deallocate (dlambdaout)
    allocate (dlambdaout(1:lmns, 1:Mvol, 0:1), stat=astat)
    dlambdaout(1:lmns, 1:Mvol, 0:1) = zero

    if (Lconstraint == 3) then
        Localconstraint = .false.
    else
        Localconstraint = .true.
    end if

end subroutine preset

