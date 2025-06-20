
subroutine ma00aa(lquad, mn, lvol, lrad)

    use constants, only: zero, half, one, two, pi, pi2

    use fileunits, only: ounit

    use inputlist, only: mpol, Wma00aa, Wmacros

    use allglobal, only: myid, ncpu, cpus, &
                         Mvol, im, in, mne, &
                         YESstellsym, NOTstellsym, &
                         gaussianweight, gaussianabscissae, &
                         DToocc, DToocs, DToosc, DTooss, &
                         TTsscc, TTsscs, TTsssc, TTssss, &
                         TDstcc, TDstcs, TDstsc, TDstss, &
                         TDszcc, TDszcs, TDszsc, TDszss, &
                         DDttcc, DDttcs, DDttsc, DDttss, &
                         DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                         DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                         ki, kija, kijs, &
                         goomne, goomno, &
                         gssmne, gssmno, &
                         gstmne, gstmno, &
                         gszmne, gszmno, &
                         gttmne, gttmno, &
                         gtzmne, gtzmno, &
                         gzzmne, gzzmno, &
                         Lcoordinatesingularity, regumm, &
                         pi2pi2nfp, pi2pi2nfphalf, Lsavedguvij, &
                         dBdX

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lquad, mn, lvol, lrad

    integer :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele

    integer :: kk, kd, kka, kks, kda, kds, Lcurvature, ideriv

    real(8) :: lss, jthweight, fee, feo, foe, foo, Tl, Dl, Tp, Dp, TlTp, TlDp, DlTp, DlDp, ikda, ikds, imn2, ilrad, lssm

    real(8) :: foocc, foocs, foosc, fooss
    real(8) :: fsscc, fsscs, fsssc, fssss
    real(8) :: fstcc, fstcs, fstsc, fstss
    real(8) :: fszcc, fszcs, fszsc, fszss
    real(8) :: fttcc, fttcs, fttsc, fttss
    real(8) :: ftzcc, ftzcs, ftzsc, ftzss
    real(8) :: fzzcc, fzzcs, fzzsc, fzzss

    real(8) :: sbar
    real(8), allocatable :: basis(:, :, :, :)

    mn2_max = mn*mn
    lp2_max = (lrad + 1)*(lrad + 1)
    imn2 = one/real(mn)
    ilrad = one/real(lrad + 1)

    DToocc = zero
    TTssss = zero
    TDstsc = zero
    TDszsc = zero
    DDttcc = zero
    DDtzcc = zero
    DDzzcc = zero

    if (NOTstellsym) then
        DToocs = zero
        DToosc = zero
        DTooss = zero

        TTsscc = zero
        TTsscs = zero
        TTsssc = zero

        TDstcc = zero
        TDstcs = zero
        TDstss = zero

        TDszcc = zero
        TDszcs = zero
        TDszss = zero

        DDttcs = zero
        DDttsc = zero
        DDttss = zero

        DDtzcs = zero
        DDtzsc = zero
        DDtzss = zero

        DDzzcs = zero
        DDzzsc = zero
        DDzzss = zero
    end if

    if (allocated(basis)) deallocate (basis)
    allocate (basis(0:lrad, 0:mpol, 0:1, lquad), stat=astat)
    basis(0:lrad, 0:mpol, 0:1, lquad) = zero

    if (dBdX%L) then; Lcurvature = 3; ideriv = 1
    else; Lcurvature = 1; ideriv = 0
    end if

    if (.not. Lsavedguvij) then
        call compute_guvijsave(lquad, lvol, ideriv, Lcurvature)
    end if
    call metrix(lquad, lvol)

    do jquad = 1, lquad
        lss = gaussianabscissae(jquad, lvol); jthweight = gaussianweight(jquad, lvol)
        sbar = (lss + one)*half
        if (Lcoordinatesingularity) then
            call get_zernike(sbar, lrad, mpol, basis(:, :, 0:1, jquad))
        else
            call get_cheby(lss, lrad, basis(:, 0, 0:1, jquad))
        end if
    end do

    do mn2 = 1, mn2_max
        ii = mod(mn2 - 1, mn) + 1
        jj = (mn2 - ii)/mn + 1

        do jquad = 1, lquad

            lss = gaussianabscissae(jquad, lvol); jthweight = gaussianweight(jquad, lvol)
            sbar = (lss + one)*half

            kks = kijs(ii, jj, 0)
            kka = kija(ii, jj, 0)
            ikds = jthweight/kijs(ii, jj, 1)
            ikda = jthweight/kija(ii, jj, 1)

            foocc = (+goomne(kks, jquad)*abs(ikds) + goomne(kka, jquad)*abs(ikda))
            fssss = (+gssmne(kks, jquad)*abs(ikds) - gssmne(kka, jquad)*abs(ikda))
            fstsc = (+gstmno(kks, jquad)*ikds + gstmno(kka, jquad)*ikda)
            fszsc = (+gszmno(kks, jquad)*ikds + gszmno(kka, jquad)*ikda)
            fttcc = (+gttmne(kks, jquad)*abs(ikds) + gttmne(kka, jquad)*abs(ikda))
            ftzcc = (+gtzmne(kks, jquad)*abs(ikds) + gtzmne(kka, jquad)*abs(ikda))
            fzzcc = (+gzzmne(kks, jquad)*abs(ikds) + gzzmne(kka, jquad)*abs(ikda))

            if (NOTstellsym) then
                foocs = (-goomno(kks, jquad)*ikds + goomno(kka, jquad)*ikda)
                foosc = (+goomno(kks, jquad)*ikds + goomno(kka, jquad)*ikda)
                fooss = (+goomne(kks, jquad)*abs(ikds) - goomne(kka, jquad)*abs(ikda))

                fsscc = (+gssmne(kks, jquad)*abs(ikds) + gssmne(kka, jquad)*abs(ikda))
                fsscs = (-gssmno(kks, jquad)*ikds + gssmno(kka, jquad)*ikda)
                fsssc = (+gssmno(kks, jquad)*ikds + gssmno(kka, jquad)*ikda)

                fstcc = (+gstmne(kks, jquad)*abs(ikds) + gstmne(kka, jquad)*abs(ikda))
                fstcs = (-gstmno(kks, jquad)*ikds + gstmno(kka, jquad)*ikda)
                fstss = (+gstmne(kks, jquad)*abs(ikds) - gstmne(kka, jquad)*abs(ikda))

                fszcc = (+gszmne(kks, jquad)*abs(ikds) + gszmne(kka, jquad)*abs(ikda))
                fszcs = (-gszmno(kks, jquad)*ikds + gszmno(kka, jquad)*ikda)
                fszss = (+gszmne(kks, jquad)*abs(ikds) - gszmne(kka, jquad)*abs(ikda))

                fttcs = (-gttmno(kks, jquad)*ikds + gttmno(kka, jquad)*ikda)
                fttsc = (+gttmno(kks, jquad)*ikds + gttmno(kka, jquad)*ikda)
                fttss = (+gttmne(kks, jquad)*abs(ikds) - gttmne(kka, jquad)*abs(ikda))

                ftzcs = (-gtzmno(kks, jquad)*ikds + gtzmno(kka, jquad)*ikda)
                ftzsc = (+gtzmno(kks, jquad)*ikds + gtzmno(kka, jquad)*ikda)
                ftzss = (+gtzmne(kks, jquad)*abs(ikds) - gtzmne(kka, jquad)*abs(ikda))

                fzzcs = (-gzzmno(kks, jquad)*ikds + gzzmno(kka, jquad)*ikda)
                fzzsc = (+gzzmno(kks, jquad)*ikds + gzzmno(kka, jquad)*ikda)
                fzzss = (+gzzmne(kks, jquad)*abs(ikds) - gzzmne(kka, jquad)*abs(ikda))
            end if

            do lp2 = 1, lp2_max
                ll = mod(lp2 - 1, lrad + 1)
                pp = (lp2 - ll - 1)/(lrad + 1)

                if (Lcoordinatesingularity) then

                    ll1 = (ll - mod(ll, 2))/2
                    pp1 = (pp - mod(pp, 2))/2

                    if (ll < im(ii)) cycle
                    if (pp < im(jj)) cycle
                    if (mod(ll + im(ii), 2) /= 0) cycle
                    if (mod(pp + im(jj), 2) /= 0) cycle

                    Tl = basis(ll, im(ii), 0, jquad)
                    Dl = basis(ll, im(ii), 1, jquad)*half

                    Tp = basis(pp, im(jj), 0, jquad)
                    Dp = basis(pp, im(jj), 1, jquad)*half

                else

                    ll1 = ll
                    pp1 = pp

                    Tl = basis(ll, 0, 0, jquad)
                    Dl = basis(ll, 0, 1, jquad)

                    Tp = basis(pp, 0, 0, jquad)
                    Dp = basis(pp, 0, 1, jquad)

                end if

                TlTp = Tl*Tp
                TlDp = Tl*Dp
                DlTp = Dl*Tp
                DlDp = Dl*Dp

                DToocc(ll1, pp1, ii, jj) = DToocc(ll1, pp1, ii, jj) + DlTp*foocc
                TTssss(ll1, pp1, ii, jj) = TTssss(ll1, pp1, ii, jj) + TlTp*fssss
                TDstsc(ll1, pp1, ii, jj) = TDstsc(ll1, pp1, ii, jj) + TlDp*fstsc
                TDszsc(ll1, pp1, ii, jj) = TDszsc(ll1, pp1, ii, jj) + TlDp*fszsc
                DDttcc(ll1, pp1, ii, jj) = DDttcc(ll1, pp1, ii, jj) + DlDp*fttcc
                DDtzcc(ll1, pp1, ii, jj) = DDtzcc(ll1, pp1, ii, jj) + DlDp*ftzcc
                DDzzcc(ll1, pp1, ii, jj) = DDzzcc(ll1, pp1, ii, jj) + DlDp*fzzcc

                if (NOTstellsym) then

                    DToocs(ll1, pp1, ii, jj) = DToocs(ll1, pp1, ii, jj) + DlTp*foocs
                    DToosc(ll1, pp1, ii, jj) = DToosc(ll1, pp1, ii, jj) + DlTp*foosc
                    DTooss(ll1, pp1, ii, jj) = DTooss(ll1, pp1, ii, jj) + DlTp*fooss

                    TTsscc(ll1, pp1, ii, jj) = TTsscc(ll1, pp1, ii, jj) + TlTp*fsscc
                    TTsscs(ll1, pp1, ii, jj) = TTsscs(ll1, pp1, ii, jj) + TlTp*fsscs
                    TTsssc(ll1, pp1, ii, jj) = TTsssc(ll1, pp1, ii, jj) + TlTp*fsssc

                    TDstcc(ll1, pp1, ii, jj) = TDstcc(ll1, pp1, ii, jj) + TlDp*fstcc
                    TDstcs(ll1, pp1, ii, jj) = TDstcs(ll1, pp1, ii, jj) + TlDp*fstcs
                    TDstss(ll1, pp1, ii, jj) = TDstss(ll1, pp1, ii, jj) + TlDp*fstss

                    TDszcc(ll1, pp1, ii, jj) = TDszcc(ll1, pp1, ii, jj) + TlDp*fszcc
                    TDszcs(ll1, pp1, ii, jj) = TDszcs(ll1, pp1, ii, jj) + TlDp*fszcs
                    TDszss(ll1, pp1, ii, jj) = TDszss(ll1, pp1, ii, jj) + TlDp*fszss

                    DDttcs(ll1, pp1, ii, jj) = DDttcs(ll1, pp1, ii, jj) + DlDp*fttcs
                    DDttsc(ll1, pp1, ii, jj) = DDttsc(ll1, pp1, ii, jj) + DlDp*fttsc
                    DDttss(ll1, pp1, ii, jj) = DDttss(ll1, pp1, ii, jj) + DlDp*fttss

                    DDtzcs(ll1, pp1, ii, jj) = DDtzcs(ll1, pp1, ii, jj) + DlDp*ftzcs
                    DDtzsc(ll1, pp1, ii, jj) = DDtzsc(ll1, pp1, ii, jj) + DlDp*ftzsc
                    DDtzss(ll1, pp1, ii, jj) = DDtzss(ll1, pp1, ii, jj) + DlDp*ftzss

                    DDzzcs(ll1, pp1, ii, jj) = DDzzcs(ll1, pp1, ii, jj) + DlDp*fzzcs
                    DDzzsc(ll1, pp1, ii, jj) = DDzzsc(ll1, pp1, ii, jj) + DlDp*fzzsc
                    DDzzss(ll1, pp1, ii, jj) = DDzzss(ll1, pp1, ii, jj) + DlDp*fzzss
                end if

            end do

        end do

    end do

    deallocate (basis, stat=astat)

    DToocc = DToocc*pi2pi2nfphalf
    TTssss = TTssss*pi2pi2nfphalf
    TDstsc = TDstsc*pi2pi2nfphalf
    TDszsc = TDszsc*pi2pi2nfphalf
    DDttcc = DDttcc*pi2pi2nfphalf
    DDtzcc = DDtzcc*pi2pi2nfphalf
    DDzzcc = DDzzcc*pi2pi2nfphalf

    if (NOTstellsym) then
        DToocs = DToocs*pi2pi2nfphalf
        DToosc = DToosc*pi2pi2nfphalf
        DTooss = DTooss*pi2pi2nfphalf

        TTsscc = TTsscc*pi2pi2nfphalf
        TTsscs = TTsscs*pi2pi2nfphalf
        TTsssc = TTsssc*pi2pi2nfphalf

        TDstcc = TDstcc*pi2pi2nfphalf
        TDstcs = TDstcs*pi2pi2nfphalf
        TDstss = TDstss*pi2pi2nfphalf

        TDszcc = TDszcc*pi2pi2nfphalf
        TDszcs = TDszcs*pi2pi2nfphalf
        TDszss = TDszss*pi2pi2nfphalf

        DDttcs = DDttcs*pi2pi2nfphalf
        DDttsc = DDttsc*pi2pi2nfphalf
        DDttss = DDttss*pi2pi2nfphalf

        DDtzcs = DDtzcs*pi2pi2nfphalf
        DDtzsc = DDtzsc*pi2pi2nfphalf
        DDtzss = DDtzss*pi2pi2nfphalf

        DDzzcs = DDzzcs*pi2pi2nfphalf
        DDzzsc = DDzzsc*pi2pi2nfphalf
        DDzzss = DDzzss*pi2pi2nfphalf

    end if

end subroutine ma00aa

