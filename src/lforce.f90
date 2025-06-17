
subroutine lforce(lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag)

    use constants, only: zero, half, one, two

    use fileunits, only: ounit

    use inputlist, only: Wlforce, Igeometry, Nvol, Lrad, gamma, pscale, adiabatic, Lcheck

    use cputiming, only: Tlforce

    use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                         Lcoordinatesingularity, Mvol, &
                         iRbc, iZbs, iRbs, iZbc, &
                         YESstellsym, NOTstellsym, &
                         mn, im, in, regumm, &
                         ijreal, ijimag, jireal, jiimag, &
                         efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                         Nt, Nz, &
                         Ate, Aze, Ato, Azo, &
                         TT, RTT, &
                         sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, &
                         mmpp, &
                         Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                         Pomn, Pemn, &
                         vvolume, &
                         build_vector_potential

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, iocons, ideriv, Ntz, iflag
    real(8) :: dAt(1:Ntz, -1:2), dAz(1:Ntz, -1:2), XX(1:Ntz), YY(1:Ntz), dRR(1:Ntz, -1:1), dZZ(1:Ntz, -1:1), DDl, MMl

    real(8) :: IIl(1:Ntz), length(1:Ntz), dLL(1:Ntz)
    integer :: Lcurvature, ii, jj, kk, ll, ifail, ivol, lnn, mi, id !, oicons
    real(8) :: dBB(1:Ntz, -1:2), lss, mfactor

    real(8) :: dAs(1:Ntz) !, dRdt(-1:1,0:1), dZdt(-1:1,0:1)
    real(8) :: lgvuij(1:Ntz, 1:3, 1:3) ! local workspace; 13 Sep 13;

    dAt(1:Ntz, -1:2) = zero ! initialize intent out; 01 Jul 14;
    dAz(1:Ntz, -1:2) = zero ! initialize intent out; 01 Jul 14;

    lss = two*iocons - one ! recall that iocons is effective local radial coordinate; 24 Apr 13;

    Lcurvature = 1

    call coords(lvol, lss, Lcurvature, Ntz, mn) ! get coordinates and derivatives wrt Rj, Zj, at specific radial location;

    call build_vector_potential(lvol, iocons, 0, 1)

    call invfft(mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz, 0), dAz(1:Ntz, 0))

    id = ideriv
    if (id == 0) then

        dBB(1:Ntz, id) = half*(dAz(1:Ntz, id)*dAz(1:Ntz, id)*guvij(1:Ntz, 2, 2, id) &
                               - two*dAz(1:Ntz, id)*dAt(1:Ntz, id)*guvij(1:Ntz, 2, 3, id) &
                               + dAt(1:Ntz, id)*dAt(1:Ntz, id)*guvij(1:Ntz, 3, 3, id))/sg(1:Ntz, 0)**2

    else
        call build_vector_potential(lvol, iocons, id, 1)
        call invfft(mn, im, in, efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nt, Nz, dAt(1:Ntz, id), dAz(1:Ntz, id))

        dBB(1:Ntz, id) = half*(dAz(1:Ntz, id)*dAz(1:Ntz, 0)*guvij(1:Ntz, 2, 2, 0) &
                               - two*dAz(1:Ntz, id)*dAt(1:Ntz, 0)*guvij(1:Ntz, 2, 3, 0) &
                               + dAt(1:Ntz, id)*dAt(1:Ntz, 0)*guvij(1:Ntz, 3, 3, 0) &
                               + dAz(1:Ntz, 0)*dAz(1:Ntz, id)*guvij(1:Ntz, 2, 2, 0) &
                               - two*dAz(1:Ntz, 0)*dAt(1:Ntz, id)*guvij(1:Ntz, 2, 3, 0) &
                               + dAt(1:Ntz, 0)*dAt(1:Ntz, id)*guvij(1:Ntz, 3, 3, 0))/sg(1:Ntz, 0)**2
    end if ! end of if( ideriv.gt.0 ) ;

    if (ideriv == -1) then
        lss = two*iocons - one; Lcurvature = 4
        call coords(lvol, lss, Lcurvature, Ntz, mn)

        dBB(1:Ntz, id) = dBB(1:Ntz, id) + &
                         half*(dAz(1:Ntz, 0)*dAz(1:Ntz, 0)*guvij(1:Ntz, 2, 2, 1) &
                               - two*dAz(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz, 2, 3, 1) &
                               + dAt(1:Ntz, 0)*dAt(1:Ntz, 0)*guvij(1:Ntz, 3, 3, 1))/sg(1:Ntz, 0)**2 &
                         - dBB(1:Ntz, 0)*two*sg(1:Ntz, 1)/sg(1:Ntz, 0)
    end if

    ijreal(1:Ntz) = adiabatic(lvol)*pscale/vvolume(lvol)**gamma + dBB(1:Ntz, 0) ! p + B^2/2; 13 Sep 13;

    if (iflag == 1) return ! iflag = 1 indicates the derivatives of the force are to be calculated; derivatives of magnetic field calculated above;

    select case (Igeometry)

    case (1:2); dLL(1:Ntz) = zero ! placeholder; 08 Feb 16;
        ; ; IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;
    case (3)

        do ivol = 0, 1

            call invfft(mn, im(1:mn), in(1:mn), iRbc(1:mn, lvol - 1 + ivol), iRbs(1:mn, lvol - 1 + ivol), &
                        iZbc(1:mn, lvol - 1 + ivol), iZbs(1:mn, lvol - 1 + ivol), &
                        Nt, Nz, iRij(1:Ntz, lvol - 1 + ivol), iZij(1:Ntz, lvol - 1 + ivol))

            call invfft(mn, im(1:mn), in(1:mn), im(1:mn)*iRbs(1:mn, lvol - 1 + ivol), -im(1:mn)*iRbc(1:mn, lvol - 1 + ivol), &
                        im(1:mn)*iZbs(1:mn, lvol - 1 + ivol), -im(1:mn)*iZbc(1:mn, lvol - 1 + ivol), &
                        Nt, Nz, tRij(1:Ntz, lvol - 1 + ivol), tZij(1:Ntz, lvol - 1 + ivol))
        end do ! end of do ivol = 0, 1 ; 18 Jul 14;

        dRij(1:Ntz, lvol) = iRij(1:Ntz, lvol) - iRij(1:Ntz, lvol - 1)
        dZij(1:Ntz, lvol) = iZij(1:Ntz, lvol) - iZij(1:Ntz, lvol - 1)

        length(1:Ntz) = sqrt(dRij(1:Ntz, lvol)**2 + dZij(1:Ntz, lvol)**2)

        dLL(1:Ntz) = (dRij(1:Ntz, lvol)*tRij(1:Ntz, lvol - 1 + iocons) + dZij(1:Ntz, lvol)*tZij(1:Ntz, lvol - 1 + iocons))/length(1:Ntz)

        if (iocons == 1) then ! include spectral condensation constraints; local to interface, i.e. no tri-diagonal structure;
            ; ; efmn(1:mn) = (mmpp(1:mn))*iRbc(1:mn, lvol)
            ; ; sfmn(1:mn) = (mmpp(1:mn))*iZbs(1:mn, lvol)
            if (NOTstellsym) then; ofmn(1:mn) = (mmpp(1:mn))*iRbs(1:mn, lvol)
                ; ; cfmn(1:mn) = (mmpp(1:mn))*iZbc(1:mn, lvol)
            else; ofmn(1:mn) = zero
                ; ; cfmn(1:mn) = zero
            end if ! end of if( NOTstellsym ) ; 20 Feb 13;

            call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                        Nt, Nz, XX(1:Ntz), YY(1:Ntz))

            if (YESstellsym) then; DDl = sum((iRbc(1:mn, lvol)**2 + iZbs(1:mn, lvol)**2))
                ; ; MMl = sum(mmpp(1:mn)*(iRbc(1:mn, lvol)**2 + iZbs(1:mn, lvol)**2))/DDl
            else; DDl = sum((iRbc(1:mn, lvol)**2 + iZbs(1:mn, lvol)**2 + iRbs(1:mn, lvol)**2 + iZbc(1:mn, lvol)**2))
                ; ; MMl = sum(mmpp(1:mn)*(iRbc(1:mn, lvol)**2 + iZbs(1:mn, lvol)**2 + iRbs(1:mn, lvol)**2 + iZbc(1:mn, lvol)**2))/DDl
            end if

            IIl(1:Ntz) = tRij(1:Ntz, lvol)*(XX(1:Ntz) - MMl*iRij(1:Ntz, lvol)) &
                         + tZij(1:Ntz, lvol)*(YY(1:Ntz) - MMl*iZij(1:Ntz, lvol))

        else ! matches if( iocons.eq.1 ) ; 11 Aug 14;

            IIl(1:Ntz) = zero ! placeholder; 11 Aug 14;

        end if ! end of if( iocons.eq.1 ) ; 20 Feb 13;

    end select ! end of select case( Igeometry ) ; 08 Feb 16;

    ; ifail = 0
    ; call tfft(Nt, Nz, ijreal(1:Ntz), IIl(1:Ntz), & ! compute force-imbalance and spectral constraints;
 mn, im(1:mn), in(1:mn), Bemn(1:mn, lvol, iocons), Bomn(1:mn, lvol, iocons), Iemn(1:mn, lvol), Iomn(1:mn, lvol), ifail)

    if (Igeometry >= 3) then ! add minimal length constraint; 18 Jul 14;

        ifail = 0; ijimag(1:Ntz) = zero

        call tfft(Nt, Nz, dLL(1:Ntz), ijimag(1:Ntz), &
                  mn, im(1:mn), in(1:mn), Semn(1:mn, lvol, iocons), Somn(1:mn, lvol, iocons), Pemn(1:mn, lvol, iocons), Pomn(1:mn, lvol, iocons), ifail)

    end if ! end of if( Igeometry.eq.3 ) ; 01 Jul 14;

end subroutine lforce

