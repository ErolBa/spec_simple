
subroutine rzaxis(Mvol, mn, inRbc, inZbs, inRbs, inZbc, ivol, LcomputeDerivatives)

    use constants, only: zero, one, half, two

    use numerical, only: vsmall

    use fileunits, only: ounit

    use inputlist, only: Wrzaxis, Igeometry, Ntor, Lcheck, Wmacros, Lreflect, Ntoraxis, Lrzaxis

    use allglobal, only: ncpu, myid, cpus, im, in, &
                         ajk, Nt, Nz, Ntz, &
                         Rij, Zij, sg, cosi, sini, &
                         ijreal, ijimag, jireal, jiimag, jkreal, jkimag, kjreal, kjimag, &
                         efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, cosi, sini, &
                         YESstellsym, NOTstellsym, Lcoordinatesingularity, &
                         dRodR, dRodZ, dZodR, dZodZ, &
                         dRadR, dRadZ, dZadR, dZadZ, &
                         iRbc, iZbs, iRbs, iZbc, &
                         dBdX

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LComputeDerivatives

    integer, intent(in) :: Mvol, mn, ivol
    real(8) :: inRbc(1:mn, 0:Mvol), inZbs(1:mn, 0:Mvol), inRbs(1:mn, 0:Mvol), inZbc(1:mn, 0:Mvol)
    real(8) :: jRbc(1:mn, 0:Mvol), jZbs(1:mn, 0:Mvol), jRbs(1:mn, 0:Mvol), jZbc(1:mn, 0:Mvol)
    real(8) :: tmpRbc(1:mn, 0:Mvol), tmpZbs(1:mn, 0:Mvol), tmpRbs(1:mn, 0:Mvol), tmpZbc(1:mn, 0:Mvol)

    real(8) :: jacbase(1:Ntz), jacbasec(1:mn), jacbases(1:mn)
    real(8) :: junkc(1:mn), junks(1:mn)

    integer :: jvol, ii, ifail, jj, id, issym, irz, imn
    integer :: idJc, idJs, idRc, idRs, idZc, idZs

    integer :: Lcurvature

    integer :: Njac, idgetrf, idgetrs
    real(8), allocatable :: jacrhs(:), djacrhs(:), jacmat(:, :), djacmat(:, :), solution(:), LU(:, :)
    integer, allocatable :: ipiv(:)

    write (*, *) "Calling rzaxis"

    jvol = 0

    Ntoraxis = min(Ntor, Ntoraxis)

    select case (Igeometry)

    case (1:2)

        inRbc(1:mn, jvol) = zero
        inRbs(1:mn, jvol) = zero

        if (Igeometry == 1 .and. Lreflect == 1) then
            inRbc(2:mn, 0) = -inRbc(2:mn, Mvol)
            if (NOTstellsym) then
                inRbs(2:mn, 0) = -inRbs(2:mn, Mvol)
            end if
        end if

    case (3)

        if (Lrzaxis == 1) then

            call invfft(mn, im(1:mn), in(1:mn), im(1:mn)*inRbs(1:mn, ivol), -im(1:mn)*inRbc(1:mn, ivol), &
                        im(1:mn)*inZbs(1:mn, ivol), -im(1:mn)*inZbc(1:mn, ivol), &
                        Nt, Nz, jkreal(1:Ntz), jkimag(1:Ntz))

            ijreal(1:Ntz) = sqrt(jkreal(1:Ntz)**2 + jkimag(1:Ntz)**2)
            ijimag(1:Ntz) = zero

            jireal(1:Ntz) = ijreal(1:Ntz)

            ifail = 0
            call tfft(Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                      mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail)

            efmn(1:mn) = efmn(1:mn)*ajk(1:mn)
            ofmn(1:mn) = ofmn(1:mn)*ajk(1:mn)
            cfmn(1:mn) = zero
            sfmn(1:mn) = zero

            call invfft(mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                        Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz))

            jiimag(1:Ntz) = ijreal(1:Ntz)

            call invfft(mn, im(1:mn), in(1:mn), inRbc(1:mn, ivol), inRbs(1:mn, ivol), &
                        inZbc(1:mn, ivol), inZbs(1:mn, ivol), &
                        Nt, Nz, kjreal(1:Ntz), kjimag(1:Ntz))

            ijreal(1:Ntz) = kjreal(1:Ntz)*jireal(1:Ntz)
            ijimag(1:Ntz) = kjimag(1:Ntz)*jireal(1:Ntz)

            ifail = 0
            call tfft(Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                      mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail)

            evmn(1:mn) = evmn(1:mn)*ajk(1:mn)
            odmn(1:mn) = odmn(1:mn)*ajk(1:mn)
            comn(1:mn) = comn(1:mn)*ajk(1:mn)
            simn(1:mn) = simn(1:mn)*ajk(1:mn)

            call invfft(mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), &
                        Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz))

            ijreal(1:Ntz) = ijreal(1:Ntz)/jiimag(1:Ntz)
            ijimag(1:Ntz) = ijimag(1:Ntz)/jiimag(1:Ntz)

            kjreal(1:Ntz) = kjreal(1:Ntz) - ijreal(1:Ntz)
            kjimag(1:Ntz) = kjimag(1:Ntz) - ijimag(1:Ntz)

            ifail = 0
            call tfft(Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                      mn, im(1:mn), in(1:mn), inRbc(1:mn, jvol), inRbs(1:mn, jvol), inZbc(1:mn, jvol), inZbs(1:mn, jvol), ifail)

        else if (Lrzaxis == 2) then

            tmpRbc = iRbc
            tmpZbs = iZbs
            tmpRbs = iRbs
            tmpZbc = iZbc

            jRbc = inRbc
            jZbs = inZbs
            jRbs = inRbs
            jZbc = inZbc

            if (YESstellsym) then
                Njac = 2*Ntoraxis + 1
            else
                Njac = 2*(2*Ntoraxis + 1)
            end if

            if (allocated(jacrhs)) deallocate (jacrhs)
            allocate (jacrhs(1:Njac), stat=astat)
            jacrhs(1:Njac) = zero
            if (allocated(jacmat)) deallocate (jacmat)
            allocate (jacmat(1:Njac, 1:Njac), stat=astat)
            jacmat(1:Njac, 1:Njac) = zero
            if (allocated(LU)) deallocate (LU)
            allocate (LU(1:Njac, 1:Njac), stat=astat)
            LU(1:Njac, 1:Njac) = zero
            if (allocated(solution)) deallocate (solution)
            allocate (solution(1:Njac), stat=astat)
            solution(1:Njac) = zero
            if (allocated(ipiv)) deallocate (ipiv)
            allocate (ipiv(1:Njac), stat=astat)
            ipiv(1:Njac) = 0

            iRbc(1:mn, 1) = jRbc(1:mn, ivol)
            iZbs(1:mn, 1) = jZbs(1:mn, ivol)
            iRbs(1:mn, 1) = jRbs(1:mn, ivol)
            iZbc(1:mn, 1) = jZbc(1:mn, ivol)

            iRbc(1:mn, 0) = zero
            iZbs(1:mn, 0) = zero
            iRbs(1:mn, 0) = zero
            iZbc(1:mn, 0) = zero

            iRbc(1:Ntor + 1, 0) = jRbc(1:Ntor + 1, ivol)
            iZbs(1:Ntor + 1, 0) = jZbs(1:Ntor + 1, ivol)
            iRbs(1:Ntor + 1, 0) = jRbs(1:Ntor + 1, ivol)
            iZbc(1:Ntor + 1, 0) = jZbc(1:Ntor + 1, ivol)

            Lcoordinatesingularity = .true.
            Lcurvature = 1
            dBdX%innout = 1

            idJc = Ntoraxis + 1
            idJs = Ntoraxis + 1 + 2*Ntoraxis + 1
            idRc = 1
            idZs = Ntoraxis + 1
            idRs = 2*Ntoraxis + 1
            idZc = 3*Ntoraxis + 2

            call coords(1, one, Lcurvature, Ntz, mn)

            jacbase = sg(1:Ntz, 0)/Rij(1:Ntz, 0, 0)

            call tfft(Nt, Nz, jacbase, Rij, &
                      mn, im(1:mn), in(1:mn), jacbasec(1:mn), jacbases(1:mn), junkc(1:mn), junks(1:mn), ifail)

            if (YESstellsym) then
                jacrhs = -jacbasec(2*(Ntor + 1) - Ntoraxis:2*(Ntor + 1) + Ntoraxis)
            else
                jacrhs(1:2*Ntoraxis + 1) = -jacbasec(2*(Ntor + 1) - Ntoraxis:2*(Ntor + 1) + Ntoraxis)
                jacrhs(2*Ntoraxis + 2:Njac) = -jacbases(2*(Ntor + 1) - Ntoraxis:2*(Ntor + 1) + Ntoraxis)
            end if

            if (YESstellsym) then

                do ii = -Ntoraxis, Ntoraxis
                    do jj = 1, Ntoraxis

                        if (ii - jj >= -Ntor) then
                            id = 2*(Ntor + 1) + ii - jj
                            jacmat(ii + Ntoraxis + 1, jj + 1) = jacmat(ii + Ntoraxis + 1, jj + 1) - jZbs(id, ivol)
                            jacmat(ii + Ntoraxis + 1, Ntoraxis + 1 + jj) = jacmat(ii + Ntoraxis + 1, Ntoraxis + 1 + jj) + jRbc(id, ivol)
                        end if

                        if (ii + jj <= Ntor) then
                            id = 2*(Ntor + 1) + ii + jj
                            jacmat(ii + Ntoraxis + 1, jj + 1) = jacmat(ii + Ntoraxis + 1, jj + 1) - jZbs(id, ivol)
                            jacmat(ii + Ntoraxis + 1, Ntoraxis + 1 + jj) = jacmat(ii + Ntoraxis + 1, Ntoraxis + 1 + jj) - jRbc(id, ivol)
                        end if

                    end do

                    id = 2*(Ntor + 1) + ii
                    jacmat(ii + Ntoraxis + 1, 1) = -two*jZbs(id, ivol)

                end do

            else

                do ii = -Ntoraxis, Ntoraxis
                    do jj = 1, Ntoraxis

                        if (ii - jj >= -Ntor) then
                            id = 2*(Ntor + 1) + ii - jj
                            jacmat(ii + idJc, jj + idRc) = jacmat(ii + idJc, jj + idRc) - jZbs(id, ivol)
                            jacmat(ii + idJc, jj + idZs) = jacmat(ii + idJc, jj + idZs) + jRbc(id, ivol)
                            jacmat(ii + idJc, jj + idRs) = jacmat(ii + idJc, jj + idRs) - jZbc(id, ivol)
                            jacmat(ii + idJc, jj + idZc) = jacmat(ii + idJc, jj + idZc) + jRbs(id, ivol)

                            jacmat(ii + idJs, jj + idRc) = jacmat(ii + idJs, jj + idRc) + jZbc(id, ivol)
                            jacmat(ii + idJs, jj + idZs) = jacmat(ii + idJs, jj + idZs) + jRbs(id, ivol)
                            jacmat(ii + idJs, jj + idRs) = jacmat(ii + idJs, jj + idRs) - jZbs(id, ivol)
                            jacmat(ii + idJs, jj + idZc) = jacmat(ii + idJs, jj + idZc) - jRbc(id, ivol)

                        end if

                        if (ii + jj <= Ntor) then
                            id = 2*(Ntor + 1) + ii + jj

                            jacmat(ii + idJc, jj + idRc) = jacmat(ii + idJc, jj + idRc) - jZbs(id, ivol)
                            jacmat(ii + idJc, jj + idZs) = jacmat(ii + idJc, jj + idZs) - jRbc(id, ivol)
                            jacmat(ii + idJc, jj + idRs) = jacmat(ii + idJc, jj + idRs) + jZbc(id, ivol)
                            jacmat(ii + idJc, jj + idZc) = jacmat(ii + idJc, jj + idZc) + jRbs(id, ivol)

                            jacmat(ii + idJs, jj + idRc) = jacmat(ii + idJs, jj + idRc) + jZbc(id, ivol)
                            jacmat(ii + idJs, jj + idZs) = jacmat(ii + idJs, jj + idZs) - jRbs(id, ivol)
                            jacmat(ii + idJs, jj + idRs) = jacmat(ii + idJs, jj + idRs) + jZbs(id, ivol)
                            jacmat(ii + idJs, jj + idZc) = jacmat(ii + idJs, jj + idZc) - jRbc(id, ivol)
                        end if

                    end do

                    id = 2*(Ntor + 1) + ii
                    jacmat(ii + idJc, idRc) = -two*jZbs(id, ivol)
                    jacmat(ii + idJc, idZc) = +two*jRbs(id, ivol)
                    jacmat(ii + idJs, idRc) = +two*jZbc(id, ivol)
                    jacmat(ii + idJs, idZc) = -two*jRbc(id, ivol)

                end do

            end if

            jacmat = jacmat*half

            LU = jacmat
            call DGETRF(Njac, Njac, LU, Njac, ipiv, idgetrf)
            solution = jacrhs
            call DGETRS('N', Njac, 1, LU, Njac, ipiv, solution, Njac, idgetrs)

            if (idgetrf < 0 .or. idgetrs < 0) then
                ; write (ounit, 1010) cput - cpus, myid, ivol, idgetrf, idgetrs, "input error ;     "
            elseif (idgetrf > 0) then
                ; write (ounit, 1010) cput - cpus, myid, ivol, idgetrf, idgetrs, "singular ;        "
            end if

1010        format("rzaxis : ", f10.2, " : myid=", i3, " ; ivol=", i3, " idgetrf idgetrs=", i3, ' ', i3, " ; "a34)

            iRbc = tmpRbc
            iZbs = tmpZbs
            iRbs = tmpRbs
            iZbc = tmpZbc

            inRbc(:, jvol) = zero
            inZbs(:, jvol) = zero
            inRbs(:, jvol) = zero
            inZbc(:, jvol) = zero

            inRbc(1:Ntoraxis + 1, jvol) = inRbc(1:Ntoraxis + 1, ivol) - solution(idRc:idRc + Ntoraxis)
            inZbs(2:Ntoraxis + 1, jvol) = inZbs(2:Ntoraxis + 1, ivol) - solution(idZs + 1:idZs + Ntoraxis)
            if (YESstellsym) then
                inRbs(1:Ntoraxis + 1, jvol) = zero
                inZbc(2:Ntoraxis + 1, jvol) = zero
            else
                inRbs(2:Ntoraxis + 1, jvol) = inRbs(2:Ntoraxis + 1, ivol) - solution(idRs + 1:idRs + Ntoraxis)
                inZbc(1:Ntoraxis + 1, jvol) = inZbc(1:Ntoraxis + 1, ivol) - solution(idZc:idZc + Ntoraxis)
            end if

            deallocate (jacrhs, stat=astat)
            deallocate (jacmat, stat=astat)
            deallocate (LU, stat=astat)
            deallocate (solution, stat=astat)
            deallocate (ipiv, stat=astat)

        end if

    end select

end subroutine rzaxis

