
subroutine dfp200(LcomputeDerivatives, vvol)

    use constants, only: zero, half, one, two

    use numerical, only: small

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wdfp200, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                         gamma, adiabatic, pscale, mu, &
                         epsilon, &
                         Lfindzero, &
                         Lconstraint, Lcheck, LHmatrix, &
                         Lextrap

    use cputiming, only: Tdfp200

    use allglobal, only: ncpu, myid, cpus, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                         Mvol, &
                         Iquad, &
                         iRbc, iZbs, iRbs, iZbc, &
                         NAdof, &
                         YESstellsym, NOTstellsym, &
                         mn, im, in, mns, Ntz, &
                         Ate, Aze, Ato, Azo, &
                         ijreal, &
                         efmn, ofmn, cfmn, sfmn, &
                         evmn, odmn, comn, simn, &
                         Nt, Nz, &
                         cosi, sini, &
                         dBdX, &
                         dtflux, dpflux, sweight, &
                         mmpp, &
                         dMA, dMB, dMD, dMG, &
                         Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                         LGdof, &
                         vvolume, dvolume, &
                         Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, &
                         diotadxup, dItGpdxtp, &
                         dFFdRZ, dBBdmp, dmupfdx, denergydrr, denergydzr, &
                         BBweight, &
                         psifactor, &
                         lmns, &
                         mn, mne, &
                         dRodR, dRodZ, dZodR, dZodZ, &
                         LocalConstraint, solution, &
                         Btemn

    use typedefns

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LComputeDerivatives
    logical :: LInnerVolume

    integer :: NN, IA, ifail, if01adf, vflag, MM, idgetrf, idgetri, Lwork, lvol, pvol
    integer :: vvol, innout, ii, jj, irz, issym, iocons, idoc, idof, imn, ll
    integer :: Lcurvature, ideriv, id
    integer :: iflag, cpu_id, cpu_id1, even_or_odd, vol_parity
    integer :: tag, tag2, req1, req2, req3, req4

    real(8) :: lastcpu, lss, lfactor, DDl, MMl
    real(8) :: det
    real(8), allocatable :: XX(:), YY(:), dBB(:, :), dII(:), dLL(:), dPP(:), length(:), dRR(:, :), dZZ(:, :), constraint(:)
    real(8), allocatable :: ddFcol1(:), ddFcol2(:), ddFcol3(:), ddFcol4(:)

    character :: packorunpack

    type(MatrixLU) :: oBI(1:Mvol)

    if (allocated(dBB)) deallocate (dBB)
    allocate (dBB(1:Ntz, -1:2), stat=astat)
    dBB(1:Ntz, -1:2) = zero
    if (allocated(XX)) deallocate (XX)
    allocate (XX(1:Ntz), stat=astat)
    XX(1:Ntz) = zero
    if (allocated(YY)) deallocate (YY)
    allocate (YY(1:Ntz), stat=astat)
    YY(1:Ntz) = zero
    if (allocated(length)) deallocate (length)
    allocate (length(1:Ntz), stat=astat)
    length(1:Ntz) = zero
    if (LocalConstraint) then

        do vvol = 1, Mvol

            if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
            else; Lcoordinatesingularity = .true.
            end if

            dBdX%vol = vvol
            ll = Lrad(vvol)
            NN = NAdof(vvol)

            vflag = 1
            call volume(vvol, vflag)

            do iocons = 0, 1

                if (vvol == 1 .and. iocons == 0) cycle
                if (vvol == Mvol .and. iocons == 1) cycle

                ideriv = 0; id = ideriv
                iflag = 0
                call lforce(vvol, iocons, ideriv, Ntz, dBB(1:Ntz, id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag)

            end do

        end do

    else

        do vvol = 1, Mvol
            if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
            else; Lcoordinatesingularity = .true.
            end if
            ll = Lrad(vvol)
            NN = NAdof(vvol)

            vflag = 1
            call volume(vvol, vflag)

            do iocons = 0, 1

                if (vvol == 1 .and. iocons == 0) cycle
                if (vvol == Mvol .and. iocons == 1) cycle

                ideriv = 0; id = ideriv
                iflag = 0
                call lforce(vvol, iocons, ideriv, Ntz, dBB(1:Ntz, id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag)

            end do
        end do

    end if

    deallocate (dBB, stat=astat)
    deallocate (XX, stat=astat)
    deallocate (YY, stat=astat)
    deallocate (length, stat=astat)

2000 continue

end subroutine dfp200
