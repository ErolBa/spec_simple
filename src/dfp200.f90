
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

    use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                         Mvol, &
                         Iquad, & ! convenience; provided to ma00aa as argument to avoid allocations;
                         iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                         NAdof, &
                         YESstellsym, NOTstellsym, &
                         mn, im, in, mns, Ntz, &
                         Ate, Aze, Ato, Azo, & ! only required for debugging;
                         ijreal, &
                         efmn, ofmn, cfmn, sfmn, &
                         evmn, odmn, comn, simn, &
                         Nt, Nz, &
                         cosi, sini, & ! FFT workspace;
                         dBdX, &
                         dtflux, dpflux, sweight, &
                         mmpp, &
                         dMA, dMB, dMD, dMG, &
                         Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                         LGdof, &
                         vvolume, dvolume, &
                         Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                         diotadxup, dItGpdxtp, &
                         dFFdRZ, dBBdmp, dmupfdx, denergydrr, denergydzr, &
                         BBweight, & ! exponential weight on force-imbalance harmonics;
                         psifactor, &
                         lmns, &
                         mn, mne, &
                         dRodR, dRodZ, dZodR, dZodZ, &
                         LocalConstraint, solution, &
                         IsMyVolume, IsMyVolumeValue, Btemn, WhichCpuID

    use typedefns

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
    logical :: LInnerVolume

    integer :: NN, IA, ifail, if01adf, vflag, MM, idgetrf, idgetri, Lwork, lvol, pvol
    integer :: vvol, innout, ii, jj, irz, issym, iocons, idoc, idof, imn, ll
    integer :: Lcurvature, ideriv, id
    integer :: iflag, cpu_id, cpu_id1, even_or_odd, vol_parity
    integer :: stat(MPI_STATUS_SIZE), tag, tag2, req1, req2, req3, req4

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

            call IsMyVolume(vvol)

            if (IsMyVolumeValue == 0) then
                cycle
            else if (IsMyVolumeValue == -1) then
                if (.true.) then
                    write (6, '("dfp200 :      fatal : myid=",i3," ; .true. ; Unassociated volume;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "dfp200 : .true. : Unassociated volume ;"
                end if
            end if

            if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
            else; Lcoordinatesingularity = .true.
            end if ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

            dBdX%vol = vvol ! Label
            ll = Lrad(vvol) ! Shorthand
            NN = NAdof(vvol) ! shorthand;

            vflag = 1
            call volume(vvol, vflag) ! compute volume;

            do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;

                if (vvol == 1 .and. iocons == 0) cycle ! fixed inner boundary (or coordinate axis);
                if (vvol == Mvol .and. iocons == 1) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;

                ideriv = 0; id = ideriv
                iflag = 0 ! XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
                call lforce(vvol, iocons, ideriv, Ntz, dBB(1:Ntz, id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag)

            end do ! end of do iocons = 0, 1;

        end do ! matches do vvol = 1, Mvol

    else ! CASE SEMI GLOBAL CONSTRAINT

        do vvol = 1, Mvol
            call IsMyVolume(vvol)

            if (IsMyVolumeValue == 0) then
                cycle
            else if (IsMyVolumeValue == -1) then
                if (.true.) then
                    write (6, '("dfp200 :      fatal : myid=",i3," ; .true. ; Unassociated volume;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "dfp200 : .true. : Unassociated volume ;"
                end if
            end if

            if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
            else; Lcoordinatesingularity = .true.
            end if ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
            ll = Lrad(vvol) ! Shorthand
            NN = NAdof(vvol) ! shorthand;

            vflag = 1
            call volume(vvol, vflag) ! compute volume;

            do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;

                if (vvol == 1 .and. iocons == 0) cycle ! fixed inner boundary (or coordinate axis);
                if (vvol == Mvol .and. iocons == 1) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;

                ideriv = 0; id = ideriv
                iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
                call lforce(vvol, iocons, ideriv, Ntz, dBB(1:Ntz, id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag)

            end do ! end of do iocons = 0, 1;
        end do

    end if ! End of if( LocalConstraint )

    deallocate (dBB, stat=astat)
    deallocate (XX, stat=astat) ! spectral constraints; not used;
    deallocate (YY, stat=astat)
    deallocate (length, stat=astat)

2000 continue

end subroutine dfp200
