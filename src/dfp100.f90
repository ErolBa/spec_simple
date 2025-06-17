
subroutine dfp100(Ndofgl, x, Fvec, LComputeDerivatives)

    use constants, only: zero, half, one, two, pi2, pi, mu0

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
                         Lconstraint, curpol

    use cputiming, only: Tdfp100

    use allglobal, only: ncpu, myid, cpus, &
                         ImagneticOK, NAdof, mn, &
                         Mvol, Iquad, &
                         dBdX, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                         IPDt, IPDtdPf, xoffset, dpflux, &
                         IconstraintOK, &
                         DToocc, DToocs, DToosc, DTooss, &
                         TTsscc, TTsscs, TTsssc, TTssss, &
                         TDstcc, TDstcs, TDstsc, TDstss, &
                         TDszcc, TDszcs, TDszsc, TDszss, &
                         DDttcc, DDttcs, DDttsc, DDttss, &
                         DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                         DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                         dMA, dMB, dMD, dMG, MBpsi, solution, &
                         Nt, Nz, Lsavedguvij, guvijsave, izbs

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer :: vvol, Ndofgl, iflag, cpu_send_one, cpu_send_two, ll, NN, ideriv, iocons
    integer :: request1, request2
    real(8) :: Fvec(1:Ndofgl), x(1:Mvol - 1), Bt00(1:Mvol, 0:1, -1:2), ldItGp(0:1, -1:2)
    logical :: LComputeDerivatives
    integer :: deriv, Lcurvature

    dpflux(2:Mvol) = x - xoffset

    do vvol = 1, Mvol

        if (Igeometry == 1 .or. vvol > 1) then; Lcoordinatesingularity = .false.
        else; Lcoordinatesingularity = .true.
        end if
        ImagneticOK(vvol) = .false.

        NN = NAdof(vvol)
        ll = Lrad(vvol)

        call allocate_geometry_matrices(vvol, LcomputeDerivatives)
        call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
        call intghs_workspace_init(vvol)

        dBdX%L = .false.

        ideriv = 0; Lcurvature = 1
        call compute_guvijsave(Iquad(vvol), vvol, ideriv, Lcurvature)
        Lsavedguvij = .true.

        call ma00aa(Iquad(vvol), mn, vvol, ll)
        call matrix(vvol, mn, ll)

        call ma02aa(vvol, NN)

        Lsavedguvij = .false.

        call intghs_workspace_destroy()
        call deallocate_Beltrami_matrices(LcomputeDerivatives)
        call deallocate_geometry_matrices(LcomputeDerivatives)

        if (Lconstraint == 3) then
            do iocons = 0, 1
                if ((Lcoordinatesingularity .and. iocons == 0) .or. (Lvacuumregion .and. iocons == 1)) cycle

                ideriv = 0
                call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)

                ideriv = 2
                call lbpol(vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)
            end do

            if (Lvacuumregion) then

                ideriv = 1
                iocons = 0
                call lbpol(Mvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons)

                iflag = 2
                call curent(Mvol, mn, Nt, Nz, iflag, ldItGp(0:1, -1:2))
            end if
        end if
    end do

    if (.not. LocalConstraint) then

        select case (Lconstraint)

        case (3)

            do vvol = 1, Mvol - 1

                IPDt(vvol) = pi2*(Bt00(vvol + 1, 0, 0) - Bt00(vvol, 1, 0))

                IPDtdPf(vvol, vvol) = pi2*Bt00(vvol + 1, 0, 2)
                if (vvol /= 1) IPDtdPf(vvol, vvol - 1) = -pi2*Bt00(vvol, 1, 2)
            end do

            if (myid == 0) then
                Fvec(1:Mvol - 1) = IPDt(1:Mvol - 1) - Isurf(1:Mvol - 1)
            end if

        case default
            if (.true.) then
                write (6, '("dfp100 :      fatal : myid=",i3," ; .true. ; Unaccepted value for Lconstraint;")') myid

                stop "dfp100 : .true. : Unaccepted value for Lconstraint ;"
            end if
        end select
    end if

end subroutine dfp100
