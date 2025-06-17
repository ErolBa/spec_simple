
subroutine ma02aa(lvol, NN)

    use constants, only: zero, half, one, ten

    use numerical, only: vsmall, small

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wma02aa, &
                         Lconstraint, mu, helicity, &
                         mupftol, mupfits, Lrad, Lcheck

    use cputiming

    use allglobal, only: ncpu, myid, cpus, MPI_COMM_SPEC, &
                         Mvol, mn, im, in, &
                         LBlinear, LBnewton, LBsequad, &
                         dMA, dMB, dMD, solution, &
                         MBpsi, Ate, &
                         ImagneticOK, &
                         lBBintegral, lABintegral, &
                         ivol, Nfielddof, &
                         dtflux, dpflux, &
                         xoffset, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: lvol, NN

    integer :: ideriv
    real(8) :: tol, dpsi(1:2), lastcpu
    character :: packorunpack

    integer :: Nxdof, Ndof, Ldfjac, iflag, maxfev, mode, LRR, nfev, njev, nprint, ihybrj
    real(8) :: Xdof(1:2), Fdof(1:2), Ddof(1:2, 1:2), oDdof(1:2, 1:2)
    real(8) :: factor, diag(1:2), RR(1:2*(2 + 1)/2), QTF(1:2), wk(1:2, 1:4)

    integer :: irevcm

    integer :: pNN

    real(8) :: xi(0:NN), Fxi(0:NN), xo(0:NN), Mxi(1:NN)

    external :: mp00ac

    integer :: ihybrj1, Ldfmuaa, lengthwork
    real(8) :: NewtonError
    real(8), allocatable :: DFxi(:, :), work(:)
    external :: df00ab

    integer :: NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, iterations, LIWk, LRWk, ie04uff
    integer, allocatable :: Istate(:), NEEDC(:), IWk(:)
    real(8) :: objectivefunction
    real(8), allocatable :: LinearConstraintMatrix(:, :), LowerBound(:), UpperBound(:)
    real(8), allocatable :: constraintfunction(:), constraintgradient(:, :), multipliers(:), objectivegradient(:), RS(:, :), RWk(:)
    character :: optionalparameter*33

    ivol = lvol

    if (LBsequad) then
        lastcpu = MPI_WTIME()

        NLinearConstraints = 0

        NNonLinearConstraints = 1

        LDA = max(1, NLinearConstraints)

        LDCJ = max(1, NNonLinearConstraints)

        LDR = NN

        if (allocated(LinearConstraintMatrix)) deallocate (LinearConstraintMatrix)
        allocate (LinearConstraintMatrix(1:LDA, 1:1), stat=astat)
        LinearConstraintMatrix(1:LDA, 1:1) = zero

        if (allocated(LowerBound)) deallocate (LowerBound)
        allocate (LowerBound(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        LowerBound(1:NN + NLinearConstraints + NNonLinearConstraints) = zero
        if (allocated(UpperBound)) deallocate (UpperBound)
        allocate (UpperBound(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        UpperBound(1:NN + NLinearConstraints + NNonLinearConstraints) = zero

        LowerBound(1:NN) = -1.0e+21
        UpperBound(1:NN) = +1.0e+21
        LowerBound(NN + 1:NN + NLinearConstraints) = -1.0e+21
        UpperBound(NN + 1:NN + NLinearConstraints) = +1.0e+21
        LowerBound(NN + NLinearConstraints + 1:NN + NLinearConstraints + NNonLinearConstraints) = helicity(lvol)
        UpperBound(NN + NLinearConstraints + 1:NN + NLinearConstraints + NNonLinearConstraints) = helicity(lvol)

        iterations = 0

        if (allocated(Istate)) deallocate (Istate)
        allocate (Istate(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        Istate(1:NN + NLinearConstraints + NNonLinearConstraints) = 0

        if (allocated(constraintfunction)) deallocate (constraintfunction)
        allocate (constraintfunction(1:NNonLinearConstraints), stat=astat)
        constraintfunction(1:NNonLinearConstraints) = zero

        if (allocated(constraintgradient)) deallocate (constraintgradient)
        allocate (constraintgradient(1:LDCJ, 1:NN), stat=astat)
        constraintgradient(1:LDCJ, 1:NN) = zero

        if (allocated(multipliers)) deallocate (multipliers)
        allocate (multipliers(1:NN + NLinearConstraints + NNonLinearConstraints), stat=astat)
        multipliers(1:NN + NLinearConstraints + NNonLinearConstraints) = zero

        objectivefunction = zero

        if (allocated(objectivegradient)) deallocate (objectivegradient)
        allocate (objectivegradient(1:NN), stat=astat)
        objectivegradient(1:NN) = zero

        if (allocated(RS)) deallocate (RS)
        allocate (RS(1:LDR, 1:NN), stat=astat)
        RS(1:LDR, 1:NN) = zero
        ideriv = 0; dpsi(1:2) = (/dtflux(lvol), dpflux(lvol)/)

        packorunpack = 'P'

        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)

        if (allocated(NEEDC)) deallocate (NEEDC)
        allocate (NEEDC(1:NNonLinearConstraints), stat=astat)
        NEEDC(1:NNonLinearConstraints) = 0

        LIWk = 3*NN + NLinearConstraints + 2*NNonLinearConstraints
        if (allocated(IWk)) deallocate (IWk)
        allocate (IWk(1:LIWk), stat=astat)
        IWk(1:LIWk) = 0

        LRWk = 2*NN**2 + NN*NLinearConstraints + 2*NN*NNonLinearConstraints + 21*NN + 11*NLinearConstraints + 22*NNonLinearConstraints + 1
        if (allocated(RWk)) deallocate (RWk)
        allocate (RWk(1:LRWk), stat=astat)
        RWk(1:LRWk) = zero

        irevcm = 0; ie04uff = 1

        MBpsi(1:NN) = matmul(dMB(1:NN, 1:2), dpsi(1:2))

        deallocate (RWk, stat=astat)
        deallocate (IWk, stat=astat)
        deallocate (NEEDC, stat=astat)
        deallocate (RS, stat=astat)
        deallocate (objectivegradient, stat=astat)
        deallocate (multipliers, stat=astat)
        deallocate (constraintgradient, stat=astat)
        deallocate (constraintfunction, stat=astat)
        deallocate (Istate, stat=astat)
        deallocate (LowerBound, stat=astat)
        deallocate (UpperBound, stat=astat)
        deallocate (LinearConstraintMatrix, stat=astat)

        packorunpack = 'U'
        call packab(packorunpack, lvol, NN, xi(1:NN), ideriv)

        lBBintegral(lvol) = half*sum(xi(1:NN)*matmul(dMA(1:NN, 1:NN), xi(1:NN))) + sum(xi(1:NN)*MBpsi(1:NN))
        lABintegral(lvol) = half*sum(xi(1:NN)*matmul(dMD(1:NN, 1:NN), xi(1:NN)))

        solution(1:NN, 0) = xi(1:NN)

    end if

    if (LBlinear) then

        lastcpu = MPI_WTIME()

        if (Lplasmaregion) then

            Xdof(1:2) = xoffset + (/mu(lvol), dpflux(lvol)/)

            select case (Lconstraint)
            case (-1); ; Nxdof = 0
                ; ; iflag = 1
            case (0); ; Nxdof = 0
                ; ; iflag = 1
            case (1); if (Lcoordinatesingularity) then; Nxdof = 1
                    ; else; Nxdof = 2
                    ; end if
            case (2); Nxdof = 1
            case (3); if (Lcoordinatesingularity) then; Nxdof = 0
                    ; else; Nxdof = 0
                    ; end if
                ; ; iflag = 2
            end select

        else

            Xdof(1:2) = xoffset + (/dtflux(lvol), dpflux(lvol)/)

            select case (Lconstraint)
            case (-1); ; Nxdof = 0
            case (0); ; Nxdof = 2
            case (1); ; Nxdof = 2
            case (2); ; Nxdof = 2
            case (3); ; Nxdof = 0
                ; iflag = 2
            end select

        end if

        select case (Nxdof)

        case (0)

            ; ; Ndof = 1; Ldfjac = Ndof; nfev = 1; njev = 0; ihybrj = 1; 
            call mp00ac(Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)

            helicity(lvol) = lABintegral(lvol)

        case (1:2)

            ; ; Ndof = Nxdof; Ldfjac = Ndof; nfev = 0; njev = 0; ihybrj = 0; 
            tol = mupftol; LRR = Ndof*(Ndof + 1)/2; mode = 0; diag(1:2) = zero; factor = one; maxfev = mupfits; nprint = 0

            if (Ndof > 2) then
                write (6, '("ma02aa :      fatal : myid=",i3," ; Ndof.gt.2 ; illegal;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "ma02aa : Ndof.gt.2 : illegal ;"
            end if

            call hybrj2(mp00ac, Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, tol, &
                        maxfev, diag(1:Ndof), mode, factor, nprint, ihybrj, nfev, njev, RR(1:LRR), LRR, QTF(1:Ndof), &
                        WK(1:Ndof, 1), WK(1:Ndof, 2), WK(1:Ndof, 3), WK(1:Ndof, 4))

            if (Lplasmaregion) then

                select case (ihybrj)
                case (0:); mu(lvol) = Xdof(1) - xoffset
                    ; ; dpflux(lvol) = Xdof(2) - xoffset
                case (:-1); Xdof(1) = mu(lvol) + xoffset
                    ; ; Xdof(2) = dpflux(lvol) + xoffset
                end select

            else

                select case (ihybrj)
                case (0:); dtflux(lvol) = Xdof(1) - xoffset
                    ; ; dpflux(lvol) = Xdof(2) - xoffset
                case (:-1); Xdof(1) = dtflux(lvol) + xoffset
                    ; ; Xdof(2) = dpflux(lvol) + xoffset
                end select

            end if

            if (Lconstraint /= 2) helicity(lvol) = lABintegral(lvol)

            if (Lconstraint == 1 .or. Lconstraint == 3 .or. (Lvacuumregion .and. Lconstraint == 0)) then

                iflag = 2; Ldfjac = Ndof

                call mp00ac(Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof), Ldfjac, iflag)

            end if

        end select

        cput = MPI_WTIME()

        select case (ihybrj)
        case (1)
            if (Wma02aa) write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "success         ", Fdof(1:Ndof)
        case (-2)
            if (Wma02aa) write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "|F| < mupftol   ", Fdof(1:Ndof)
        case (-1)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "Beltrami fail   ", Fdof(1:Ndof)
        case (0)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "input error     ", Fdof(1:Ndof)
        case (2)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "consider restart", Fdof(1:Ndof)
        case (3)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "xtol too small  ", Fdof(1:Ndof)
        case (4:5)
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "bad progress    ", Fdof(1:Ndof)
        case default
            ; write (ounit, 1040) cput - cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput - lastcpu, "illegal ifail   ", Fdof(1:Ndof)
            if (.true.) then
                write (6, '("ma02aa :      fatal : myid=",i3," ; .true. ; illegal ifail returned by hybrj;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "ma02aa : .true. : illegal ifail returned by hybrj ;"
            end if
        end select

    end if

1010 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; SQP    : ie04uff=", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ":, a36)
1020 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; Newton : ihybrj1=", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; " &
           "error="es7.0" ; ":, a18)
1040 format("ma02aa : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; Linear : ihybrj =", i3, " hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; " &
           :, a16" ; F="2es08.0)

end subroutine ma02aa

