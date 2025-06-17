
subroutine mp00ac(Ndof, Xdof, Fdof, Ddof, Ldfjac, iflag)

    use constants, only: zero, half, one

    use numerical, only: small, machprec

    use fileunits, only: ounit

    use inputlist, only: Wmacros, Wmp00ac, Wtr00ab, Wcurent, Wma02aa, &
                         mu, helicity, iota, oita, curtor, curpol, Lrad, Ntor, &
                         Lconstraint, mupftol, &
                         NiterGMRES, epsGMRES, LGMRESprec, epsILU

    use cputiming, only: Tmp00ac

    use allglobal, only: myid, ncpu, cpus, ivol, MPI_COMM_SPEC, &
                         YESstellsym, NOTstellsym, &
                         Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                         mn, im, in, mns, &
                         Nt, Nz, &
                         NAdof, &
                         dMA, dMB, dMD, dMG, &
                         Adotx, Ddotx, &
                         NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, &
                         solution, GMRESlastsolution, &
                         dtflux, dpflux, &
                         diotadxup, dItGpdxtp, &
                         lBBintegral, lABintegral, &
                         xoffset, &
                         ImagneticOK, &
                         Ate, Aze, Ato, Azo, Mvol, Iquad, &
                         GMRESlastsolution, ext

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: Ndof, Ldfjac
    real(8), intent(in) :: Xdof(1:Ndof)
    real(8) :: Fdof(1:Ndof), Ddof(1:Ldfjac, 1:Ndof)
    integer :: iflag

    integer, parameter :: NB = 4

    integer :: lvol, NN, MM, ideriv, lmns, ii, jj, nnz, Lwork

    integer :: idgetrf(0:1), idgetrs(0:1), idgerfs(0:1), idgecon(0:1)

    real(8) :: lmu, dpf, dtf, dpsi(1:2), tpsi(1:2), ppsi(1:2), lcpu, test(2, 2)

    real(8) :: anorm, rcond, ferr(2), berr(2), signfactor

    character :: packorunpack

    integer, allocatable :: ipiv(:), Iwork(:)

    real(8), allocatable :: matrix(:, :), rhs(:, :), LU(:, :)

    real(8), allocatable :: RW(:), RD(:, :)

    real(8), allocatable :: matrixC(:, :)

    integer, parameter :: nrestart = 5
    integer :: maxfil

    integer :: NS, itercount, Nbilut

    real(8), allocatable :: matrixS(:), bilut(:)
    integer, allocatable :: ibilut(:), jbilut(:)

    integer, parameter :: ipar_SIZE = 128
    integer :: ipar(ipar_SIZE), iluierr, RCI_REQUEST, nw, t1, t2, t3
    real(8) :: fpar(ipar_SIZE), v1
    real(8), allocatable :: wk(:)
    integer, allocatable :: jw(:), iperm(:)

    lvol = ivol

    if (Lplasmaregion) then

        ; ; lmu = Xdof(1) - xoffset
        ; ; dtf = dtflux(lvol)
        if (Ndof == 2) then; dpf = Xdof(2) - xoffset
        else; dpf = dpflux(lvol)
        end if

    else

        ; ; lmu = zero

        ; ; dtf = Xdof(1) - xoffset
        if (Ndof == 2) then; dpf = Xdof(2) - xoffset
        else; dpf = dpflux(lvol)
        end if

    end if

    dpsi(1:2) = (/dtf, dpf/)
    tpsi(1:2) = (/one, zero/)
    ppsi(1:2) = (/zero, one/)

    diotadxup(0:1, -1:2, lvol) = zero
    dItGpdxtp(0:1, -1:2, lvol) = zero

    NN = NAdof(lvol)

    if (allocated(rhs)) deallocate (rhs)
    allocate (rhs(1:NN, 0:2), stat=astat)
    rhs(1:NN, 0:2) = zero
    if (allocated(matrix)) deallocate (matrix)
    allocate (matrix(1:NN, 1:NN), stat=astat)
    matrix(1:NN, 1:NN) = zero

    solution(1:NN, -1:2) = zero

    Lwork = NB*NN

    if (allocated(RW)) deallocate (RW)
    allocate (RW(1:Lwork), stat=astat)
    RW(1:Lwork) = zero
    if (allocated(RD)) deallocate (RD)
    allocate (RD(1:NN, 0:2), stat=astat)
    RD(1:NN, 0:2) = zero
    if (allocated(LU)) deallocate (LU)
    allocate (LU(1:NN, 1:NN), stat=astat)
    LU(1:NN, 1:NN) = zero
    if (allocated(ipiv)) deallocate (ipiv)
    allocate (ipiv(1:NN), stat=astat)
    ipiv(1:NN) = 0
    if (allocated(Iwork)) deallocate (Iwork)
    allocate (Iwork(1:NN), stat=astat)
    Iwork(1:NN) = 0
    nw = (NN + 3)*(nrestart + 2) + (nrestart + 1)*nrestart
    if (allocated(wk)) deallocate (wk)
    allocate (wk(1:nw), stat=astat)
    wk(1:nw) = zero
    if (allocated(jw)) deallocate (jw)
    allocate (jw(1:2*NN), stat=astat)
    jw(1:2*NN) = 0
    if (allocated(iperm)) deallocate (iperm)
    allocate (iperm(1:2*NN), stat=astat)
    iperm(1:2*NN) = 0

    idgetrf(0:1) = 0
    idgetrs(0:1) = 0
    idgerfs(0:1) = 0
    idgecon(0:1) = 0

    do ideriv = 0, 1

        if (iflag == 1 .and. ideriv == 1) cycle

        if (Lcoordinatesingularity) then

            ; matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN) - lmu*dMD(1:NN, 1:NN)

            ; select case (ideriv)
                ; case (0); rhs(1:NN, 0) = -matmul(dMB(1:NN, 1:2), dpsi(1:2))
                ; case (1)
                ; ; ; rhs(1:NN, 1) = -matmul(-one*dMD(1:NN, 1:NN), solution(1:NN, 0))
                ; ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2))
                ; end select

        else

            if (Lplasmaregion) then

                matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN) - lmu*dMD(1:NN, 1:NN)

                ; select case (ideriv)
                    ; case (0); rhs(1:NN, 0) = -matmul(dMB(1:NN, 1:2), dpsi(1:2))
                    ; case (1)
                    ; ; ; rhs(1:NN, 1) = -matmul(-one*dMD(1:NN, 1:NN), solution(1:NN, 0))

                    ; ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2))
                    ; end select

            else

                matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN)

                select case (ideriv)
                case (0); rhs(1:NN, 0) = -dMG(1:NN) - matmul(dMB(1:NN, 1:2), dpsi(1:2))
                case (1); rhs(1:NN, 1) = -matmul(dMB(1:NN, 1:2), tpsi(1:2))
                    ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2))
                end select

            end if

        end if

        select case (ideriv)

        case (0)

            MM = 1
            call DCOPY(NN*NN, matrix, 1, LU, 1)
            solution(1:NN, 0) = rhs(:, 0)
            call DGETRF(NN, NN, LU, NN, ipiv, idgetrf(ideriv))

            anorm = maxval(sum(abs(matrix), 1))
            call DGECON('I', NN, LU, NN, anorm, rcond, RW, Iwork, idgecon(ideriv))

            call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN, 0), NN, idgetrs(ideriv))
            call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1, 0), NN, solution(1, 0), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv))

        case (1)

            MM = 2
            solution(1:NN, 1:MM) = rhs(:, 1:MM)
            call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN, 1:MM), NN, idgetrs(ideriv))
            call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1, 1), NN, solution(1, 1), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv))

        end select

        cput = MPI_WTIME()

        if (idgetrf(ideriv) == 0 .and. idgetrs(ideriv) == 0 .and. idgerfs(ideriv) == 0 .and. rcond >= machprec) then
            if (Wmp00ac) write (ounit, 1010) cput - cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "success ;         ", cput - lcpu
        elseif (idgetrf(ideriv) < 0 .or. idgetrs(ideriv) < 0 .or. idgerfs(ideriv) < 0) then
            ; write (ounit, 1010) cput - cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "input error ;     "
        elseif (idgetrf(ideriv) > 0) then
            ; write (ounit, 1010) cput - cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "singular ;        "
        elseif (rcond <= machprec) then
            ; write (ounit, 1010) cput - cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "ill conditioned ; "
        else
            ; write (ounit, 1010) cput - cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "invalid error ; "
        end if

1010    format("mp00ac : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; ideriv="i2" ; "a23"=", i3, ' ', i3, ' ', i3, " ; "a34, :" time=", f10.2, " ;")

    end do

    do ideriv = 0, 2

        if (iflag == 1 .and. ideriv > 0) cycle

        packorunpack = 'U'
        call packab(packorunpack, lvol, NN, solution(1:NN, ideriv), ideriv)

    end do

    lBBintegral(lvol) = half*sum(solution(1:NN, 0)*matmul(dMA(1:NN, 1:NN), solution(1:NN, 0))) &
                        + sum(solution(1:NN, 0)*matmul(dMB(1:NN, 1:2), dpsi(1:2)))

    lABintegral(lvol) = half*sum(solution(1:NN, 0)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 0)))

    deallocate (matrix, stat=astat)
    deallocate (rhs, stat=astat)
    deallocate (RW, stat=astat)
    deallocate (RD, stat=astat)
    deallocate (LU, stat=astat)
    deallocate (ipiv, stat=astat)
    deallocate (Iwork, stat=astat)

    idgetrf(0:1) = abs(idgetrf(0:1)) + abs(idgetrs(0:1)) + abs(idgerfs(0:1)) + abs(idgecon(0:1))
    if (idgetrf(0) /= 0 .or. idgetrf(1) /= 0) then

        ImagneticOK(lvol) = .false.

        if (iflag == 1) Fdof(1:Ndof) = zero
        if (iflag == 2) Ddof(1:Ndof, 1:Ndof) = zero

        iflag = -1

        return

    else

        ImagneticOK(lvol) = .true.

    end if

    if (YESstellsym) then; lmns = 1 + (mns - 1)
    else; lmns = 1 + (mns - 1) + (mns - 1)
    end if

    select case (Lconstraint)

    case (-1)

        if (Lplasmaregion) then

            if (Wtr00ab) then
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero
            Ddof(1:Ndof, 1:Ndof) = zero

        else

            if (Wtr00ab) then
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            if (Wcurent) then
                call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))
                curtor = dItGpdxtp(0, 0, lvol)
                curpol = dItGpdxtp(1, 0, lvol)
            end if

            Fdof(1:Ndof) = zero
            Ddof(1:Ndof, 1:Ndof) = zero

        end if

    case (0)

        if (Lplasmaregion) then

            if (Wtr00ab) then
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero
            Ddof(1:Ndof, 1:Ndof) = zero

        else

            call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))

            if (iflag == 1) Fdof(1:2) = (/dItGpdxtp(0, 0, lvol) - curtor, dItGpdxtp(1, 0, lvol) - curpol/)
            if (iflag == 2) Ddof(1:2, 1) = (/dItGpdxtp(0, 1, lvol), dItGpdxtp(1, 1, lvol)/)
            if (iflag == 2) Ddof(1:2, 2) = (/dItGpdxtp(0, 2, lvol), dItGpdxtp(1, 2, lvol)/)

        end if

    case (1)

        call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))

        if (Lplasmaregion) then

            if (Lcoordinatesingularity) then
                if (iflag == 1) Fdof(1) = diotadxup(1, 0, lvol) - iota(lvol)
                if (iflag == 2) Ddof(1, 1) = diotadxup(1, 1, lvol)
            end if

            if (Ndof == 2) then
                if (iflag == 1) Fdof(1:2) = diotadxup(0:1, 0, lvol) - (/oita(lvol - 1), iota(lvol)/)
                if (iflag == 2) Ddof(1:2, 1) = diotadxup(0:1, 1, lvol)
                if (iflag == 2) Ddof(1:2, 2) = diotadxup(0:1, 2, lvol)
            end if

        else

            call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))

            curtor = dItGpdxtp(0, 0, lvol)

            if (iflag == 1) Fdof(1:2) = (/diotadxup(0, 0, lvol) - oita(lvol - 1), dItGpdxtp(1, 0, lvol) - curpol/)
            if (iflag == 2) Ddof(1:2, 1) = (/diotadxup(0, 1, lvol), dItGpdxtp(1, 1, lvol)/)
            if (iflag == 2) Ddof(1:2, 2) = (/diotadxup(0, 2, lvol), dItGpdxtp(1, 2, lvol)/)

        end if

    case (2)

        if (iflag == 1) Fdof(1) = lABintegral(lvol) - helicity(lvol)

        if (iflag == 2) Ddof(1, 1) = half*sum(solution(1:NN, 1)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 0))) &
                                     + half*sum(solution(1:NN, 0)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 1)))

    case (3)

        if (Lplasmaregion) then

            if (Wtr00ab) then
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero
            Ddof(1:Ndof, 1:Ndof) = zero

        else

            Fdof(1:Ndof) = zero
            Ddof(1:Ndof, 1:Ndof) = zero

        end if

    end select

    if (Wmp00ac .or. Wma02aa) then

        cput = MPI_WTIME()

        if (Lplasmaregion) then
            select case (iflag)
            case (0); write (ounit, 3000) cput - cpus, myid, lvol, lmu, dpf, iflag
            case (1); write (ounit, 3000) cput - cpus, myid, lvol, lmu, dpf, iflag, Fdof(1:Ndof)
            case (2); write (ounit, 3010) cput - cpus, myid, lvol, lmu, dpf, iflag, Ddof(1:Ndof, 1:Ndof)
            case default
                if (.true.) then
                    write (6, '("mp00ac :      fatal : myid=",i3," ; .true. ; illegal iflag on entry;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "mp00ac : .true. : illegal iflag on entry ;"
                end if
            end select
        else
            select case (iflag)
            case (0); write (ounit, 3001) cput - cpus, myid, lvol, dtf, dpf, iflag
            case (1); write (ounit, 3001) cput - cpus, myid, lvol, dtf, dpf, iflag, Fdof(1:Ndof)
            case (2); write (ounit, 3011) cput - cpus, myid, lvol, dtf, dpf, iflag, Ddof(1:Ndof, 1:Ndof)
            case default
                if (.true.) then
                    write (6, '("mp00ac :      fatal : myid=",i3," ; .true. ; illegal iflag on entry;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "mp00ac : .true. : illegal iflag on entry ;"
                end if
            end select
        end if

3000    format("mp00ac : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3010    format("mp00ac : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")
3001    format("mp00ac : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3011    format("mp00ac : ", f10.2, " : myid=", i3, " ; lvol=", i3, " ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")

    end if

    if (iflag == 1) then

        if (sum(abs(Fdof(1:Ndof)))/Ndof < mupftol) then

            if (Lplasmaregion) then; mu(lvol) = lmu; ; dpflux(lvol) = dpf

            else; mu(lvol) = zero; dtflux(lvol) = dtf; dpflux(lvol) = dpf

            end if

            iflag = -2

        end if

    end if

end subroutine mp00ac
