
subroutine mp00ac(Ndof, Xdof, Fdof, Ddof, Ldfjac, iflag) ! argument list is fixed by NAG; ma02aa calls mp00ac through C05PCF;

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
                         Nt, Nz, & ! only required to pass through as arguments to tr00ab;
                         NAdof, &
                         dMA, dMB, dMD, dMG, &
                         Adotx, Ddotx, &
                         NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, & ! preconditioning matrix
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

    integer, parameter :: NB = 4 ! optimal workspace block size for LAPACK:DGECON;

    integer :: lvol, NN, MM, ideriv, lmns, ii, jj, nnz, Lwork

    integer :: idgetrf(0:1), idgetrs(0:1), idgerfs(0:1), idgecon(0:1)

    real(8) :: lmu, dpf, dtf, dpsi(1:2), tpsi(1:2), ppsi(1:2), lcpu, test(2, 2)

    real(8) :: anorm, rcond, ferr(2), berr(2), signfactor

    character :: packorunpack

    integer, allocatable :: ipiv(:), Iwork(:)

    real(8), allocatable :: matrix(:, :), rhs(:, :), LU(:, :)

    real(8), allocatable :: RW(:), RD(:, :)

    real(8), allocatable :: matrixC(:, :)

    integer, parameter :: nrestart = 5 ! do GMRES restart after nrestart iterations
    integer :: maxfil ! bandwidth for ILU subroutines, will be estimated

    integer :: NS, itercount, Nbilut

    real(8), allocatable :: matrixS(:), bilut(:)
    integer, allocatable :: ibilut(:), jbilut(:)

    integer, parameter :: ipar_SIZE = 128
    integer :: ipar(ipar_SIZE), iluierr, RCI_REQUEST, nw, t1, t2, t3
    real(8) :: fpar(ipar_SIZE), v1
    real(8), allocatable :: wk(:)
    integer, allocatable :: jw(:), iperm(:)

    lvol = ivol ! recall that ivol is global;

    if (Lplasmaregion) then

        ; ; lmu = Xdof(1) - xoffset
        ; ; dtf = dtflux(lvol)
        if (Ndof == 2) then; dpf = Xdof(2) - xoffset
        else; dpf = dpflux(lvol)
        end if

    else ! Lvacuumregion;

        ; ; lmu = zero ! restrict attention to strict vacuum field;

        ; ; dtf = Xdof(1) - xoffset
        if (Ndof == 2) then; dpf = Xdof(2) - xoffset
        else; dpf = dpflux(lvol)
        end if

    end if ! end of if( Lplasmaregion ) ;

    dpsi(1:2) = (/dtf, dpf/) ! enclosed poloidal fluxes and their derivatives;
    tpsi(1:2) = (/one, zero/) ! enclosed toroidal fluxes and their derivatives;
    ppsi(1:2) = (/zero, one/) ! enclosed toroidal fluxes and their derivatives;

    diotadxup(0:1, -1:2, lvol) = zero ! rotational-transform, and its derivatives with respect to lmu and dpf, or toroidal current, on the inner/outer interface;
    dItGpdxtp(0:1, -1:2, lvol) = zero ! plasma and linking currents;

    NN = NAdof(lvol) ! shorthand;

    if (allocated(rhs)) deallocate (rhs)
    allocate (rhs(1:NN, 0:2), stat=astat)
    rhs(1:NN, 0:2) = zero
    if (allocated(matrix)) deallocate (matrix)
    allocate (matrix(1:NN, 1:NN), stat=astat)
    matrix(1:NN, 1:NN) = zero

    solution(1:NN, -1:2) = zero ! this is a global array allocated in dforce;

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

    idgetrf(0:1) = 0 ! error flags;
    idgetrs(0:1) = 0 ! error flags;
    idgerfs(0:1) = 0 ! error flags;
    idgecon(0:1) = 0 ! error flags;

    do ideriv = 0, 1 ! loop over derivatives;

        if (iflag == 1 .and. ideriv == 1) cycle ! only need to return function; recall the derivative estimate requires function evaluation;

        if (Lcoordinatesingularity) then

            ; matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN) - lmu*dMD(1:NN, 1:NN)

            ; select case (ideriv)
                ; case (0); rhs(1:NN, 0) = -matmul(dMB(1:NN, 1:2), dpsi(1:2))
                ; case (1) ! construct dMD*solution
                ; ; ; rhs(1:NN, 1) = -matmul(-one*dMD(1:NN, 1:NN), solution(1:NN, 0))
                ; ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2))
                ; end select

        else ! .not.Lcoordinatesingularity;

            if (Lplasmaregion) then

                matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN) - lmu*dMD(1:NN, 1:NN)

                ; select case (ideriv)
                    ; case (0); rhs(1:NN, 0) = -matmul(dMB(1:NN, 1:2), dpsi(1:2))
                    ; case (1) ! construct dMD*solution
                    ; ; ; rhs(1:NN, 1) = -matmul(-one*dMD(1:NN, 1:NN), solution(1:NN, 0))

                    ; ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2))
                    ; end select

            else ! Lvacuumregion ;

                matrix(1:NN, 1:NN) = dMA(1:NN, 1:NN) ! - lmu * dMD(1:NN,1:NN) ;

                select case (ideriv)
                case (0); rhs(1:NN, 0) = -dMG(1:NN) - matmul(dMB(1:NN, 1:2), dpsi(1:2)) ! perhaps there is an lmu term missing here;
                case (1); rhs(1:NN, 1) = -matmul(dMB(1:NN, 1:2), tpsi(1:2)) ! perhaps there is an lmu term missing here;
                    ; ; rhs(1:NN, 2) = -matmul(dMB(1:NN, 1:2), ppsi(1:2)) ! perhaps there is an lmu term missing here;
                end select

            end if ! end of if( Lplasmaregion ) ;

        end if ! end of if( Lcoordinatesingularity ) ;

        select case (ideriv)

        case (0) ! ideriv=0;

            MM = 1
            call DCOPY(NN*NN, matrix, 1, LU, 1) ! BLAS version
            solution(1:NN, 0) = rhs(:, 0)
            call DGETRF(NN, NN, LU, NN, ipiv, idgetrf(ideriv)) ! LU factorization

            anorm = maxval(sum(abs(matrix), 1))
            call DGECON('I', NN, LU, NN, anorm, rcond, RW, Iwork, idgecon(ideriv)) ! estimate the condition number

            call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN, 0), NN, idgetrs(ideriv)) ! sovle linear equation
            call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1, 0), NN, solution(1, 0), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv)) ! refine the solution

        case (1) ! ideriv=1;

            MM = 2
            solution(1:NN, 1:MM) = rhs(:, 1:MM)
            call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN, 1:MM), NN, idgetrs(ideriv))
            call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1, 1), NN, solution(1, 1), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv))

        end select ! ideriv;

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

    end do ! end of do ideriv;

    do ideriv = 0, 2

        if (iflag == 1 .and. ideriv > 0) cycle

        packorunpack = 'U'
        call packab(packorunpack, lvol, NN, solution(1:NN, ideriv), ideriv) ! unpacking; this assigns oAt, oAz through common;

    end do ! do ideriv = 0, 2;

    lBBintegral(lvol) = half*sum(solution(1:NN, 0)*matmul(dMA(1:NN, 1:NN), solution(1:NN, 0))) &
                        + sum(solution(1:NN, 0)*matmul(dMB(1:NN, 1:2), dpsi(1:2))) !

    lABintegral(lvol) = half*sum(solution(1:NN, 0)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 0))) !

    deallocate (matrix, stat=astat)
    deallocate (rhs, stat=astat)
    deallocate (RW, stat=astat)
    deallocate (RD, stat=astat)
    deallocate (LU, stat=astat)
    deallocate (ipiv, stat=astat)
    deallocate (Iwork, stat=astat)

    idgetrf(0:1) = abs(idgetrf(0:1)) + abs(idgetrs(0:1)) + abs(idgerfs(0:1)) + abs(idgecon(0:1))
    if (idgetrf(0) /= 0 .or. idgetrf(1) /= 0) then ! failed to construct Beltrami/vacuum field and/or derivatives;

        ImagneticOK(lvol) = .false. ! set error flag;

        if (iflag == 1) Fdof(1:Ndof) = zero ! provide dummy intent out;
        if (iflag == 2) Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        iflag = -1 ! this value will be returned by C05PCF to ma02aa;

        return

    else

        ImagneticOK(lvol) = .true. ! set error flag; used in dforce;

    end if

    if (YESstellsym) then; lmns = 1 + (mns - 1) ! number of independent degrees of freedom in angle transformation;
    else; lmns = 1 + (mns - 1) + (mns - 1) ! only required for dense, Fourier angle transformation;
    end if

    select case (Lconstraint)

    case (-1) ! Lconstraint=-1;

        if (Lplasmaregion) then

            if (Wtr00ab) then ! compute rotational transform only for diagnostic purposes;
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero ! provide dummy intent out; Lconstraint=-1 indicates no iterations over mu   , dpflux are required;
            Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        else ! Lvacuumregion

            if (Wtr00ab) then ! compute rotational transform only for diagnostic purposes;
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            if (Wcurent) then ! compute enclosed currents    only for diagnostic purposes;
                call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))
                curtor = dItGpdxtp(0, 0, lvol) ! icurrent(0) ! update input variables;
                curpol = dItGpdxtp(1, 0, lvol) ! gcurrent(0)
            end if

            Fdof(1:Ndof) = zero ! provide dummy intent out;Lconstraint=-1 indicates no iterations over dtflux, dpflux are required;
            Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        end if ! end of if( Lplasmaregion) ;

    case (0) ! Lconstraint= 0;

        if (Lplasmaregion) then

            if (Wtr00ab) then ! compute rotational transform only for diagnostic purposes;
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero ! provide dummy intent out; Lconstraint= 0 indicates no iterations over mu, dpflux are required;
            Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        else ! Lvacuumregion

            call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))

            if (iflag == 1) Fdof(1:2) = (/dItGpdxtp(0, 0, lvol) - curtor, dItGpdxtp(1, 0, lvol) - curpol/)
            if (iflag == 2) Ddof(1:2, 1) = (/dItGpdxtp(0, 1, lvol), dItGpdxtp(1, 1, lvol)/)
            if (iflag == 2) Ddof(1:2, 2) = (/dItGpdxtp(0, 2, lvol), dItGpdxtp(1, 2, lvol)/)

        end if ! end of if( Lplasmaregion) ;

    case (1) ! Lconstraint= 1;

        call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol)) ! required for both plasma and vacuum region;

        if (Lplasmaregion) then

            if (Lcoordinatesingularity) then ! Ndof = 1;
                if (iflag == 1) Fdof(1) = diotadxup(1, 0, lvol) - iota(lvol)
                if (iflag == 2) Ddof(1, 1) = diotadxup(1, 1, lvol) ! derivative of outer rotational-transform wrt helicity multiplier
            end if

            if (Ndof == 2) then
                if (iflag == 1) Fdof(1:2) = diotadxup(0:1, 0, lvol) - (/oita(lvol - 1), iota(lvol)/)
                if (iflag == 2) Ddof(1:2, 1) = diotadxup(0:1, 1, lvol)
                if (iflag == 2) Ddof(1:2, 2) = diotadxup(0:1, 2, lvol)
            end if

        else ! Lvacuumregion

            call curent(lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1, -1:2, lvol))

            curtor = dItGpdxtp(0, 0, lvol) ! update input variables; 08 Jun 16;

            if (iflag == 1) Fdof(1:2) = (/diotadxup(0, 0, lvol) - oita(lvol - 1), dItGpdxtp(1, 0, lvol) - curpol/)
            if (iflag == 2) Ddof(1:2, 1) = (/diotadxup(0, 1, lvol), dItGpdxtp(1, 1, lvol)/)
            if (iflag == 2) Ddof(1:2, 2) = (/diotadxup(0, 2, lvol), dItGpdxtp(1, 2, lvol)/)

        end if ! end of if( Lplasmaregion) ;

    case (2)

        if (iflag == 1) Fdof(1) = lABintegral(lvol) - helicity(lvol)

        if (iflag == 2) Ddof(1, 1) = half*sum(solution(1:NN, 1)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 0))) &
                                     + half*sum(solution(1:NN, 0)*matmul(dMD(1:NN, 1:NN), solution(1:NN, 1)))

    case (3)

        if (Lplasmaregion) then

            if (Wtr00ab) then ! compute rotational transform only for diagnostic purposes;
                call tr00ab(lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1, -1:2, lvol))
            end if

            Fdof(1:Ndof) = zero ! provide dummy intent out; no iteration other mu and psip locally
            Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        else ! Lvacuumregion

            Fdof(1:Ndof) = zero ! provide dummy intent out; no iteration other mu and psip locally
            Ddof(1:Ndof, 1:Ndof) = zero ! provide dummy intent out;

        end if ! end of if( Lplasmaregion) ;

    end select ! end of select case( Lconstraint ) ;

    if (Wmp00ac .or. Wma02aa) then ! the following is screen output;

        cput = MPI_WTIME()

        if (Lplasmaregion) then
            select case (iflag)
            case (0); write (ounit, 3000) cput - cpus, myid, lvol, lmu, dpf, iflag ! this is impossible by above logic;
            case (1); write (ounit, 3000) cput - cpus, myid, lvol, lmu, dpf, iflag, Fdof(1:Ndof)
            case (2); write (ounit, 3010) cput - cpus, myid, lvol, lmu, dpf, iflag, Ddof(1:Ndof, 1:Ndof)
            case default
                if (.true.) then
                    write (6, '("mp00ac :      fatal : myid=",i3," ; .true. ; illegal iflag on entry;")') myid
                    call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                    stop "mp00ac : .true. : illegal iflag on entry ;"
                end if
            end select
        else ! Lvacuumregion
            select case (iflag)
            case (0); write (ounit, 3001) cput - cpus, myid, lvol, dtf, dpf, iflag ! this is impossible by above logic;
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

    end if ! end of if( Wmp00ac .or. Wma02aa ) ;

    if (iflag == 1) then ! only in this case is Fdof defined;

        if (sum(abs(Fdof(1:Ndof)))/Ndof < mupftol) then ! satisfactory;

            if (Lplasmaregion) then; mu(lvol) = lmu; ; dpflux(lvol) = dpf

            else; mu(lvol) = zero; dtflux(lvol) = dtf; dpflux(lvol) = dpf

            end if

            iflag = -2 ! return "acceptance" flag through to ma02aa via ifail; early termination;

        end if ! end of if( sum(Fdof) ) ;

    end if ! end of if( iflag.eq.1 ) ;

end subroutine mp00ac
