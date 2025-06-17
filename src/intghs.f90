
module intghs_module

    type intghs_workspace
        real(8), allocatable :: efmn(:, :)
        real(8), allocatable :: ofmn(:, :)
        real(8), allocatable :: cfmn(:, :)
        real(8), allocatable :: sfmn(:, :)
        real(8), allocatable :: evmn(:, :)
        real(8), allocatable :: odmn(:, :)
        real(8), allocatable :: ijreal(:, :)
        real(8), allocatable :: jireal(:, :)
        real(8), allocatable :: jkreal(:, :)
        real(8), allocatable :: kjreal(:, :)
        real(8), allocatable :: Bloweremn(:, :, :)
        real(8), allocatable :: Bloweromn(:, :, :)
        real(8), allocatable :: gBupper(:, :, :)
        real(8), allocatable :: Blower(:, :, :)
        real(8), allocatable :: basis(:, :, :, :)
    end type

    type(intghs_workspace) :: wk

end module intghs_module

subroutine intghs_workspace_init(lvol)

    use constants, only: zero
    use inputlist, only: Mpol, Lrad, Wmacros, Wintghs
    use fileunits, only: ounit

    use allglobal, only: Ntz, mn, Iquad, myid, ncpu, cpus
    use intghs_module

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(IN) :: lvol
    integer :: lquad

    lquad = Iquad(lvol)

    if (allocated(wk%gBupper)) deallocate (wk%gBupper)
    allocate (wk%gBupper(1:Ntz, 3, lquad), stat=astat)
    wk%gBupper(1:Ntz, 3, lquad) = zero
    if (allocated(wk%Blower)) deallocate (wk%Blower)
    allocate (wk%Blower(1:Ntz, 3, lquad), stat=astat)
    wk%Blower(1:Ntz, 3, lquad) = zero
    if (allocated(wk%Bloweremn)) deallocate (wk%Bloweremn)
    allocate (wk%Bloweremn(1:mn, 3, lquad), stat=astat)
    wk%Bloweremn(1:mn, 3, lquad) = zero
    if (allocated(wk%Bloweromn)) deallocate (wk%Bloweromn)
    allocate (wk%Bloweromn(1:mn, 3, lquad), stat=astat)
    wk%Bloweromn(1:mn, 3, lquad) = zero
    if (allocated(wk%efmn)) deallocate (wk%efmn)
    allocate (wk%efmn(1:mn, lquad), stat=astat)
    wk%efmn(1:mn, lquad) = zero
    if (allocated(wk%ofmn)) deallocate (wk%ofmn)
    allocate (wk%ofmn(1:mn, lquad), stat=astat)
    wk%ofmn(1:mn, lquad) = zero
    if (allocated(wk%evmn)) deallocate (wk%evmn)
    allocate (wk%evmn(1:mn, lquad), stat=astat)
    wk%evmn(1:mn, lquad) = zero
    if (allocated(wk%odmn)) deallocate (wk%odmn)
    allocate (wk%odmn(1:mn, lquad), stat=astat)
    wk%odmn(1:mn, lquad) = zero
    if (allocated(wk%cfmn)) deallocate (wk%cfmn)
    allocate (wk%cfmn(1:mn, lquad), stat=astat)
    wk%cfmn(1:mn, lquad) = zero
    if (allocated(wk%sfmn)) deallocate (wk%sfmn)
    allocate (wk%sfmn(1:mn, lquad), stat=astat)
    wk%sfmn(1:mn, lquad) = zero
    if (allocated(wk%ijreal)) deallocate (wk%ijreal)
    allocate (wk%ijreal(1:mn, lquad), stat=astat)
    wk%ijreal(1:mn, lquad) = zero
    if (allocated(wk%jkreal)) deallocate (wk%jkreal)
    allocate (wk%jkreal(1:mn, lquad), stat=astat)
    wk%jkreal(1:mn, lquad) = zero
    if (allocated(wk%jireal)) deallocate (wk%jireal)
    allocate (wk%jireal(1:mn, lquad), stat=astat)
    wk%jireal(1:mn, lquad) = zero
    if (allocated(wk%kjreal)) deallocate (wk%kjreal)
    allocate (wk%kjreal(1:mn, lquad), stat=astat)
    wk%kjreal(1:mn, lquad) = zero
    if (allocated(wk%basis)) deallocate (wk%basis)
    allocate (wk%basis(0:Lrad(lvol), 0:mpol, 0:1, lquad), stat=astat)
    wk%basis(0:Lrad(lvol), 0:mpol, 0:1, lquad) = zero

end subroutine intghs_workspace_init

subroutine intghs_workspace_destroy()

    use inputlist, only: Wmacros, Wintghs
    use fileunits, only: ounit

    use allglobal, only: myid, ncpu, cpus
    use intghs_module

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    deallocate (wk%gBupper, stat=astat)
    deallocate (wk%Blower, stat=astat)
    deallocate (wk%Bloweremn, stat=astat)
    deallocate (wk%Bloweromn, stat=astat)
    deallocate (wk%efmn, stat=astat)
    deallocate (wk%ofmn, stat=astat)
    deallocate (wk%evmn, stat=astat)
    deallocate (wk%odmn, stat=astat)
    deallocate (wk%cfmn, stat=astat)
    deallocate (wk%sfmn, stat=astat)
    deallocate (wk%ijreal, stat=astat)
    deallocate (wk%jkreal, stat=astat)
    deallocate (wk%jireal, stat=astat)
    deallocate (wk%kjreal, stat=astat)
    deallocate (wk%basis, stat=astat)

end subroutine intghs_workspace_destroy
