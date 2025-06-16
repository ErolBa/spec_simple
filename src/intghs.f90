





















module intghs_module

  type intghs_workspace
    REAL, allocatable   :: efmn(:,:)        !< This is efmn.
    REAL, allocatable   :: ofmn(:,:)        !< This is ofmn.
    REAL, allocatable   :: cfmn(:,:)        !<
    REAL, allocatable   :: sfmn(:,:)        !<
    REAL, allocatable   :: evmn(:,:)        !<
    REAL, allocatable   :: odmn(:,:)        !<
    REAL, allocatable   :: ijreal(:,:)      !<
    REAL, allocatable   :: jireal(:,:)      !<
    REAL, allocatable   :: jkreal(:,:)      !<
    REAL, allocatable   :: kjreal(:,:)      !<
    REAL, allocatable   :: Bloweremn(:,:,:) !<
    REAL, allocatable   :: Bloweromn(:,:,:) !<
    REAL, allocatable   :: gBupper(:,:,:)   !<
    REAL, allocatable   :: Blower(:,:,:)    !<
    REAL, allocatable   :: basis(:,:,:,:)   !<
  end type

  TYPE(intghs_workspace) :: wk !< This is an instance of the intghs_workspace type.

end module intghs_module

subroutine intghs_workspace_init(lvol)

  use constants, only : zero
  use inputlist, only : Mpol, Lrad, Wmacros, Wintghs
  use fileunits, only : ounit
  use cputiming, only : Tintghs
  use allglobal, only : Ntz, mn, Iquad, myid, ncpu, cpus, MPI_COMM_SPEC
  use intghs_module

  LOCALS

  INTEGER, INTENT(IN) :: lvol
  INTEGER             :: lquad

  BEGIN(intghs)

  lquad = Iquad(lvol)

  SALLOCATE(wk%gBupper,   (1:Ntz,3,lquad), zero)
  SALLOCATE(wk%Blower,    (1:Ntz,3,lquad), zero)
  SALLOCATE(wk%Bloweremn, (1:mn,3,lquad), zero)
  SALLOCATE(wk%Bloweromn, (1:mn,3,lquad), zero)
  SALLOCATE(wk%efmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%ofmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%evmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%odmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%cfmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%sfmn,      (1:mn,lquad), zero)
  SALLOCATE(wk%ijreal,    (1:mn,lquad), zero)
  SALLOCATE(wk%jkreal,    (1:mn,lquad), zero)
  SALLOCATE(wk%jireal,    (1:mn,lquad), zero)
  SALLOCATE(wk%kjreal,    (1:mn,lquad), zero)
  SALLOCATE(wk%basis,     (0:Lrad(lvol),0:mpol,0:1, lquad), zero)

  RETURN(intghs)

end subroutine intghs_workspace_init

subroutine intghs_workspace_destroy()

  use inputlist, only : Wmacros, Wintghs
  use fileunits, only : ounit
  use cputiming, only : Tintghs
  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC
  use intghs_module

  LOCALS

  BEGIN(intghs)

  DALLOCATE(wk%gBupper)
  DALLOCATE(wk%Blower)
  DALLOCATE(wk%Bloweremn)
  DALLOCATE(wk%Bloweromn)
  DALLOCATE(wk%efmn)
  DALLOCATE(wk%ofmn)
  DALLOCATE(wk%evmn)
  DALLOCATE(wk%odmn)
  DALLOCATE(wk%cfmn)
  DALLOCATE(wk%sfmn)
  DALLOCATE(wk%ijreal)
  DALLOCATE(wk%jkreal)
  DALLOCATE(wk%jireal)
  DALLOCATE(wk%kjreal)
  DALLOCATE(wk%basis)

  RETURN(intghs)

end subroutine intghs_workspace_destroy
