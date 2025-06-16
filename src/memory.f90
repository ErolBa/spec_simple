!> \file
!> \brief memory management module

!> \brief allocate Beltrami matrices
!>
!> @param vvol
!> @param LcomputeDerivatives
subroutine allocate_Beltrami_matrices(vvol, LcomputeDerivatives)

  use fileunits

  use inputlist, only: Wmemory, Wmacros

  use allglobal

  use cputiming

  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  REAL      :: cput, cpui, cpuo=0

  INTEGER, intent(in) :: vvol
  LOGICAL, intent(in) :: LcomputeDerivatives
  INTEGER     :: NN

  

  NN = NAdof(vvol) ! shorthand;

  SALLOCATE( dMA, (0:NN,0:NN), zero ) ! required for both plasma region and vacuum region;
  SALLOCATE( dMD, (0:NN,0:NN), zero )

  ! we will need the rest even with or without matrix-free
  SALLOCATE( dMB, (0:NN,0: 2), zero )
  SALLOCATE( dMG, (0:NN     ), zero )

  SALLOCATE( solution, (1:NN,-1:2), zero ) ! this will contain the vector potential from the linear solver and its derivatives;

  SALLOCATE( MBpsi, (1:NN), zero )



end subroutine allocate_Beltrami_matrices



!> \brief deallocate Beltrami matrices
!>
!> @param LcomputeDerivatives
subroutine deallocate_Beltrami_matrices(LcomputeDerivatives)

  use fileunits

  use inputlist, only: Wmemory, Wmacros

  use allglobal

  use cputiming

  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  REAL      :: cput, cpui, cpuo=0

  LOGICAL, intent(in) :: LcomputeDerivatives

  

  deallocate(dMA,stat=astat)
  deallocate(dMD,stat=astat)

  deallocate(dMB,stat=astat)

  deallocate(dMG,stat=astat)

  deallocate(solution,stat=astat)

  deallocate(MBpsi,stat=astat)



end subroutine deallocate_Beltrami_matrices



!> \brief allocate geometry matrices
!>
!> @param vvol
!> @param LcomputeDerivatives
subroutine allocate_geometry_matrices(vvol, LcomputeDerivatives)

! Allocate all geometry dependent matrices for a given ll

  use constants, only: zero

  use fileunits

  use inputlist, only:  Wmemory, Wmacros, Mpol, Lrad

  use allglobal

  use cputiming



  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  REAL      :: cput, cpui, cpuo=0

  INTEGER         :: vvol

  LOGICAL, intent(in) :: LcomputeDerivatives

  INTEGER         :: ll, lldof, jjdof, iidof

  

  ll = Lrad(vvol)

  if (Lcoordinatesingularity) then ! different radial dof for Zernike; 02 Jul 19
    lldof = (Lrad(vvol) - mod(Lrad(vvol),2)) / 2

      ! we need full-size matrices
      iidof = mn
      jjdof = mn

  else
    lldof = Lrad(vvol)

      iidof = mn
      jjdof = mn

  end if

  SALLOCATE( guvijsave, (1:Ntz,1:3,1:3,1:Iquad(vvol)), zero)

  SALLOCATE( DToocc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TTssss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TDstsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( TDszsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDttcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDtzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
  SALLOCATE( DDzzcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

  SALLOCATE( Tss, (0:lldof,1:mn), zero )
  SALLOCATE( Dtc, (0:lldof,1:mn), zero )
  SALLOCATE( Dzc, (0:lldof,1:mn), zero )
  SALLOCATE( Ttc, (0:lldof,1:mn), zero )
  SALLOCATE( Tzc, (0:lldof,1:mn), zero )

  if (NOTstellsym) then

    SALLOCATE( DToocs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DToosc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DTooss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TTsscc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsscs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TTsssc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDstcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDstss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( TDszcc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( TDszss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDttcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDttss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDtzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDtzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( DDzzcs, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzsc, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )
    SALLOCATE( DDzzss, (0:lldof,0:lldof,1:iidof,1:jjdof), zero )

    SALLOCATE( Tsc, (0:lldof,1:mn), zero )
    SALLOCATE( Dts, (0:lldof,1:mn), zero )
    SALLOCATE( Dzs, (0:lldof,1:mn), zero )
    SALLOCATE( Tts, (0:lldof,1:mn), zero )
    SALLOCATE( Tzs, (0:lldof,1:mn), zero )

  end if !NOTstellsym



end subroutine allocate_geometry_matrices



!> \brief deallocate geometry matrices
!>
!> @param LcomputeDerivatives
subroutine deallocate_geometry_matrices(LcomputeDerivatives)

! Deallocate all geometry dependent matrices
  use constants, only: zero

  use fileunits

  use inputlist, only: Wmemory, Wmacros

  use allglobal

  use cputiming

  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  REAL      :: cput, cpui, cpuo=0

  LOGICAL, intent(in) :: LcomputeDerivatives

  

  Lsavedguvij = .false.
  deallocate(guvijsave,stat=astat)

  deallocate(DToocc,stat=astat)
  deallocate(TTssss,stat=astat)
  deallocate(TDstsc,stat=astat)
  deallocate(TDszsc,stat=astat)
  deallocate(DDttcc,stat=astat)
  deallocate(DDtzcc,stat=astat)
  deallocate(DDzzcc,stat=astat)

  deallocate(Tss,stat=astat)
  deallocate(Dtc,stat=astat)
  deallocate(Dzc,stat=astat)
  deallocate(Ttc,stat=astat)
  deallocate(Tzc,stat=astat)

  if (NOTstellsym) then

    deallocate(DToocs,stat=astat)
    deallocate(DToosc,stat=astat)
    deallocate(DTooss,stat=astat)

    deallocate(TTsscc,stat=astat)
    deallocate(TTsscs,stat=astat)
    deallocate(TTsssc,stat=astat)

    deallocate(TDstcc,stat=astat)
    deallocate(TDstcs,stat=astat)
    deallocate(TDstss,stat=astat)

    deallocate(TDszcc,stat=astat)
    deallocate(TDszcs,stat=astat)
    deallocate(TDszss,stat=astat)

    deallocate(DDttcs,stat=astat)
    deallocate(DDttsc,stat=astat)
    deallocate(DDttss,stat=astat)

    deallocate(DDtzcs,stat=astat)
    deallocate(DDtzsc,stat=astat)
    deallocate(DDtzss,stat=astat)

    deallocate(DDzzcs,stat=astat)
    deallocate(DDzzsc,stat=astat)
    deallocate(DDzzss,stat=astat)

    deallocate(Tsc,stat=astat)
    deallocate(Dts,stat=astat)
    deallocate(Dzs,stat=astat)
    deallocate(Tts,stat=astat)
    deallocate(Tzs,stat=astat)

  endif



end subroutine deallocate_geometry_matrices
