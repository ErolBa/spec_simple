
subroutine brcast( lvol )



  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wbrcast, Wcurent, MNvol, Nvol, mu, Lrad, &
                        curtor, curpol, Lconstraint, Lfindzero, helicity

  use cputiming, only : Tbrcast

  use allglobal, only : myid, cpus, ncpu, MPI_COMM_SPEC, &
                        dtflux, dpflux, Ntz, mn, Mvol, &
                        diotadxup, dItGpdxtp, &
                        Ate, Aze, Ato, Azo, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, Pomn, Pemn, &
                        ImagneticOK, &
                        LGdof, dFFdRZ, dBBdmp, dmupfdx, &
                        lBBintegral, lABintegral, &
                        vvolume, &
                        NOTstellsym, LocalConstraint, &
						IsMyVolume, IsMyVolumeValue



  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  real(8)      :: cput, cpui, cpuo=0

  INTEGER, intent(in) :: lvol

  INTEGER             :: llmodnp, io, iRZl, ii, ideriv, Nbc

  






if( lvol.le.0 .or. lvol.gt.Mvol ) then
     write(6,'("brcast :      fatal : myid=",i3," ; lvol.le.0 .or. lvol.gt.Mvol ; error;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "brcast : lvol.le.0 .or. lvol.gt.Mvol : error ;"
   endif



  llmodnp = modulo(lvol-1,ncpu) ! identify which node contains data; this must be consistent with previous looping / parallelization;



call MPI_BCAST(mu(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(dtflux(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(dpflux(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(helicity(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

call MPI_BCAST(vvolume(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(lBBintegral(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(lABintegral(lvol), 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

call MPI_BCAST(diotadxup(0:1,-1:2,lvol), 8, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(dItGpdxtp(0:1,-1:2,lvol), 8, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

  call MPI_BCAST(ImagneticOK(lvol),1,MPI_LOGICAL,llmodnp,MPI_COMM_SPEC,ierr)



call MPI_BCAST(Bemn(1:mn,lvol,0:1), 2*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Iomn(1:mn,lvol    ), mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Somn(1:mn,lvol,0:1), 2*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Pomn(1:mn,lvol,0:2), 3*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)

  if( NOTstellsym ) then

call MPI_BCAST(Bomn(1:mn,lvol,0:1), 2*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Iemn(1:mn,lvol    ), mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Semn(1:mn,lvol,0:1), 2*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(Pemn(1:mn,lvol,0:2), 3*mn, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
  endif ! end of if( NOTstellsym) ; 11 Aug 14;



  if( lvol.gt.Nvol                         .and. Wcurent ) then ! 27 Feb 17;
call MPI_BCAST(curtor, 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
call MPI_BCAST(curpol, 1, MPI_DOUBLE_PRECISION, llmodnp, MPI_COMM_SPEC, ierr)
  endif






end subroutine brcast


