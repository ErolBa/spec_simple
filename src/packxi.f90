
subroutine packxi( NGdof, position, Mvol, mn, iRbc, iZbs, iRbs, iZbc, packorunpack, LComputeDerivatives, LComputeAxis )



  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wpackxi, Igeometry, Ntor, Nvol, Lfindzero

  use cputiming, only : Tpackxi

  use allglobal, only : ncpu, myid, cpus, im, in, MPI_COMM_SPEC, &
                        YESstellsym, NOTstellsym, &
                        ajk, Nt, Nz, Ntz, iRij, iZij, tRij, tZij, &
                        ijreal, ijimag, jireal, jiimag, efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        psifactor, Rscale



  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  real(8)      :: cput, cpui, cpuo=0

  LOGICAL, intent(in)    :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  LOGICAL, intent(in)    :: LComputeAxis        ! if to recompute the axis

  INTEGER, intent(in)    :: NGdof, Mvol, mn
  real(8)                   :: position(0:NGdof), iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol)
  CHARACTER              :: packorunpack

  INTEGER                :: lvol, jj, kk, irz, issym, idof, ifail, ivol

  



  idof = 0 ! initialize counter; 14 Jan 13;



  do lvol = 1, Mvol-1 ! loop over internal interfaces;

   do jj = 1, mn ! loop over Fourier harmonics;

    do irz = 0, 1 ! loop over R & Z;

     if( Igeometry.lt.3 .and. irz.eq.1 ) cycle ! no dependence on Z; 14 Jan 13;

     do issym = 0, 1 ! loop over even & odd;

      if( YESstellsym .and. issym.eq.1 ) cycle

      if( issym.eq.0 .and. irz.eq.1 .and. jj.eq.1 ) cycle ! no dependence on Zbs_{0,0}; 14 Jan 13;
      if( issym.eq.1 .and. irz.eq.0 .and. jj.eq.1 ) cycle ! no dependence on Rbs_{0,0}; 14 Jan 13;

      idof = idof + 1

      select case( packorunpack )

      case( 'P' ) !   pack vector of unknowns;

       if( irz.eq.0 .and. issym.eq.0 ) position(idof) = iRbc(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) position(idof) = iZbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) position(idof) = iRbs(jj,lvol) / psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) position(idof) = iZbc(jj,lvol) / psifactor(jj,lvol)

      case( 'U' ) ! unpack vector of unknowns;

       if( irz.eq.0 .and. issym.eq.0 ) iRbc(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.0 ) iZbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.0 .and. issym.eq.1 ) iRbs(jj,lvol) = position(idof) * psifactor(jj,lvol)
       if( irz.eq.1 .and. issym.eq.1 ) iZbc(jj,lvol) = position(idof) * psifactor(jj,lvol)

      end select

     enddo ! end of do issym;

    enddo ! end of do irz;

   enddo ! end of do jj;

  enddo ! end of do lvol;



  if( YESstellsym ) then ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1:mn,0:Mvol) = zero
   ;                     ; iZbc(1:mn,0:Mvol) = zero
  else                   ! iRbc(    ,0:Mvol) = zero
   ;                     ; iZbs(1   ,0:Mvol) = zero
   ;                     ; iRbs(1   ,0:Mvol) = zero
   ;                     ! iZbc(    ,0:Mvol) = zero
  endif


  select case( packorunpack )

  case( 'P' )

  case( 'U' )

   ivol = 1 ! take care with ivol: this variable name might be a global variable, but here it is local; 19 Jul 16;


  end select






end subroutine packxi


