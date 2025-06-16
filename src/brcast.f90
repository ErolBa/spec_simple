
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
                        Lhessianallocated, LGdof, dFFdRZ, dBBdmp, dmupfdx, &
                        denergydrr,denergydzr,Lhessian3Dallocated, &
                        lBBintegral, lABintegral, &
                        vvolume, &
                        NOTstellsym, LocalConstraint, &
						IsMyVolume, IsMyVolumeValue



  LOCALS

  INTEGER, intent(in) :: lvol

  INTEGER             :: llmodnp, io, iRZl, ii, ideriv, Nbc

  BEGIN(brcast)






  FATAL( brcast, lvol.le.0 .or. lvol.gt.Mvol, error )



  llmodnp = modulo(lvol-1,ncpu) ! identify which node contains data; this must be consistent with previous looping / parallelization;



  RlBCAST(     mu(lvol), 1, llmodnp )
  RlBCAST( dtflux(lvol), 1, llmodnp )
  RlBCAST( dpflux(lvol), 1, llmodnp )
  RlBCAST( helicity(lvol), 1, llmodnp)

  RlBCAST(     vvolume(lvol), 1, llmodnp )
  RlBCAST( lBBintegral(lvol), 1, llmodnp )
  RlBCAST( lABintegral(lvol), 1, llmodnp )

  RlBCAST( diotadxup(0:1,-1:2,lvol), 8, llmodnp )
  RlBCAST( dItGpdxtp(0:1,-1:2,lvol), 8, llmodnp )



  if( Lhessianallocated ) then


   if( LocalConstraint ) then
 	  Nbc =             LGdof*       2*  LGdof*  2
 	  RlBCAST( dFFdRZ(1:LGdof,0:1,1:LGdof,0:1,lvol), Nbc, llmodnp )
    

	  Nbc =             LGdof*       2*  2
	  RlBCAST( dBBdmp(1:LGdof,lvol,0:1,1:2), Nbc, llmodnp )

	  Nbc =                   2*  LGdof*  2
	  RlBCAST( dmupfdx(lvol,1:1   ,1:2,1:LGdof,0:1), Nbc, llmodnp ) ! why is this broadcast; 02 Sep 14;
   endif


  endif ! end of if( Lhessianallocated ) ; 12 Sep 16;

  if (Lhessian3Dallocated) then

      Nbc =             LGdof*       2*  LGdof*  2
      RlBCAST(denergydrr(1:LGdof,lvol,0:1,1:LGdof,0:1), Nbc, llmodnp )
      RlBCAST(denergydzr(1:LGdof,lvol,0:1,1:LGdof,0:1), Nbc, llmodnp )
  endif
  


  LlBCAST( ImagneticOK(lvol), 1, llmodnp )



  RlBCAST( Bemn(1:mn,lvol,0:1), 2*mn, llmodnp ) ! perhaps all these should be re-ordered; 18 Jul 14;
  RlBCAST( Iomn(1:mn,lvol    ),   mn, llmodnp )
  RlBCAST( Somn(1:mn,lvol,0:1), 2*mn, llmodnp )
  RlBCAST( Pomn(1:mn,lvol,0:2), 3*mn, llmodnp ) ! 15 Sep 15;

  if( NOTstellsym ) then

      RlBCAST( Bomn(1:mn,lvol,0:1), 2*mn, llmodnp )
      RlBCAST( Iemn(1:mn,lvol    ),   mn, llmodnp )
      RlBCAST( Semn(1:mn,lvol,0:1), 2*mn, llmodnp )
      RlBCAST( Pemn(1:mn,lvol,0:2), 3*mn, llmodnp )
  endif ! end of if( NOTstellsym) ; 11 Aug 14;



  if( lvol.gt.Nvol                         .and. Wcurent ) then ! 27 Feb 17;
   RlBCAST( curtor, 1, llmodnp )
   RlBCAST( curpol, 1, llmodnp )
  endif


  RETURN(brcast)



end subroutine brcast


