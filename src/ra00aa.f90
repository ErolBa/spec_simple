
subroutine ra00aa( writeorread )



  use constants, only : zero

  use numerical, only :

  use fileunits, only : ounit, aunit

  use inputlist, only : Wmacros, Wra00aa, Nfp, Mpol, Ntor, Lrad

  use cputiming, only : Tra00aa

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, ext, Mvol, mn, im, in, Ate, Aze, Ato, Azo

  use sphdf5,    only : write_vector_potential



  LOCALS

  CHARACTER, intent(in) :: writeorread

  LOGICAL               :: exist

  INTEGER               :: vvol, oldMvol, oldMpol, oldNtor, oldmn, oldNfp, oldLrad, ii, jj, minLrad, llmodnp, ideriv, sumLrad
  INTEGER, allocatable  :: oldim(:), oldin(:)
  REAL   , allocatable  :: oldAte(:), oldAze(:), oldAto(:), oldAzo(:)
  REAL   , allocatable  :: allAte(:,:), allAze(:,:), allAto(:,:), allAzo(:,:)


  



  ideriv = 0 ! write vector potential and not derivative (?)



  select case( writeorread )



  case( 'W' ) ! write vector potential harmonics to file;

   sumLrad = sum(Lrad(1:Mvol)+1)

   allocate(allAte(1:sumLrad,1:mn))
   allocate(allAze(1:sumLrad,1:mn))
   allocate(allAto(1:sumLrad,1:mn))
   allocate(allAzo(1:sumLrad,1:mn))

   sumLrad = 1
   do vvol = 1, Mvol


    do ii = 1, mn ! loop over Fourier harmonics;
     allAte(sumLrad:sumLrad+Lrad(vvol),ii) = Ate(vvol,ideriv,ii)%s(0:Lrad(vvol))
     allAze(sumLrad:sumLrad+Lrad(vvol),ii) = Aze(vvol,ideriv,ii)%s(0:Lrad(vvol))
     allAto(sumLrad:sumLrad+Lrad(vvol),ii) = Ato(vvol,ideriv,ii)%s(0:Lrad(vvol))
     allAzo(sumLrad:sumLrad+Lrad(vvol),ii) = Azo(vvol,ideriv,ii)%s(0:Lrad(vvol))
    enddo ! end of do ii;  6 Feb 13;
    sumLrad = sumLrad+Lrad(vvol)+1
   enddo ! end of do vvol;  6 Feb 13;

   sumLrad = sum(Lrad(1:Mvol)+1)

   WCALL( ra00aa, write_vector_potential, (sumLrad, allAte, allAze, allAto, allAzo) )

   deallocate(allAte)
   deallocate(allAze)
   deallocate(allAto)
   deallocate(allAzo)



  case( 'R' ) ! read potential from file; interpolate onto new radial grid;


   if( myid.eq.0 ) then

    inquire(file="."//trim(ext)//".sp.A",exist=exist)

    if( .not.exist ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; .ext.sp.A does not exist ;")') cput-cpus, myid ; goto 9998
    endif

    open(aunit,file="."//trim(ext)//".sp.A",status="old",form="unformatted",iostat=ios) ! this will contain initial guess for vector potential;

    if( ios.ne.0 ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; opening .ext.sp.A ;")') cput-cpus, myid ; goto 9997
    endif

    read(aunit,iostat=ios) oldMvol, oldMpol, oldNtor, oldmn, oldNfp ! these are the "old" resolution parameters;

    if( ios.ne.0 ) then
     write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; reading oldMvol, oldMpol, oldNtor, oldmn, oldNfp;")') cput-cpus, myid
     goto 9997
    endif

    if( oldNfp .ne.Nfp  ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Nfp ; ")') cput-cpus, myid ; goto 9997
    endif
    if( oldMvol.ne.Mvol ) then ; write(ounit,'("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Mvol ;")') cput-cpus, myid ; goto 9997
    endif

    SALLOCATE( oldim, (1:oldmn), 0 )
    SALLOCATE( oldin, (1:oldmn), 0 )

    read(aunit,iostat=ios) oldim(1:oldmn)
    read(aunit,iostat=ios) oldin(1:oldmn)

    do vvol = 1, oldMvol

     read(aunit,iostat=ios) oldLrad

     minLrad = min(oldLrad,Lrad(vvol))

     SALLOCATE( oldAte, (0:oldLrad), zero )
     SALLOCATE( oldAze, (0:oldLrad), zero )
     SALLOCATE( oldAto, (0:oldLrad), zero )
     SALLOCATE( oldAzo, (0:oldLrad), zero )

     do jj = 1, oldmn

      read(aunit,iostat=ios) oldAte(0:oldLrad)
      read(aunit,iostat=ios) oldAze(0:oldLrad)
      read(aunit,iostat=ios) oldAto(0:oldLrad)
      read(aunit,iostat=ios) oldAzo(0:oldLrad)

      do ii = 1, mn ! compare Fourier harmonic with old; 26 Feb 13;
       if( im(ii).eq.oldim(jj) .and. in(ii).eq.oldin(jj) ) then ; Ate(vvol,ideriv,ii)%s(0:minLrad) = oldAte(0:minLrad)
        ;                                                       ; Aze(vvol,ideriv,ii)%s(0:minLrad) = oldAze(0:minLrad)
        ;                                                       ; Ato(vvol,ideriv,ii)%s(0:minLrad) = oldAto(0:minLrad)
        ;                                                       ; Azo(vvol,ideriv,ii)%s(0:minLrad) = oldAzo(0:minLrad)
       endif
      enddo ! end of do ii; 26 Feb 13;

     enddo ! end of do jj; 26 Feb 13;

     DALLOCATE(oldAte)
     DALLOCATE(oldAze)
     DALLOCATE(oldAto)
     DALLOCATE(oldAzo)

    enddo ! end of do vvol; 26 Feb 13;

    DALLOCATE(oldim)
    DALLOCATE(oldin)

9997 continue

    close(aunit)

9998 continue

   endif  ! end of if( myid.eq.0 ) ; 26 Feb 13;

   do vvol = 1, Mvol

    llmodnp = 0 ! this node contains the information that is to be broadcast; 26 Feb 13;

    do ii = 1, mn
     RlBCAST( Ate(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
     RlBCAST( Aze(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
    enddo
    do ii = 1, mn
     RlBCAST( Ato(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
     RlBCAST( Azo(vvol,ideriv,ii)%s(0:Lrad(vvol)), Lrad(vvol)+1, llmodnp )
    enddo

   enddo ! end of do vvol; 26 Feb 13;



  case default

   FATAL(ra00aa, .true., invalid writeorread flag supplied on input )

  end select







end subroutine ra00aa




