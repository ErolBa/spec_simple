










subroutine bnorml( mn, Ntz, efmn, ofmn )



  use constants, only : zero, half, one, two, pi, pi2, ten

  use numerical, only : small

  use fileunits, only : ounit, lunit

  use inputlist, only : Wmacros, Wbnorml, Igeometry, Lcheck, vcasingtol, vcasingper, Lrad

  use cputiming, only : Tbnorml

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, pi2nfp, Mvol, &
                        Nt, Nz, &
                        Rij, Zij, guvij, sg, TT, &
                        NOTstellsym, Lcoordinatesingularity, &
                        im, in, Ate, Aze, Ato, Azo, &
                        Nt, Nz, cfmn, sfmn, &
                        ijreal, ijimag, jireal, jiimag, &
                        globaljk, tetazeta, virtualcasingfactor, gteta, gzeta, Dxyz, Nxyz



  LOCALS

  INTEGER, intent(in)  :: mn, Ntz
  REAL   , intent(out) :: efmn(1:mn), ofmn(1:mn)

  INTEGER              :: lvol, Lcurvature, Lparallel, ii, jj, kk, jk, ll, kkmodnp, jkmodnp, ifail, id01daf, nvccalls, icasing, ideriv
  REAL                 :: lss, zeta, teta, cszeta(0:1), tetalow, tetaupp, absacc, gBn
  REAL                 :: Jxyz(1:Ntz,1:3), Bxyz(1:Ntz,1:3), dAt(1:Ntz), dAz(1:Ntz), distance(1:Ntz)


  BEGIN(bnorml)



  Lparallel = 1 ! controls choice of parallelization; see below;



  ijreal(1:Ntz) = zero ! normal plasma field; 15 Oct 12;




  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz

   if( Igeometry.eq.3 ) then ; cszeta(0:1) = (/ cos(zeta), sin(zeta) /)
   endif

   do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt

    globaljk = jk ! this is global; passed through to vcintegrand & casing;

    select case( Lparallel ) ! perform in parallel;
    case( 0 ) ! Lparallel = 0 ; 09 Mar 17;
     if( myid.ne.modulo(kk,ncpu) ) cycle
    case( 1 ) ! Lparallel = 1 ; 09 Mar 17;
     if( myid.ne.modulo(jk-1,ncpu) ) cycle ! 11 Oct 12; this is a weird parallelization, but perhaps better exploits all available cpus;
    case default ! Lparallel; 09 Mar 17;
     FATAL( bnorml, .true., invalid Lparallel in parallelization loop )
    end select ! end of select case( Lparallel ) ; 09 Mar 17;

    tetazeta(1:2) = (/ teta, zeta /) ! this is global; passed through to zetalow & zetaupp; 14 Apr 17;


    WCALL( bnorml, casing, ( teta, zeta, gBn, icasing ) ) ! tetazeta is global; 26 Apr 17;

    ijreal(jk) = gBn


1000 format("bnorml : ", 10x ," : myid=",i3," : \z =",f6.3," ; \t =",f6.3," ; B . x_t x x_z =",2f22.15," ; ":"err =",es13.5," ;")

   enddo ! end of do jj;
  enddo ! end of do kk;



1001 format("bnorml : ", 10x ," : "a1" : (t,z) = ("f8.4","f8.4" ) ; gBn=",f23.15," ; ":" error =",f23.15" ;")



  if( myid.eq.0 .and. Lcheck.eq.6 ) then ! THIS WAS CORRUPTED; see before 14 Apr 17 for complete source;
   close(lunit)
  endif



  do kk = 0, Nz-1

   kkmodnp = modulo(kk,ncpu)

   select case( Lparallel )

   case( 0 ) ! Lparallel = 0 ; 09 Mar 17;

    RlBCAST(ijreal(1+kk*Nt:Nt+kk*Nt),Nt,kkmodnp) ! plasma; 03 Apr 13;


   case( 1 ) ! Lparallel = 1 ; 09 Mar 17;

    do jj = 0, Nt-1

     jk = 1 + jj + kk*Nt

     jkmodnp = modulo(jk-1,ncpu)

     RlBCAST(ijreal(jk),1,jkmodnp) ! plasma; 03 Apr 13;


    enddo

   case default ! Lparallel; 09 Mar 17;

    FATAL( bnorml, .true., invalid Lparallel for broadcasting )

   end select ! end of select case( Lparallel ) ; 09 Mar 17;

  enddo ! 11 Oct 12;



  ijreal(1:Ntz) = ijreal(1:Ntz) * virtualcasingfactor
  ijimag(1:Ntz) = zero

  call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
             mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier decompose normal field;



  RETURN(bnorml)



end subroutine bnorml
