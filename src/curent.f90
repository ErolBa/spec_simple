
subroutine curent( lvol, mn, Nt, Nz, iflag, ldItGp )



  use constants, only : zero, one, two, pi2

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wcurent, Lrad

  use cputiming, only : Tcurent

  use allglobal, only : ncpu, cpus, myid, MPI_COMM_SPEC, &
                        Mvol, im, in, mne, ime, ine, &
                        YESstellsym, NOTstellsym, &
                        sg, guvij, &
                        Ntz, ijreal, ijimag, jireal, jiimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                        Ate, Aze, Ato, Azo, TT, &
                        build_vector_potential



  LOCALS

  INTEGER, intent(in)  :: lvol, mn, Nt, Nz, iflag
  REAL   , intent(out) :: ldItGp(0:1,-1:2)

  INTEGER              :: innout, ideriv, ii, ll, Lcurvature, ifail
  REAL                 :: lss
  REAL                 :: Bsupt(1:Nt*Nz,-1:2), Bsupz(1:Nt*Nz,-1:2)

  










  innout = 0 ; lss = two * innout - one ! lvol = Mvol ! internal shorthand/variables; 08 Feb 16; note that lvol = Mvol is hardwired; 12 Sep 16;



  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;

   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;

   call build_vector_potential(lvol, innout, ideriv, 1)

   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), &
                Nt, Nz, Bsupz(1:Ntz,ideriv), Bsupt(1:Ntz,ideriv) ) ! map to real space;

  enddo ! end of do ideriv; 31 Jan 13;



  select case( iflag )
  case( -1 ) ; Lcurvature = 3 ! will need derivatives of metrics w.r.t. interface deformation; 15 Sep 16;
  case(  1 ) ; Lcurvature = 1
  case(  2 ) ; Lcurvature = 1
  end select

  WCALL( curent, coords,( lvol, lss, Lcurvature, Ntz, mn ) ) ! get "lower" metric elements evaluated on innout interface;



  do ideriv = -1, 2 ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;

   if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of currents                                                   are not required; 20 Jun 14;
   if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of currents  wrt geometry                                     is  not required; 20 Jun 14;
   if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of currents  wrt enclosed toroidal and enclosed poloidal flux are not required; 20 Jun 14;

   ijreal(1:Ntz) =                 ( - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,2,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) =                 ( - Bsupt(1:Ntz,ideriv) * guvij(1:Ntz,2,3,0) + Bsupz(1:Ntz,ideriv) * guvij(1:Ntz,3,3,0) ) / sg(1:Ntz,0)
   if( ideriv.eq.-1 ) then ! add derivatives of metrics with respect to interface geometry; 15 Sep 16;
   ijreal(1:Ntz) = ijreal(1:Ntz) + ( - Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,2,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,2,3,1) ) / sg(1:Ntz,0)
   ijimag(1:Ntz) = ijimag(1:Ntz) + ( - Bsupt(1:Ntz,     0) * guvij(1:Ntz,2,3,1) + Bsupz(1:Ntz,     0) * guvij(1:Ntz,3,3,1) ) / sg(1:Ntz,0)
   endif

   ifail = 0
   call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
              mne, ime(1:mne), ine(1:mne), efmn(1:mne), ofmn(1:mne), cfmn(1:mne), sfmn(1:mne), ifail )

   ldItGp(0,ideriv) = efmn(1) * pi2 ! "toroidal", plasma  current; 12 Sep 16;
   ldItGp(1,ideriv) = cfmn(1) * pi2 ! "poloidal", linking current; 12 Sep 16;



  enddo ! end of do ideriv; 31 Jan 13;






end subroutine curent


