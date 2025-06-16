





















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



subroutine intghs( lquad, mn, lvol, lrad, idx )



  use constants, only : zero, half, one, two, pi, pi2

  use numerical, only : vsmall, small, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : mpol, Wintghs, Wmacros

  use cputiming, only : Tintghs

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, &
                        Mvol, im, in, mne, Ntz, Nt, Nz, &
                        YESstellsym, NOTstellsym, &
                        gaussianweight, gaussianabscissae, &
                        Tsc, Tss, Dtc, Dts, Dzc, Dzs, &
                        Ttc, Tts, Tzc, Tzs, &
                        Lcoordinatesingularity, &
                        pi2pi2nfp, pi2pi2nfphalf, dBdX, &
                        Ntz, NOTstellsym, dBdX, Lsavedguvij, &
                        Ate, Ato, Aze, Azo, &
                        sg, guvij, guvijsave

  use intghs_module


  LOCALS

  INTEGER, intent(in) :: lquad, mn, lvol, lrad, idx

  INTEGER             :: jquad, ll, pp, ll1, pp1, uv, ii, jj, io, mn2, lp2, mn2_max, lp2_max, nele, ideriv, ifail, Lcurvature
  INTEGER             :: mi, ni, bid

  REAL                :: lss, jthweight, Tl, Dl, sbar, dfactor, ik, w(1:lquad)

  BEGIN( intghs )



  Tss = zero
  Dtc = zero
  Dzc = zero

  if (.not.dBdX%L) then
    Ttc = zero
    Tzc = zero
  endif

  if (NOTstellsym) then
    Tsc = zero
    Dts = zero
    Dzs = zero
    if (.not.dBdX%L) then
      Tts = zero
      Tzs = zero
    endif

  endif !NOTstellsym

  if( dBdX%L ) then ; Lcurvature = 3 ; ideriv = 1
  else              ; Lcurvature = 1 ; ideriv = 0
  endif

  if (.not. Lsavedguvij) then
    WCALL( intghs, compute_guvijsave, (lquad, lvol, ideriv, Lcurvature) )
  endif

  do jquad = 1, lquad
    lss = gaussianabscissae(jquad,lvol) ; jthweight = gaussianweight(jquad,lvol)
    sbar = (lss + one) * half
    if (Lcoordinatesingularity) then
      call get_zernike(sbar, lrad, mpol, wk%basis(:,:,0:1,jquad)) ! use Zernike polynomials 29 Jun 19;
    else
      call get_cheby(lss, lrad, wk%basis(:,0,0:1,jquad))
    endif
  enddo

  wk%efmn = zero ; wk%ofmn = zero ; wk%cfmn = zero ; wk%sfmn = zero
  wk%evmn = zero ; wk%odmn = zero ;

  wk%ijreal = zero ; wk%jireal = zero ; wk%jkreal = zero ; wk%kjreal = zero

  wk%gBupper = zero; wk%Blower = zero;


  do jquad = 1, lquad


    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;

      if (Lcoordinatesingularity) then
        do ll = mi, lrad, 2 ! loop over Zernike polynomials; Lrad is the radial resolution;
        ;                      ; wk%efmn(ii,jquad) = wk%efmn(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
        ;                      ; wk%cfmn(ii,jquad) = wk%cfmn(ii,jquad) - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
        ;                      ; wk%odmn(ii,jquad) = wk%odmn(ii,jquad) - Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * ni &
                                                    - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * mi
        ;                      ; wk%ijreal(ii,jquad) = wk%ijreal(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        ;                      ; wk%jkreal(ii,jquad) = wk%jkreal(ii,jquad) + Aze(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        if( NOTstellsym ) then ; wk%ofmn(ii,jquad) = wk%ofmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
          ;                    ; wk%sfmn(ii,jquad) = wk%sfmn(ii,jquad) - Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,1,jquad)*half)
          ;                    ; wk%evmn(ii,jquad) = wk%evmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * ni &
                                                    + Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,mi,0,jquad)) * mi
          ;                    ; wk%jireal(ii,jquad) = wk%jireal(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
          ;                    ; wk%kjreal(ii,jquad) = wk%kjreal(ii,jquad) + Azo(lvol,idx,ii)%s(ll) * wk%basis(ll,mi,0,jquad)
        endif
        enddo ! end of do ll; 20 Feb 13;
      else
        do ll = 0, lrad ! loop over Chebyshev polynomials; Lrad is the radial resolution;
        ;                      ; wk%efmn(ii,jquad) = wk%efmn(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
        ;                      ; wk%cfmn(ii,jquad) = wk%cfmn(ii,jquad) - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
        ;                      ; wk%odmn(ii,jquad) = wk%odmn(ii,jquad) - Ate(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * ni &
                                                    - Aze(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * mi
        ;                      ; wk%ijreal(ii,jquad) = wk%ijreal(ii,jquad) + Ate(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        ;                      ; wk%jkreal(ii,jquad) = wk%jkreal(ii,jquad) + Aze(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        if( NOTstellsym ) then ; wk%ofmn(ii,jquad) = wk%ofmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
          ;                    ; wk%sfmn(ii,jquad) = wk%sfmn(ii,jquad) - Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,1,jquad))
          ;                    ; wk%evmn(ii,jquad) = wk%evmn(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * ni &
                                                    + Azo(lvol,idx,ii)%s(ll) * ( wk%basis(ll,0,0,jquad)) * mi
          ;                    ; wk%jireal(ii,jquad) = wk%jireal(ii,jquad) + Ato(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
          ;                    ; wk%kjreal(ii,jquad) = wk%kjreal(ii,jquad) + Azo(lvol,idx,ii)%s(ll) * wk%basis(ll,0,0,jquad)
        endif
        enddo ! end of do ll; 20 Feb 13;
      end if ! Lcoordinatesingularity; 01 Jul 19
    enddo ! end of do ii; 20 Feb 13;

    call invfft( mn, im, in, wk%efmn(1:mn,jquad), wk%ofmn(1:mn,jquad), wk%cfmn(1:mn,jquad), wk%sfmn(1:mn,jquad), Nt, Nz, wk%gBupper(1:Ntz,3,jquad), wk%gBupper(1:Ntz,2,jquad) )
    call invfft( mn, im, in, wk%evmn(1:mn,jquad), wk%odmn(1:mn,jquad), wk%cfmn(1:mn,jquad), wk%sfmn(1:mn,jquad), Nt, Nz, wk%gBupper(1:Ntz,1,jquad), wk%gBupper(1:Ntz,2,jquad) )

    do ii = 1, 3
      do jj = 1, 3
        wk%Blower(:,ii,jquad) = wk%Blower(:,ii,jquad) + wk%gBupper(:,jj,jquad) * guvijsave(1:Ntz,jj,ii,jquad)
      enddo
    enddo



    ifail = 0

    call tfft( Nt, Nz, wk%Blower(1:Ntz,2,jquad), wk%Blower(1:Ntz,3,jquad), &
               mn, im(1:mn), in(1:mn), wk%Bloweremn(1:mn,2,jquad), wk%Bloweromn(1:mn,2,jquad), wk%Bloweremn(1:mn,3,jquad), wk%Bloweromn(1:mn,3,jquad), ifail )
    call tfft( Nt, Nz, wk%Blower(1:Ntz,1,jquad), wk%Blower(1:Ntz,3,jquad), &
               mn, im(1:mn), in(1:mn), wk%Bloweromn(1:mn,1,jquad), wk%Bloweremn(1:mn,1,jquad), wk%Bloweremn(1:mn,3,jquad), wk%Bloweromn(1:mn,3,jquad), ifail )
    wk%Bloweremn(1,1  ,jquad) = zero
    wk%Bloweromn(1,2:3,jquad) = zero

  end do

w(1:lquad) = gaussianweight(1:lquad,lvol)

do ii = 1, mn

  if (ii==1) then ;ik = two
  else            ;ik = one
  endif

    do ll = 0, lrad

      if (Lcoordinatesingularity) then
        if (ll.lt.im(ii) .or. mod(ll+im(ii),2).ne.0) cycle
        ll1 = (ll - mod(ll,2))/2 ! shrinked dof for Zernike; 02 Jul 19
        bid = im(ii)
        dfactor = half
      else
        ll1 = ll
        bid = 0
        dfactor = one
      endif

      Tss(ll1,ii) = sum(wk%basis(ll, bid, 0, :) * wk%Bloweremn(ii,1,:) * w) * ik
      Dtc(ll1,ii) = sum(wk%basis(ll, bid, 1, :) * wk%Bloweremn(ii,2,:) * w) * dfactor * ik
      Dzc(ll1,ii) = sum(wk%basis(ll, bid, 1, :) * wk%Bloweremn(ii,3,:) * w) * dfactor * ik

      if (NOTstellsym) then
        Tsc(ll1,ii) =  sum(wk%basis(ll, bid, 0, :) * wk%Bloweromn(ii,1,:) * w) * ik
        Dts(ll1,ii) =  sum(wk%basis(ll, bid, 1, :) * wk%Bloweromn(ii,2,:) * w) * dfactor * ik
        Dzs(ll1,ii) =  sum(wk%basis(ll, bid, 1, :) * wk%Bloweromn(ii,3,:) * w) * dfactor * ik
      endif

      if (dBdX%L) cycle ! dMD matrix does not depend on geometry
      Ttc(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%cfmn(ii,:) * w) + sum(wk%basis(ll, bid, 1, :) * wk%jkreal(ii,:) * w) * dfactor) * ik
      Tzc(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%efmn(ii,:) * w) - sum(wk%basis(ll, bid, 1, :) * wk%ijreal(ii,:) * w) * dfactor) * ik

      if (NOTstellsym) then
        Tts(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%sfmn(ii,:) * w) + sum(wk%basis(ll, bid, 1, :) * wk%kjreal(ii,:) * w) * dfactor) * ik
        Tzs(ll1,ii) = (sum(wk%basis(ll, bid, 0, :) * wk%ofmn(ii,:) * w) - sum(wk%basis(ll, bid, 1, :) * wk%jireal(ii,:) * w) * dfactor) * ik
      endif

    enddo !ll
  enddo !ii



  Tss = Tss * pi2pi2nfphalf
  Dtc = Dtc * pi2pi2nfphalf
  Dzc = Dzc * pi2pi2nfphalf

  if (.not.dBdX%L) then
    Ttc = Ttc * pi2pi2nfphalf
    Tzc = Tzc * pi2pi2nfphalf
  endif

  if (NOTstellsym) then

    Tsc = Tsc * pi2pi2nfphalf
    Dts = Dts * pi2pi2nfphalf
    Dzs = Dzs * pi2pi2nfphalf

    if (.not.dBdX%L) then
      Tts = Tts * pi2pi2nfphalf
      Tzs = Tzs * pi2pi2nfphalf
    endif

  endif



  RETURN( intghs )



end subroutine intghs



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
