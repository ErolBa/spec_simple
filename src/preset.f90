
subroutine preset



  use constants, only : zero, one, mu0

  use numerical, only : sqrtmachprec, vsmall, small

  use fileunits, only : ounit

  use inputlist

  use cputiming, only : Tpreset

  use allglobal

  use fftw_interface



  LOCALS

  INTEGER   :: innout, idof, jk, ll, ii, ifail, ideriv, vvol, mi, ni, mj, nj, mk, nk, mimj, ninj, mkmj, nknj, jj, kk, lvol, mm, nn, imn
  INTEGER   :: lquad, igauleg, maxIquad, Mrad, jquad, Lcurvature, zerdof, iret, work1, work2
  REAL      :: teta, zeta, arg, lss, cszeta(0:1), error
  LOGICAL   :: LComputeAxis

  LOGICAL              :: Lchangeangle
  INTEGER              :: nb, ix, ij, ip, idx_mode
  REAL                 :: xx




  





  select case( Istellsym )
  case( 0 )    ; YESstellsym = .false. ; NOTstellsym = .true.
  case( 1 )    ; YESstellsym = .true.  ; NOTstellsym = .false.
  case default ;
   FATAL( readin, .true., illegal Istellsym )
  end select

  Mvol = Nvol + Lfreebound



  SALLOCATE( beltramierror,(1:Mvol,1:9), zero)




  mn = 1 + Ntor +  Mpol * ( 2 *  Ntor + 1 ) ! Fourier resolution of interface geometry & vector potential;

  SALLOCATE( im, (1:mn), 0 )
  SALLOCATE( in, (1:mn), 0 )

  call gi00ab(  Mpol,  Ntor, Nfp, mn, im(1:mn), in(1:mn) ) ! this sets the im and in mode identification arrays;




  SALLOCATE( halfmm, (1:mn), im(1:mn) * half )
  SALLOCATE( regumm, (1:mn), im(1:mn) * half )

  if( Mregular.ge.2 ) then

   where( im.gt.Mregular ) regumm = Mregular * half

  endif





  lMpol = 4*Mpol ; lNtor = 4*Ntor ! extra-enhanced resolution for metrics;

  mne = 1 + lNtor + lMpol * ( 2 * lNtor + 1 ) ! resolution of metrics; enhanced resolution; see metrix;

  SALLOCATE( ime, (1:mne), 0 )
  SALLOCATE( ine, (1:mne), 0 )

  call gi00ab( lMpol, lNtor, Nfp, mne, ime(1:mne), ine(1:mne) )




  sMpol = iMpol ; sNtor = iNtor

  if( iMpol.le.0 ) sMpol = Mpol - iMpol
  if( iNtor.le.0 ) sNtor = Ntor - iNtor
  if(  Ntor.eq.0 ) sNtor = 0

  mns = 1 + sNtor + sMpol * ( 2 * sNtor + 1 ) ! resolution of straight-field line transformation on interfaces; see tr00ab; soon to be redundant;

  SALLOCATE( ims, (1:mns), 0 )
  SALLOCATE( ins, (1:mns), 0 )

  call gi00ab( sMpol, sNtor, Nfp, mns, ims(1:mns), ins(1:mns) ) ! note that the field periodicity factor is included in ins;













  if( Lcheck.eq.5 ) then ; forcetol = 1.0e+12 ; nPpts = 0 ! will check Hessian using finite-differences;
  endif





  SALLOCATE( iRbc, (1:mn,0:Mvol), zero ) ! interface Fourier harmonics;
  SALLOCATE( iZbs, (1:mn,0:Mvol), zero )
  SALLOCATE( iRbs, (1:mn,0:Mvol), zero )
  SALLOCATE( iZbc, (1:mn,0:Mvol), zero )

  if( Lperturbed.eq.1 ) then
  SALLOCATE( dRbc, (1:mn,0:Mvol), zero ) ! interface Fourier harmonics;
  SALLOCATE( dZbs, (1:mn,0:Mvol), zero )
  SALLOCATE( dRbs, (1:mn,0:Mvol), zero )
  SALLOCATE( dZbc, (1:mn,0:Mvol), zero )
  endif



  SALLOCATE( iVns, (1:mn), zero )
  SALLOCATE( iBns, (1:mn), zero )
  SALLOCATE( iVnc, (1:mn), zero )
  SALLOCATE( iBnc, (1:mn), zero )






  SALLOCATE( ajk, (1:mn), zero ) ! this must be allocated & assigned now, as it is used in readin; primarily used in packxi; 02 Jan 15;

  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)

   if( mk.eq.0 ) ajk(kk) = pi2

  enddo ! end of do kk;



  if( myid.eq.0 ) then ! read plasma boundary & computational boundary; initialize interface geometry;

   if( Igeometry.eq.3 .and. Rbc(0,+1)+Rbc(0,-1).gt.zero .and. Zbs(0,+1)-Zbs(0,-1).gt.zero ) then ; Lchangeangle = .true.
   else                                                                                          ; Lchangeangle = .false.
   endif

   if( Lchangeangle ) write(ounit,'("readin : " 10x " : CHANGING ANGLE ;")')

   do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp ! set plasma boundary, computational boundary; 29 Apr 15;

    if( Lchangeangle ) then ; jj = -1 ; kk = -nn ! change sign of poloidal angle;
    else                    ; jj = +1 ; kk = +nn
    endif

    if( mm.eq.0 .and. nn.eq.0 ) then

     ;iRbc(ii,Nvol) = Rbc( nn, mm)                         ! plasma        boundary is ALWAYS given by namelist Rbc & Zbs;
     ;iZbs(ii,Nvol) = zero
      if( NOTstellsym ) then
     ;iRbs(ii,Nvol) = zero
     ;iZbc(ii,Nvol) = Zbc( nn, mm)
      else
     ;iRbs(ii,Nvol) = zero
     ;iZbc(ii,Nvol) = zero
      endif

    else ! if( mm.eq.0 .and. nn.eq.0 ) then ; matches

     ;iRbc(ii,Nvol) =   Rbc( kk, mm) + Rbc(-kk,-mm)        ! plasma        boundary is ALWAYS given by namelist Rbc & Zbs;
     ;iZbs(ii,Nvol) = ( Zbs( kk, mm) - Zbs(-kk,-mm) ) * jj
      if( NOTstellsym ) then
     ;iRbs(ii,Nvol) = ( Rbs( kk, mm) - Rbs(-kk,-mm) ) * jj
     ;iZbc(ii,Nvol) =   Zbc( kk, mm) + Zbc(-kk,-mm)
      else
     ;iRbs(ii,Nvol) =   zero
     ;iZbc(ii,Nvol) =   zero
      endif

    endif ! end of if( mm.eq.0 .and. nn.eq.0 ) ;

   enddo ! end of do ii = 1, mn;


   select case( Linitialize ) ! 24 Oct 12;

   case( :0 ) ! Linitialize=0 ; initial guess for geometry of the interior surfaces is given in the input file;

    if( Lchangeangle ) then ; jj = -1  ! change sign of poloidal angle; Loizu Nov 18;
    else                    ; jj = +1
    endif

    do idx_mode=1, num_modes! will read in Fourier harmonics until the end of file is reached;
     mm = mmRZRZ(idx_mode)
     nn = nnRZRZ(idx_mode)

     do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over harmonics within range;
      if( mm.eq.0 .and. mi.eq.0 .and. nn*Nfp.eq.ni ) then
       iRbc(ii,1:Nvol-1) = allRZRZ(1,1:Nvol-1, idx_mode) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = allRZRZ(2,1:Nvol-1, idx_mode) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = allRZRZ(3,1:Nvol-1, idx_mode) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = allRZRZ(4,1:Nvol-1, idx_mode) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      elseif( mm.eq.mi .and. nn*Nfp.eq.jj*ni ) then
       iRbc(ii,1:Nvol-1) = allRZRZ(1,1:Nvol-1, idx_mode) ! select relevant harmonics;
       iZbs(ii,1:Nvol-1) = jj*allRZRZ(2,1:Nvol-1, idx_mode) ! select relevant harmonics;
       if( NOTstellsym ) then
        iRbs(ii,1:Nvol-1) = jj*allRZRZ(3,1:Nvol-1, idx_mode) ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = allRZRZ(4,1:Nvol-1, idx_mode) ! select relevant harmonics;
       else
        iRbs(ii,1:Nvol-1) = zero             ! select relevant harmonics;
        iZbc(ii,1:Nvol-1) = zero             ! select relevant harmonics;
       endif
      endif
     enddo ! end of do ii;

    enddo ! end of do;

   end select ! end select case( Linitialize );

   if( Igeometry.eq.3 ) then
    if( Rac(0).gt.zero ) then ! user has supplied logically possible coordinate axis;
     iRbc(1:Ntor+1,0) = Rac(0:Ntor)
     iZbs(1:Ntor+1,0) = Zas(0:Ntor)
     iRbs(1:Ntor+1,0) = Ras(0:Ntor)
     iZbc(1:Ntor+1,0) = Zac(0:Ntor)
    else ! see preset for poloidal-average specification of coordinate axis and geometrical initialization;
    endif ! end of if( Igeometry.eq.3 ) then ;
   endif

  endif ! end of if myid.eq.0 loop; only the master will read the input file; all variables need to be broadcast;



  ; RlBCAST( iRbc(1:mn,0:Mvol), (Mvol+1)*mn, 0 )
  if( Igeometry.eq.3 ) then
   ;RlBCAST( iZbs(1:mn,0:Mvol), (Mvol+1)*mn, 0 ) ! only required for ii > 1 ;
  endif
  if( NOTstellsym ) then
   ;RlBCAST( iRbs(1:mn,0:Mvol), (Mvol+1)*mn, 0 ) ! only required for ii > 1 ;
   if( Igeometry.eq.3 ) then
    RlBCAST( iZbc(1:mn,0:Mvol), (Mvol+1)*mn, 0 )
   endif
  endif

  if( Igeometry.eq.1 .or. Igeometry.eq.2 ) then
   ;iRbc(1:mn,0) = zero ! innermost volume must be trivial; this is used in volume; innermost interface is coordinate axis;
   if( NOTstellsym ) then
    iRbs(1:mn,0) = zero ! innermost volume must be trivial; this is used in volume;
   endif
  endif

  if( Igeometry.eq.3 ) then
   iZbs(1,0:Mvol) = zero ! Zbs_{m=0,n=0} is irrelevant;
  endif
  if( NOTstellsym) then
   iRbs(1,0:Mvol) = zero ! Rbs_{m=0,n=0} is irrelevant;
  endif

  if ( Igeometry.eq.1 .and. Lreflect.eq.1) then ! reflect upper and lower bound in slab, each take half the amplitude
    iRbc(2:mn,Mvol) = iRbc(2:mn,Mvol) * half
    iRbc(2:mn,0) = -iRbc(2:mn,Mvol)
   if( NOTstellsym ) then
    iRbs(2:mn,Mvol) = iRbs(2:mn,Mvol) * half
    iRbs(2:mn,0) = -iRbs(2:mn,Mvol)
   endif
  endif



  Rscale = iRbc(1,Mvol) ! this will be used to normalize the geometrical degrees-of-freedom;

  if( myid.eq.0 ) write(ounit,'("readin : ", 10x ," : myid=",i3," ; Rscale=",es22.15," ;")') myid, Rscale




  call random_seed()



  pi2nfp         = pi2 / Nfp

  pi2pi2nfp      = pi2 * pi2nfp
  pi2pi2nfphalf  = pi2 * pi2nfp * half
  pi2pi2nfpquart = pi2 * pi2nfp * quart

  Mrad  = maxval( Lrad(1:Mvol) )

  if( myid.eq.0 ) write(ounit,'("preset : ",10x," : myid=",i3," ; Mrad=",i3," : Lrad=",257(i3,",",:))') myid, Mrad, Lrad(1:Mvol)




  select case( Igeometry )
  case( 1:2)
   if( YESstellsym ) LGdof = mn
   if( NOTstellsym ) LGdof = mn        + mn-1
  case(   3)
   if( YESstellsym ) LGdof = mn + mn-1
   if( NOTstellsym ) LGdof = mn + mn-1 + mn-1 + mn
  end select

  NGdof = ( Mvol-1 ) * LGdof

  if( Wpreset ) then ; cput = GETTIME ; write(ounit,'("preset : ",f10.2," : myid=",i3," ; NGdof=",i9," ;")') cput-cpus, myid, NGdof
  endif




  do vvol = 0, Nvol

   if( ql(vvol).eq.0 .and. qr(vvol).eq.0 ) then ; iota(vvol) = iota(vvol)
   else                                         ; iota(vvol) = ( pl(vvol) + goldenmean * pr(vvol) ) / ( ql(vvol) + goldenmean * qr(vvol) )
   endif

   if( lq(vvol).eq.0 .and. rq(vvol).eq.0 ) then ; oita(vvol) = oita(vvol)
   else                                         ; oita(vvol) = ( lp(vvol) + goldenmean * rp(vvol) ) / ( lq(vvol) + goldenmean * rq(vvol) )
   endif

   if( Wpreset .and. myid.eq.0 ) write(ounit,1002) vvol, pl(vvol), ql(vvol), pr(vvol), qr(vvol), iota(vvol), lp(vvol), lq(vvol), rp(vvol), rq(vvol), oita(vvol)

1002 format("preset : ",10x," :      ",3x," ; transform : ",i3," : (",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ",&
                                                                  "(",i3," /",i3," ) * (",i3," /",i3," ) = ",f18.15," ; ")

  enddo




  SALLOCATE( dtflux, (1:Mvol), zero )
  SALLOCATE( dpflux, (1:Mvol), zero )

  select case( Igeometry )
  case( 1   ) ; dtflux(1) = tflux(1) ; dpflux(1) = pflux(1) ! Cartesian              ; this is the "inverse" operation defined in xspech; 09 Mar 17;
  case( 2:3 ) ; dtflux(1) = tflux(1) ; dpflux(1) = zero     ! cylindrical or toroidal;
  end select

  dtflux(2:Mvol) = tflux(2:Mvol) - tflux(1:Mvol-1)
  dpflux(2:Mvol) = pflux(2:Mvol) - pflux(1:Mvol-1)

  dtflux(1:Mvol) = dtflux(1:Mvol) * phiedge / pi2
  dpflux(1:Mvol) = dpflux(1:Mvol) * phiedge / pi2




if (Lconstraint.EQ.3) then

  mu(1) = Ivolume(1) / (tflux(1) * phiedge)

  do vvol = 2, Mvol
    mu(vvol) = (Ivolume(vvol) - Ivolume(vvol-1)) / ((tflux(vvol) - tflux(vvol-1)) * phiedge)
  end do

endif




  SALLOCATE( sweight, (1:Mvol), zero )
  do vvol = 1, Mvol ; sweight(vvol) = upsilon * (vvol*one/Nvol)**wpoloidal ! 11 July 18;
  enddo





  SALLOCATE( TT, (0:Mrad,0:1,0:1), zero )
  SALLOCATE(RTT, (0:Lrad(1),0:Mpol,0:1,0:1), zero )
  SALLOCATE(RTM, (0:Lrad(1),0:Mpol), zero )

  call get_cheby( -one, Mrad, TT(:,0,:))
  call get_cheby( one , Mrad, TT(:,1,:))

  call get_zernike( zero, Lrad(1), Mpol, RTT(:,:,0,:))
  call get_zernike( one, Lrad(1), Mpol, RTT(:,:,1,:))
  call get_zernike_rm(zero, Lrad(1), Mpol, RTM(:,:))




  SALLOCATE( ImagneticOK, (1:Mvol), .false. )


  SALLOCATE( ki, (1:mn,0:1), 0 )
  SALLOCATE( kija, (1:mn,1:mn,0:1), 0 )
  SALLOCATE( kijs, (1:mn,1:mn,0:1), 0 )

  do ii = 1, mn  ; mi =  im(ii) ; ni =  in(ii)

    call getimn(lMpol, lNtor, Nfp, mi, ni, kk)
    if (kk.gt.0) then
      if( mi.eq.0 .and. ni.eq.0 ) then ; ki(ii,0:1) = (/ kk, 1 /)
      else                             ; ki(ii,0:1) = (/ kk, 2 /)
      endif
    endif

    do jj = 1, mn  ; mj =  im(jj) ; nj =  in(jj) ; mimj = mi + mj ; ninj = ni + nj !   adding   ; 17 Dec 15;

      call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
      if (kk.gt.0) then
        if( mimj.eq.0 .and. ninj.eq.0 ) then ; kija(ii,jj,0:1) = (/ kk, 1 /)
        else                                 ; kija(ii,jj,0:1) = (/ kk, 2 /)
        endif
      endif
      ;                                           ; mimj = mi - mj ; ninj = ni - nj ! subtracting; 17 Dec 15;

      if( mimj.gt.0 .or. ( mimj.eq.0 .and. ninj.ge.0 ) ) then
        call getimn(lMpol, lNtor, Nfp, mimj, ninj, kk)
        if (kk.gt.0) then
          if( mimj.eq.0 .and. ninj.eq.0 ) then ; kijs(ii,jj,0:1) = (/ kk, 1 /)
          else                                 ; kijs(ii,jj,0:1) = (/ kk, 2 /)
          endif
        endif
      else
        call getimn(lMpol, lNtor, Nfp, -mimj, -ninj, kk)
        if (kk.gt.0) then
          ;                                    ; kijs(ii,jj,0:1) = (/ kk , - 2 /) ! only the sine modes need the sign factor; 17 Dec 15;
        endif
      endif

    enddo ! end of do jj; 29 Jan 13;

  enddo ! end of do ii; 29 Jan 13;



  if( Igeometry.eq.2 ) then ! standard cylindrical; 04 Dec 14;

   SALLOCATE( djkp, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;
   SALLOCATE( djkm, (1:mn,1:mn), 0 ) ! only used in volume; trigonometric identities; 04 Dec 14;

   do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
    do jj = 1, mn ; mj = im(jj) ; nj = in(jj)
     if( mi-mj.eq.0 .and. ni-nj.eq.0 ) djkp(ii,jj) = 1
     if( mi+mj.eq.0 .and. ni+nj.eq.0 ) djkm(ii,jj) = 1
    enddo
   enddo

  endif ! end of if( Igeometry.eq.2 ) ; 04 Dec 14;




  SALLOCATE( iotakkii, (1:mn      ), 0 ) ! used to identify matrix elements in straight-field-line angle transformation;

  SALLOCATE( iotaksub, (1:mn,1:mns), 0 )
  SALLOCATE( iotaksgn, (1:mn,1:mns), 0 )
  SALLOCATE( iotakadd, (1:mn,1:mns), 0 )

  do kk = 1, mn ; mk = im(kk) ; nk = in(kk)

    call getimn(sMpol, sNtor, Nfp, mk, nk, ii)
    if (ii.gt.0) iotakkii(kk) = ii

    do jj = 1, mns ; mj = ims(jj) ; nj = ins(jj)

      mkmj = mk - mj ; nknj = nk - nj

      if( mkmj.gt.0 .or. ( mkmj.eq.0 .and. nknj.ge.0 ) ) then

        call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
        if (ii.gt.0) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  1
        endif

      else

        call getimn(sMpol, sNtor, Nfp, -mkmj, -nknj, ii)
        if (ii.gt.0) then ; iotaksub(kk,jj) = ii ; iotaksgn(kk,jj) =  -1
        endif

      endif

      mkmj = mk + mj ; nknj = nk + nj

      call getimn(sMpol, sNtor, Nfp, mkmj, nknj, ii)
      if (ii.gt.0) then ; iotakadd(kk,jj) = ii
      endif

    enddo ! end of do jj; 29 Jan 13;

  enddo ! end of do kk; 29 Jan 13;






  SALLOCATE( IPDt, (1:Mvol), zero)
  SALLOCATE( IPDtDpf, (1:Mvol-1, 1:Mvol-1), zero)

  SALLOCATE( cheby, (0:Mrad,0:2), zero )
  SALLOCATE( zernike, (0:Lrad(1), 0:Mpol, 0:2), zero )




  SALLOCATE( Iquad, (1:Mvol), 0 ) ! 16 Jan 13;

  do vvol = 1, Mvol

   LREGION(vvol)

   if( Nquad.gt.0 ) then ;            Iquad(vvol) =                         Nquad
   else
    if(      Lcoordinatesingularity ) Iquad(vvol) = Mpol + 2 * Lrad(vvol) - Nquad ! NEED TO REVISE REGULARIZATION FACTORS; 26 Feb 13;
    if( .not.Lcoordinatesingularity ) Iquad(vvol) =        2 * Lrad(vvol) - Nquad
   endif

  enddo ! end of do vvol; 18 Feb 13;


  maxIquad = maxval(Iquad(1:Mvol))

  SALLOCATE( gaussianweight   , (1:maxIquad,1:Mvol), zero ) ! perhaps it would be neater to make this a structure; 26 Jan 16;
  SALLOCATE( gaussianabscissae, (1:maxIquad,1:Mvol), zero )

  do vvol = 1, Mvol

   lquad = Iquad(vvol)

   call gauleg( lquad, gaussianweight(1:lquad,vvol), gaussianabscissae(1:lquad,vvol), igauleg ) ! JAB; 28 Jul 17

   if( myid.eq.0 ) then
    cput= GETTIME
    select case( igauleg ) !                                                  123456789012345
    case( 0 )    ; if( Wpreset ) write(ounit,1000) cput-cpus, vvol, igauleg, "success        ", gaussianabscissae(1:lquad,vvol)
    case( 1 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "failed         ", gaussianabscissae(1:lquad,vvol)
    case( 2 )    ;               write(ounit,1000) cput-cpus, vvol, igauleg, "input error    ", gaussianabscissae(1:lquad,vvol)
    case default ;               write(ounit,1000) cput-cpus, vvol, igauleg, "weird          ", gaussianabscissae(1:lquad,vvol)
     FATAL( preset, .true., weird ifail returned by gauleg )
    end select
    ;            ; if( Wpreset ) write(ounit,1001)                                              gaussianweight(1:lquad,vvol)
   endif

1000 format("preset : ",f10.2," : lvol=",i3," ; igauleg=",i5," ; ",a15," ; abscissae ="99f09.05)
1001 format("preset : ", 10x ," :      ",3x,"           ",5x,"   ",15x," ; weights   ="99f09.05)

  enddo ! end of do vvol;  7 Mar 13;




  LBsequad = .false.
  LBnewton = .false.
  LBlinear = .false.

  if( LBeltrami.eq.1 .or. LBeltrami.eq.3 .or. LBeltrami.eq.5 .or. LBeltrami.eq.7 ) LBsequad = .true.
  if( LBeltrami.eq.2 .or. LBeltrami.eq.3 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBnewton = .true.
  if( LBeltrami.eq.4 .or. LBeltrami.eq.5 .or. LBeltrami.eq.6 .or. LBeltrami.eq.7 ) LBlinear = .true.

  if (Lconstraint .eq. 2) then
    FATAL( preset, Lfreebound.eq.1, The combination of helicity constraint and free boundary is under construction )
    if (Igeometry .eq. 3 .and. myid.eq.0) then
      write(ounit, *) 'WARNING: The Hessian matrix needs further review for Igeometry = 3'
      write(ounit, *) '         However, it can still serve the purpose of Lfindzero = 2'
    endif
  endif

  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("preset : ",f10.2," : LBsequad="L2" , LBnewton="L2" , LBlinear="L2" ;")')cput-cpus, LBsequad, LBnewton, LBlinear
  endif




  SALLOCATE( BBweight, (1:mn), opsilon * exp( - escale * ( im(1:mn)**2 + (in(1:mn)/Nfp)**2 ) ) )

  if( myid.eq.0 .and. escale.gt.small ) then
   do ii = 1, mn ; write(ounit,'("preset : " 10x " : myid="i3" ; ("i3","i3") : BBweight="es13.5" ;")') myid, im(ii), in(ii)/Nfp, BBweight(ii)
   enddo
  endif




  SALLOCATE( mmpp, (1:mn), zero )

  do ii = 1, mn ; mi = im(ii)

   if( mi.eq.0 ) then ; mmpp(ii) = zero
   else               ; mmpp(ii) = mi**pcondense
   endif ! end of if( mi.eq.0 ) ; 11 Aug 14;

  enddo ! end of do ii; 08 Nov 13;




  SALLOCATE( NAdof, (1:Mvol          ), 0 ) ! Beltrami degrees-of-freedom in each annulus;
  SALLOCATE( Nfielddof,(1:Mvol       ), 0 ) ! Beltrami degrees-of-freedom in each annulus, field only;
  SALLOCATE( NdMASmax, (1:Mvol       ), 0 ) ! The maximum size of sparse matrix for GMRES preconditioning;
  SALLOCATE( NdMAS   , (1:Mvol       ), 0 ) ! The actual size of sparse matrix for GMRES preconditioning;

  NALLOCATE( Ate  , (1:Mvol,-2:2,1:mn)    ) ! recall that this is type:sub-grid; 31 Jan 13;
  NALLOCATE( Aze  , (1:Mvol,-2:2,1:mn)    ) ! -2 : for use of matrix-free solver ; -1 : for use of force gradient
  NALLOCATE( Ato  , (1:Mvol,-2:2,1:mn)    ) !  0 : normal data
  NALLOCATE( Azo  , (1:Mvol,-2:2,1:mn)    ) ! 1:2: use to compute derivative w.r.t. fluxes

  SALLOCATE( Fso  , (1:Mvol,     1:mn), 0 ) ! these will become redundant if/when Lagrange multipliers are used to enforce bounday constraints; 26 Jan 16;
  SALLOCATE( Fse  , (1:Mvol,     1:mn), 0 )

  SALLOCATE( Lma  , (1:Mvol,     1:mn), 0 ) ! degree of freedom index; for Lagrange multiplier; 08 Feb 16;
  SALLOCATE( Lmb  , (1:Mvol,     1:mn), 0 )
  SALLOCATE( Lmc  , (1:Mvol,     1:mn), 0 ) ! only need Lmc(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmd  , (1:Mvol,     1:mn), 0 ) ! only need Lmd(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lme  , (1:Mvol,     1:mn), 0 ) ! only need Lme(2:mn) ;
  SALLOCATE( Lmf  , (1:Mvol,     1:mn), 0 ) ! only need Lmf(2:mn) ; only for NOTstellsym; 08 Feb 16;
  SALLOCATE( Lmg  , (1:Mvol,     1:mn), 0 ) ! only need Lmg(1   ) ;
  SALLOCATE( Lmh  , (1:Mvol,     1:mn), 0 ) ! only need Lmh(1   ) ;

  SALLOCATE( Lmavalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmbvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmcvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmdvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmevalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmfvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmgvalue, (1:Mvol,     1:mn), zero )
  SALLOCATE( Lmhvalue, (1:Mvol,     1:mn), zero )

  do vvol = 1, Mvol

   LREGION(vvol)

   if( Lcoordinatesingularity ) then
    zerdof = 0                                       ! count Zernike degree of freedom 30 Jun 19
    do ii = 2, Mpol                                  ! for m>1
     do jj = ii, Lrad(vvol), 2
      zerdof = zerdof + 2 * ntor + 1                 ! plus and minus sign for n>1, unique for n==0
      if( NOTstellsym ) zerdof = zerdof + 2*ntor + 1 ! plus and minus sign for n
     enddo
    enddo
    zerdof = zerdof * 2                              ! we have one for At and one for Az

    do jj = 0, Lrad(vvol), 2                         ! for m==0
     zerdof = zerdof + ntor + 1                      ! minus sign for n, Aze
     if (jj .ge. 2) zerdof = zerdof + ntor + 1       ! minus sign for n, Ate, without l=0 due to recombination

     if( NOTstellsym ) then
      zerdof = zerdof + ntor                         ! sin component minus sign for n, Azo
      if (jj .ge. 2) zerdof = zerdof + ntor          ! minus sign for n, Ato, without l=0 due to recombination
     endif
    enddo

    if (Mpol .ge. 1) then ! for m==1
      do jj = 1, Lrad(vvol), 2
        zerdof = zerdof + 2 * ntor + 1                  ! minus and plus sign for n, Aze
        if (jj .ge. 2) zerdof = zerdof + 2 * ntor + 1   ! minus sign for n, Ate, without l=0 due to recombination

        if( NOTstellsym ) then
          zerdof = zerdof + 2 * ntor + 1                 ! sin component minus and plus sign for n, Azo
          if (jj .ge. 2) zerdof = zerdof + 2 * ntor + 1  ! minus and plus sign for n, Ato, without l=0 due to recombination
        endif
      enddo
    endif

    Nfielddof(vvol) = zerdof
    if( YESstellsym ) NAdof(vvol) = zerdof                               + mn        + Ntor+1        + mn-1        + 1 + 0
    if( NOTstellsym ) NAdof(vvol) = zerdof                               + mn + mn-1 + Ntor+1 + Ntor + mn-1 + mn-1 + 1 + 0 ! this is broken at the moment

    NAdof(vvol) = NAdof(vvol) - (ntor + 1)
    if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - ntor

    if (Mpol .ge. 1) then
      NAdof(vvol) = NAdof(vvol) - (2 * ntor + 1)
      if (NOTstellsym) NAdof(vvol) = NAdof(vvol) - (2 * ntor + 1)
    endif

   else ! .not.Lcoordinatesingularity;                                     a    c      b        d      e      f      g   h
    if( YESstellsym ) NAdof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)    )                            + mn-1        + 1 + 1
    if( NOTstellsym ) NAdof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)    )                            + mn-1 + mn-1 + 1 + 1

    if( YESstellsym ) Nfielddof(vvol) = 2 * ( mn        ) * ( Lrad(vvol)    )
    if( NOTstellsym ) Nfielddof(vvol) = 2 * ( mn + mn-1 ) * ( Lrad(vvol)    )


   endif ! end of if( Lcoordinatesingularity );

   do ii = 1, mn ! loop over Fourier harmonics;

    do ideriv = -2, 2 ! loop over derivatives; 14 Jan 13;

     SALLOCATE( Ate(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Aze(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Ato(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )
     SALLOCATE( Azo(vvol,ideriv,ii)%s, (0:Lrad(vvol)), zero )

    enddo ! end of do ideriv;

    ;  ideriv =  0

     SALLOCATE( Ate(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 ) ! degree of freedom index; 17 Jan 13;
     SALLOCATE( Aze(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Ato(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )
     SALLOCATE( Azo(vvol,ideriv,ii)%i, (0:Lrad(vvol)), 0 )

   enddo ! end of do ii;


   select case( Linitgues ) ! for iterative solver of the Beltrami fields, an initial guess is required; 11 Mar 16;
   case( 0 )    ;
   case( 1 )    ; Ate(vvol,0,1)%s(0:1) = dtflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    ;           ; Aze(vvol,0,1)%s(0:1) = dpflux(vvol) * half ! this is an integrable approximation; NEEDS CHECKING; 26 Feb 13;
    if (Lcoordinatesingularity) then
    ;           ; Ate(vvol,0,1)%s(2) = dtflux(vvol) * half * half
    endif
   case( 2 )    ;                                            ! will call ra00aa below to read initial vector potential from file;
   case( 3 )    ;                                            ! the initial guess will be randomized, maximum is maxrndgues; 5 Mar 19;
    do ii = 1, mn ! loop over Fourier harmonics;

     do ideriv = -2, 2 ! loop over derivatives; 14 Jan 13;

      call random_number(Ate(vvol,ideriv,ii)%s)
      call random_number(Aze(vvol,ideriv,ii)%s)
      Ate(vvol,ideriv,ii)%s = Ate(vvol,ideriv,ii)%s * maxrndgues
      Aze(vvol,ideriv,ii)%s = Aze(vvol,ideriv,ii)%s * maxrndgues
      if (.not. YESstellsym) then
       call random_number(Ato(vvol,ideriv,ii)%s)
       call random_number(Azo(vvol,ideriv,ii)%s)
       Ato(vvol,ideriv,ii)%s = Ato(vvol,ideriv,ii)%s * maxrndgues
       Azo(vvol,ideriv,ii)%s = Azo(vvol,ideriv,ii)%s * maxrndgues
      endif

     enddo ! end of do ideriv;

    enddo ! end of do ii;

   end select

   idof = 0 ! degree of freedom index; reset to 0 in each volume;

   if( Lcoordinatesingularity ) then

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

     do ll = 0, Lrad(vvol)
      if (ll>=mi .and. mod(mi+ll,2)==0)then
      if (.not.((ll==0.and.mi==0).or.(ll==1.and.mi==1))) then
                                            ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof ! Zernike 30 Jun 19
      endif
      ;                                     ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( NOTstellsym .and. ii.gt.1 ) then
        if (.not.((ll==0.and.mi==0).or.(ll==1.and.mi==1))) then
                                            ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof ! Zernike 30 Jun 19
        endif
       ;                                    ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif ! NOTstellsym
      endif ! Zernike
     enddo ! end of do ll; 17 Jan 13;

    enddo ! end of do ii

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
     if ( mi.ne.0 .and. mi.ne.1     )  then ; idof = idof + 1 ; Lma(vvol,  ii)       = idof
     endif
     if(  mi.eq.0                   ) then ; idof = idof + 1 ; Lmb(vvol,  ii)       = idof ! 18 May 16;
     endif
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
     endif
     if( NOTstellsym ) then
      if(  mi.ne.0 .and. mi.ne.1    ) then ; idof = idof + 1 ; Lmc(vvol,  ii)       = idof ! 18 May 16;
      endif
      if(  ii.gt.1                  ) then ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof ! 18 May 16;
     endif
     if(  ii.gt.1 .and. mi.eq.0     ) then ; idof = idof + 1 ; Lmd(vvol,  ii)       = idof ! 18 May 16;
     endif
     endif ! end of if( NOTstellsym ) ; 19 Jul 16;

    enddo ! end of do ii; 25 Jan 13;

    FATAL( preset, idof.ne.NAdof(vvol), need to count Beltrami degrees-of-freedom more carefully  for coordinate singularity )
    FATAL( preset, (idof+1)**2.ge.HUGE(idof)), NAdof too big, should be smaller than maximum of int32 type )

   else ! .not.Lcoordinatesingularity;

    do ii = 1, mn
     do ll = 1, Lrad(vvol)                 ; idof = idof + 1 ; Ate(vvol,0,ii)%i(ll) = idof
      ;                                    ; idof = idof + 1 ; Aze(vvol,0,ii)%i(ll) = idof
      if( ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Ato(vvol,0,ii)%i(ll) = idof
       ;                                   ; idof = idof + 1 ; Azo(vvol,0,ii)%i(ll) = idof
      endif
     enddo ! end of do ll; 08 Feb 16;
    enddo

    do ii = 1, mn
     if(  ii.gt.1                   ) then ; idof = idof + 1 ; Lme(vvol,  ii)       = idof
     endif
     if(  ii.gt.1 .and. NOTstellsym ) then ; idof = idof + 1 ; Lmf(vvol,  ii)       = idof
     endif
     if(  ii.eq.1                   ) then ; idof = idof + 1 ; Lmg(vvol,  ii)       = idof
      ;                                    ; idof = idof + 1 ; Lmh(vvol,  ii)       = idof
     endif
    enddo ! end of do ii; 25 Jan 13;


    FATAL( preset, idof.ne.NAdof(vvol), need to count degrees-of-freedom more carefully for new matrix )
    FATAL( preset, (idof+1)**2.ge.HUGE(idof)), NAdof too big, should be smaller than maximum of int32 type )

   endif ! end of if( Lcoordinatesingularity ) ;

   FATAL( preset, idof.ne.NAdof(vvol), impossible logic )

   do ii = 1, mn
      do jj = 0, Lrad(vvol)
        if (Ate(vvol,0,ii)%i(jj) == 0) Ate(vvol,0,ii)%s(jj) = zero
        if (Aze(vvol,0,ii)%i(jj) == 0) Aze(vvol,0,ii)%s(jj) = zero
        if (.not. YESstellsym) then
          if (Ato(vvol,0,ii)%i(jj) == 0) Azo(vvol,0,ii)%s(jj) = zero
          if (Azo(vvol,0,ii)%i(jj) == 0) Azo(vvol,0,ii)%s(jj) = zero
        end if
      end do !jj
   end do !ii


  enddo ! end of do vvol = 1, Nvol loop;



  if( Linitgues.eq.2 ) then ; WCALL( preset, ra00aa, ('R') )  ! read initial guess for Beltrami field from file; 02 Jan 15;
  endif



  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')
   write(ounit,'("preset : ",f10.2," : Nquad="i4" ; mn="i5" ; NGdof="i6" ; NAdof="16(i6",")" ...")') cput-cpus, Nquad, mn, NGdof, NAdof(1:min(Mvol,16))
  endif





  Nt = max( Ndiscrete*4*Mpol, 1 ) ; Nz = max( Ndiscrete*4*Ntor, 1 ) ; Ntz = Nt*Nz ; soNtz = one / sqrt( one*Ntz ) ! exaggerated discrete resolution;


  ;                  ; hNt = Nt / 2
  if( Nz.gt.1 ) then ; hNz = Nz / 2
  else               ; hNz = 0
  endif

  if( myid.eq.0 ) then ! 17 Oct 12;
   cput = GETTIME
   write(ounit,'("preset : ", 10x ," : ")')
   write(ounit,'("preset : ",f10.2," : Nt="i6" ; Nz="i6" ; Ntz="i9" ;")') cput-cpus, Nt, Nz, Ntz
  endif

  SALLOCATE( iRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; ! 18 Jul 14;
  SALLOCATE( iZij, (1:Ntz,0:Mvol), zero ) !
  SALLOCATE( dRij, (1:Ntz,1:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( dZij, (1:Ntz,1:Mvol), zero )
  SALLOCATE( tRij, (1:Ntz,0:Mvol), zero ) ! interface geometry in real space; poloidal derivative; ! 18 Jul 14;
  SALLOCATE( tZij, (1:Ntz,0:Mvol), zero )

  SALLOCATE(   Rij, (1:Ntz,0:3,0:3    ), zero ) ! these are used for inverse fft to reconstruct real space geometry from interpolated Fourier harmonics;
  SALLOCATE(   Zij, (1:Ntz,0:3,0:3    ), zero )
  SALLOCATE(   sg , (1:Ntz,0:3        ), zero )
  SALLOCATE( guvij, (1:Ntz,0:3,0:3,-1:3), zero ) ! need this on higher resolution grid for accurate Fourier decomposition;
  SALLOCATE( gvuij, (1:Ntz,0:3,0:3    ), zero ) ! need this on higher resolution grid for accurate Fourier decomposition; 10 Dec 15;

  if ((Lfindzero .eq. 2) .or. (Lcheck.eq.5 .or. LHevalues .or. LHevectors .or. LHmatrix .or. Lperturbed.eq.1)) then
    SALLOCATE( dRadR, (1:mn,0:1,0:1,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRadZ, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadR, (1:mn,0:1,0:1,1:mn), zero )
    SALLOCATE( dZadZ, (1:mn,0:1,0:1,1:mn), zero )

    SALLOCATE( dRodR, (1:Ntz,0:3,1:mn), zero ) ! calculated in rzaxis; 19 Sep 16;
    SALLOCATE( dRodZ, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodR, (1:Ntz,0:3,1:mn), zero )
    SALLOCATE( dZodZ, (1:Ntz,0:3,1:mn), zero )
  endif



  SALLOCATE( goomne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( goomno, (0:mne, maxIquad), zero )
  SALLOCATE( gssmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gssmno, (0:mne, maxIquad), zero )
  SALLOCATE( gstmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gstmno, (0:mne, maxIquad), zero )
  SALLOCATE( gszmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gszmno, (0:mne, maxIquad), zero )
  SALLOCATE( gttmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gttmno, (0:mne, maxIquad), zero )
  SALLOCATE( gtzmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gtzmno, (0:mne, maxIquad), zero )
  SALLOCATE( gzzmne, (0:mne, maxIquad), zero ) ! workspace for Fourier decomposition of metric terms;
  SALLOCATE( gzzmno, (0:mne, maxIquad), zero )



  SALLOCATE( ijreal, (1:Ntz), zero ) ! real space grid;
  SALLOCATE( ijimag, (1:Ntz), zero )
  SALLOCATE( jireal, (1:Ntz), zero )
  SALLOCATE( jiimag, (1:Ntz), zero )

  SALLOCATE( jkreal, (1:Ntz), zero )
  SALLOCATE( jkimag, (1:Ntz), zero )
  SALLOCATE( kjreal, (1:Ntz), zero )
  SALLOCATE( kjimag, (1:Ntz), zero )

  SALLOCATE( cplxin,  (1:Nt,1:Nz,1), zero )
  SALLOCATE( cplxout, (1:Nt,1:Nz,1), zero )

  planf = fftw_plan_dft_2d( Nz, Nt, cplxin(:,:,1), cplxout(:,:,1), FFTW_FORWARD,  FFTW_MEASURE + FFTW_DESTROY_INPUT )
  planb = fftw_plan_dft_2d( Nz, Nt, cplxin(:,:,1), cplxout(:,:,1), FFTW_BACKWARD, FFTW_MEASURE + FFTW_DESTROY_INPUT )

  SALLOCATE( efmn, (1:mne), zero ) ! Fourier harmonics workspace; 24 Apr 13;
  SALLOCATE( ofmn, (1:mne), zero )
  SALLOCATE( cfmn, (1:mne), zero )
  SALLOCATE( sfmn, (1:mne), zero )
  SALLOCATE( evmn, (1:mne), zero )
  SALLOCATE( odmn, (1:mne), zero )
  SALLOCATE( comn, (1:mne), zero )
  SALLOCATE( simn, (1:mne), zero )




  SALLOCATE( gteta, (1:Ntz), zero )
  SALLOCATE( gzeta, (1:Ntz), zero )

  SALLOCATE( cosi, (1:Ntz,1:mn), zero )
  SALLOCATE( sini, (1:Ntz,1:mn), zero )

  FATAL( preset, Nz.eq.0, illegal division )
  FATAL( preset, Nt.eq.0, illegal division )

  do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier harmonics;

  do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
    do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt ; arg = mi * teta - ni * zeta
    gteta(jk) = teta
    gzeta(jk) = zeta
    cosi(jk,ii) = cos(arg)
    sini(jk,ii) = sin(arg)
    enddo
  enddo

  enddo ! end of do ii; 13 May 13;



  if( Igeometry.eq.3 .and. iRbc(1,0).lt.small ) then ! have not yet assigned coordinate axis; see global;readin for user-supplied Rac, Zas, etc. ; 19 Jul 16;

    write(*,*) "Finding initial axis..."
   select case( Linitialize )
   case( :-1 ) ; vvol = Nvol + Linitialize
   case(   0 ) ; vvol =    1 ! this is really a dummy; no interpolation of interface geometry is required; packxi calls rzaxis with lvol=1; 19 Jul 16;
   case(   1 ) ; vvol = Nvol
   case(   2 ) ; vvol = Mvol
   end select

   WCALL( preset, rzaxis, ( Mvol, mn, iRbc(1:mn,0:Mvol), iZbs(1:mn,0:Mvol), iRbs(1:mn,0:Mvol), iZbc(1:mn,0:Mvol), vvol, .false. ) ) ! set coordinate axis; 19 Jul 16;

  endif ! end of if( Igeometry.eq.3 ) then ; 19 Jul 16;




  SALLOCATE( psifactor, (1:mn,1:Mvol), zero )
  SALLOCATE( inifactor, (1:mn,1:Mvol), zero )

  psifactor(1:mn,1:Mvol) = one
  inifactor(1:mn,1:Mvol) = one

  select case( Igeometry )

  case( 1 )

   psifactor(1:mn,1:Nvol) = one

  case( 2 )

   do vvol = 1, Nvol
    do ii = 1, mn
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = tflux(vvol)**(          +half) ! 28 Jan 15;
     else                   ; psifactor(ii,vvol) = tflux(vvol)**(halfmm(ii)-half) ! 28 Jan 15;
     endif
    enddo
   enddo

  case( 3 )

   do vvol = 1, Nvol
    do ii = 1, mn
     if( im(ii).eq.0 ) then ; psifactor(ii,vvol) = Rscale * tflux(vvol)**zero       ! 08 Feb 16;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**half       ! 17 Dec 18;
     else                   ; psifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 29 Apr 15;
                            ; inifactor(ii,vvol) = Rscale * tflux(vvol)**halfmm(ii) ! 17 Dec 18
     endif
    enddo
   enddo

  case default

   FATAL( readin, .true., invalid Igeometry for construction of psifactor )

  end select



  if( Linitialize.ne.0 ) then ! interpolate / extrapolate interior interface geometry; 19 Jul 16;

   select case( Igeometry )

   case( 1 ) ! Cartesian; 29 Apr 14;


    do vvol = 1, Nvol
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Mvol) * tflux(vvol) / tflux(Mvol) ! 14 Apr 17;
     endif
    enddo

   case( 2 ) ! cylindrical - standard; 20 Apr 13;

    FATAL( preset, Linitialize.ne.1, geometrical initialization under construction for cylindrical )

    do vvol = 1, Nvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,Nvol) * psifactor(1:mn,vvol)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,Nvol) * psifactor(2:mn,vvol)
     endif
    enddo

   case( 3 ) ! toroidal; 20 Apr 13;

    FATAL( preset, Linitialize.lt.0, geometrical initialization under construction for toroidal ) ! see commented-out source below; 19 Jul 16;

    lvol = Nvol-1 + Linitialize

    do vvol = 1, lvol-1
     ;iRbc(1:mn,vvol) = iRbc(1:mn,0) + ( iRbc(1:mn,lvol) - iRbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     ;iZbs(2:mn,vvol) = iZbs(2:mn,0) + ( iZbs(2:mn,lvol) - iZbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
     if( NOTstellsym ) then
      iRbs(2:mn,vvol) = iRbs(2:mn,0) + ( iRbs(2:mn,lvol) - iRbs(2:mn,0) ) * ( inifactor(2:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(2:mn)
      iZbc(1:mn,vvol) = iZbc(1:mn,0) + ( iZbc(1:mn,lvol) - iZbc(1:mn,0) ) * ( inifactor(1:mn,vvol) / Rscale ) / tflux(lvol)**halfmm(1:mn)
     endif
    enddo


   end select ! matches select case( Igeometry ); 19 Jul 16;

  endif ! matches if( Linitialize.ne.0 ) then; 19 Jul 16;




  SALLOCATE( Bsupumn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \theta on boundary; required for virtual casing;
  SALLOCATE( Bsupvmn, (1:Nvol,0:1,1:mn), zero ) ! Fourier components of {\bf B}\cdot\nabla \zeta  on boundary;




  SALLOCATE( diotadxup, (0:1,-1:2,1:Mvol), zero ) ! measured rotational transform on inner/outer interfaces in each annulus;
  SALLOCATE( dItGpdxtp, (0:1,-1:2,1:Mvol), zero ) ! measured plasma and linking currents                                   ;

  SALLOCATE( glambda, (1:Ntz+1,0:2,0:1,1:Mvol), zero ) ! save initial guesses for iterative calculation of rotational-transform; 21 Apr 13;




  SALLOCATE( Bemn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Bomn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Iomn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Iemn, (1:mn,1:Mvol    ), zero )
  SALLOCATE( Somn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Semn, (1:mn,1:Mvol,0:1), zero )
  SALLOCATE( Pomn, (1:mn,1:Mvol,0:2), zero )
  SALLOCATE( Pemn, (1:mn,1:Mvol,0:2), zero )

  SALLOCATE( BBe , (1:Mvol-1), zero )
  SALLOCATE( IIo , (1:Mvol-1), zero )
  SALLOCATE( BBo , (1:Mvol-1), zero )
  SALLOCATE( IIe , (1:Mvol-1), zero )



  SALLOCATE( Btemn, (1:mn,0:1,1:Mvol), zero ) ! these are declared in global, calculated in sc00aa, broadcast in xspech, and written to file in hdfint;
  SALLOCATE( Bzemn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Btomn, (1:mn,0:1,1:Mvol), zero )
  SALLOCATE( Bzomn, (1:mn,0:1,1:Mvol), zero )

  SALLOCATE( Bloweremn, (1:mn, 3), zero) ! these are declared in global, calculated in getbco, used in mtrxhs
  SALLOCATE( Bloweromn, (1:mn, 3), zero)



  SALLOCATE( vvolume    , (1:Mvol), zero ) ! volume integral of \sqrt g;
  SALLOCATE( lBBintegral, (1:Mvol), zero ) ! volume integral of B.B    ;
  SALLOCATE( lABintegral, (1:Mvol), zero ) ! volume integral of A.B    ;



  if( YESstellsym ) lmns = 1 + (mns-1)           ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13;
  if( NOTstellsym ) lmns = 1 + (mns-1) + (mns-1) ! number of independent degrees-of-freedom in angle transformation; 30 Jan 13;

  SALLOCATE( dlambdaout, (1:lmns,1:Mvol,0:1), zero )

  if( Lconstraint .EQ. 3) then
    Localconstraint = .false.
  else
    Localconstraint = .true.
  endif
  




end subroutine preset


