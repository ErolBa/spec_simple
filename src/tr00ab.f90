
subroutine tr00ab( lvol, mn, NN, Nt, Nz, iflag, ldiota ) ! construct straight-field line magnetic coordinates;



  use constants, only : zero, third, half, one, two, pi2, goldenmean

  use numerical, only : vsmall, small, machprec, sqrtmachprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wtr00ab, Nvol, Lrad, Mpol, Ntor, &
                        Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol

  use cputiming, only : Ttr00ab

  use allglobal, only : ncpu, cpus, myid, MPI_COMM_SPEC, &
                        pi2nfp, &
                        Mvol, im, in, mns, ims, ins, &
                        YESstellsym, NOTstellsym, &
                        glambda, & ! global lambda: initial guesses will be saved; 21 Apr 13;
                        Ntz, hNt, hNz, &
                        iotakkii, iotaksub, iotakadd, iotaksgn, &
                        Ate, Aze, Ato, Azo, TT, RTT, &
                        Lcoordinatesingularity, Lvacuumregion, regumm, dlambdaout



  LOCALS

  INTEGER, intent(in)  :: lvol, mn, NN, Nt, Nz, iflag

  REAL, intent(inout)  :: ldiota(0:1,-1:2)

  INTEGER              :: innout, ll, ii, jj, kk, jb, kb, mj, nj, ideriv, jderiv, id, MM, ielement, nelements, Lcurvature, idof, icon, mi, ni, imupf

  REAL                 :: lcpu, mfactor, lss, Dteta, Dzeta, rfac, tol, rnorm, omega, diotaerror!, sparsedenseerror

  REAL                 :: lAte(0:mn,-1:2), lAze(0:mn,-1:2), lAto(0:mn,-1:2), lAzo(0:mn,-1:2)



  INTEGER              :: IA, if04aaf, idgesvx, ipiv(1:NN), iwork4(1:NN)
  REAL  , allocatable  :: dmatrix(:,:,:), omatrix(:,:), FAA(:,:)
  REAL                 :: drhs(1:NN,-1:2), dlambda(1:NN,-1:2)
  REAL                 :: Rdgesvx(1:NN), Cdgesvx(1:NN), work4(1:4*NN), rcond, ferr(1), berr(1), ferr2(1:2), berr2(1:2)
  CHARACTER            :: equed

  INTEGER              :: maxitn, reqdits, extralength, lrwork, integerwork(1:2*Nt*Nz+2+1), if11def, if11zaf, if11xaf
  INTEGER              :: IAA, if04atf, if04arf
  INTEGER              :: Ndof, label(-3:Nt+2,-3:Nz+2), isym

  INTEGER              :: idgelsd, Lwork, Liwork, Irank, nlvl
  REAL                 :: sval(1:NN)
  REAL   , allocatable :: work(:)

  REAL                 ::                      Bsupt(1:Nt*Nz,-1:2), Bsupz(1:Nt*Nz,-1:2), tdot(1:Nt*Nz)
  REAL                 :: Bsubs(1:Nt*Nz,-1:2), Bsubt(1:Nt*Nz,-1:2), Bsubz(1:Nt*Nz,-1:2)

  REAL                 :: dotteta, dotzeta

  REAL   , allocatable :: rmatrix(:,:,:), rrhs(:,:), rlambda(:,:), wks1(:), wks2(:), AA(:,:)

  INTEGER              :: inz(-1:2), lnz
  INTEGER, allocatable :: irow(:,:), jcol(:,:), istr(:), iwork(:)
  REAL   , allocatable :: smatrix(:,:), srhs(:,:), slambda(:,:), swork(:)
  CHARACTER            :: duplicate*1, zeros*1, method*8, precon*1, trans*1, check*1 ! logical control of sparse routines; 20 Apr 13;

  


  do innout = 0, 1 ! loop over inner and outer interfaces;



   if( Lcoordinatesingularity .and. innout.eq.0 ) cycle ! transform on coordinate    axis     is not required              ; 20 Apr 13;
   if( Lvacuumregion          .and. innout.eq.1 ) cycle ! transform on computational boundary is not required (not defined); 20 Apr 13;



   lAte(0:mn,-1:2) = zero ! radial derivatives of vector potential evaluated at interfaces; 20 Apr 13;
   lAze(0:mn,-1:2) = zero
   lAto(0:mn,-1:2) = zero
   lAzo(0:mn,-1:2) = zero

   do ideriv = -1, 2 ; id = ideriv ! labels derivative of magnetic field wrt enclosed fluxes; 20 Apr 13;

    if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives of transform                                                        are not required; 20 Jun 14;
    if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives of transform wrt geometry                                           is  not required; 20 Jun 14;
    if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivatives of transform wrt helicity multiplier and differential poloidal flux are not required; 20 Jun 14;

    do ii = 1, mn ! loop over Fourier harmonics; 20 Apr 13;

     mi = im(ii)

     if (Lcoordinatesingularity) then
      do ll = mi, Lrad(lvol), 2 ! loop over Zernike polynomials; 01 Jul 19;

        ;lAte(ii,id) = lAte(ii,id) + Ate(lvol,id,ii)%s(ll) * RTT(ll,mi,innout,1) * half
        ;lAze(ii,id) = lAze(ii,id) - Aze(lvol,id,ii)%s(ll) * RTT(ll,mi,innout,1) * half
        if( NOTstellsym ) then
        lAto(ii,id) = lAto(ii,id) + Ato(lvol,id,ii)%s(ll) * RTT(ll,mi,innout,1) * half
        lAzo(ii,id) = lAzo(ii,id) - Azo(lvol,id,ii)%s(ll) * RTT(ll,mi,innout,1) * half
        endif

      enddo ! end of do ll; 30 Jan 13;
     else
      do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; 20 Apr 13;

        ;lAte(ii,id) = lAte(ii,id) + Ate(lvol,id,ii)%s(ll) * TT(ll,innout,1) ! compute radial derivative of vector potential; 20 Apr 13;
        ;lAze(ii,id) = lAze(ii,id) - Aze(lvol,id,ii)%s(ll) * TT(ll,innout,1)
        if( NOTstellsym ) then
        lAto(ii,id) = lAto(ii,id) + Ato(lvol,id,ii)%s(ll) * TT(ll,innout,1)
        lAzo(ii,id) = lAzo(ii,id) - Azo(lvol,id,ii)%s(ll) * TT(ll,innout,1)
        endif

      enddo ! end of do ll; 30 Jan 13;
     end if ! Lcoordinatesingularity; 01 Jul 19;
    enddo ! end of do ii; 30 Jan 13; ! Fourier harmonics, and their derivatives, have been constructed; 20 Apr 13;

    if( Lsparse.gt.0 ) then ! will construct transformation to straight-field line angle in real space; 24 Apr 13;
     call invfft( mn, im(1:mn), in(1:mn), lAte(1:mn,id), lAto(1:mn,id), lAze(1:mn,id), lAzo(1:mn,id), &
                  Nt, Nz, Bsupz(1:Ntz,id), Bsupt(1:Ntz,id) ) ! map to real space;
    endif

   enddo ! end of do ideriv; 31 Jan 13;




   if( Lsparse.eq.0 .or. Lsparse.eq.3 ) then
    SALLOCATE( dmatrix, (1:NN,1:NN,-1:2), zero )
    SALLOCATE( omatrix, (1:NN,1:NN), zero )
    SALLOCATE( FAA, (1:NN,1:NN), zero )
   endif

   if( Lsparse.eq.0 .or. Lsparse.eq.3 ) then ! Fourier transformation; 24 Apr 13;

    drhs(1:NN,-1:2) = zero

    dmatrix(1:NN,1:NN,-1:2) = zero ! initialize summation; 30 Jan 13;

    do ideriv = -1, 2

     if( iflag.eq. 1 .and. ideriv.ne.0 ) cycle ! derivatives                                                        not required; 20 Jun 14;
     if( iflag.eq. 2 .and. ideriv.lt.0 ) cycle ! derivatives wrt helicity multiplier and differential poloidal flux are required; 20 Jun 14;
     if( iflag.eq.-1 .and. ideriv.gt.0 ) cycle ! derivative  wrt geometry                                               required; 20 Jun 14;
     do kk = 1, mn

      ii = iotakkii(kk)

      ;drhs(ii      ,ideriv) = + lAze(kk,ideriv)
      if( NOTstellsym .and. kk.gt.0 ) then
       drhs(ii+mns-1,ideriv) = + lAzo(kk,ideriv)
      endif

      ;dmatrix(ii      ,1,ideriv) = lAte(kk,ideriv)
      if( NOTstellsym .and. kk.gt.0 ) then
       dmatrix(ii+mns-1,1,ideriv) = lAto(kk,ideriv)
      endif

      do jj = 2, mns ; mj = ims(jj) ; nj = ins(jj) ! this seems to ignore the non-stellarator symmetric mode; 02 Sep 14;

       ii = iotakadd(kk,jj)

       if( ii.lt.1 ) cycle

       dmatrix(ii      ,jj      ,ideriv) = dmatrix(ii      ,jj      ,ideriv) + ( - mj * lAze(kk,ideriv) + nj * lAte(kk,ideriv) ) * half
       if( NOTstellsym) then
        FATAL( tr00ab,ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN, illegal subscript ) ! THIS CAN BE DELETED EVENTUALLY;
        dmatrix(ii+mns-1,jj      ,ideriv) = dmatrix(ii+mns-1,jj      ,ideriv) + ( - mj * lAzo(kk,ideriv) + nj * lAto(kk,ideriv) ) * half
        dmatrix(ii      ,jj+mns-1,ideriv) = dmatrix(ii      ,jj+mns-1,ideriv) - ( + mj * lAzo(kk,ideriv) - nj * lAto(kk,ideriv) ) * half
        dmatrix(ii+mns-1,jj+mns-1,ideriv) = dmatrix(ii+mns-1,jj+mns-1,ideriv) + ( + mj * lAze(kk,ideriv) - nj * lAte(kk,ideriv) ) * half
       endif

      enddo ! end of do jj; 30 Jan 13;


      do jj = 2, mns ; mj = ims(jj) ; nj = ins(jj)

       ii = iotaksub(kk,jj)

       if( ii.lt.1 ) cycle

       FATAL( tr00ab,ii.gt.NN .or. jj.gt.NN, illegal subscript ) ! THIS CAN BE DELETED EVENTUALLY; 02 Sep 14;
       dmatrix(ii      ,jj      ,ideriv) = dmatrix(ii      ,jj      ,ideriv) + ( - mj * lAze(kk,ideriv) + nj * lAte(kk,ideriv) ) * half
       if( NOTstellsym) then
        FATAL( tr00ab,ii+mns-1.lt.1 .or. ii+mns-1.gt.NN .or. jj+mns-1.lt.1 .or. jj+mns-1.gt.NN, illegal subscript ) ! THIS CAN BE DELETED EVENTUALLY;
        dmatrix(ii+mns-1,jj      ,ideriv) = dmatrix(ii+mns-1,jj      ,ideriv) + ( - mj * lAzo(kk,ideriv) + nj * lAto(kk,ideriv) ) * half * iotaksgn(kk,jj)
        dmatrix(ii      ,jj+mns-1,ideriv) = dmatrix(ii      ,jj+mns-1,ideriv) + ( + mj * lAzo(kk,ideriv) - nj * lAto(kk,ideriv) ) * half
        dmatrix(ii+mns-1,jj+mns-1,ideriv) = dmatrix(ii+mns-1,jj+mns-1,ideriv) - ( + mj * lAze(kk,ideriv) - nj * lAte(kk,ideriv) ) * half * iotaksgn(kk,jj)
       endif

      enddo ! end of do jj; 30 Jan 13;

     enddo ! end of do kk; 30 Jan 13;
    enddo ! end of ideriv; 30 Jan 13;



    call DCOPY(NN*NN, dmatrix(1,1,0), 1, omatrix(1,1), 1) ! BLAS version 21 Jul 19

    do jderiv = 0, 1

     if( iflag.eq. 1 .and. jderiv.ne.0 ) cycle ! derivatives of rotational transform (wrt either enclosed-fluxes/currents/geometry) are not required;

     select case( jderiv )
     case( 0 ) ;!             drhs(1:NN, 0) = drhs(1:NN, 0)
     case( 1 ) ;
      if( iflag.eq. 2) then ; call DGEMV('N',NN,NN,-one,dmatrix(1,1, 1),NN,dlambda(1,0),1,one,drhs(1, 1),1)     ! BLAS version 21 Jul 19
       ;                    ; call DGEMV('N',NN,NN,-one,dmatrix(1,1, 2),NN,dlambda(1,0),1,one,drhs(1, 2),1)     ! BLAS version 21 Jul 19
      endif
      if( iflag.eq.-1) then ; call DGEMV('N',NN,NN,-one,dmatrix(1,1,-1),NN,dlambda(1,0),1,one,drhs(1,-1),1)     ! BLAS version 21 Jul 19
      endif
     case default
      FATAL( tr00ab, .true., invalid jderiv )
     end select

     lcpu = GETTIME ! record time taken in dgesvx; 09 Nov 17;

     select case( Lsvdiota )

     case( 0 ) ! Lsvdiota = 0; use linear solver to invert linear equations that define the straight fieldline angle; 01 Jul 14;

      if04aaf = 1

      select case( jderiv )

      case( 0 ) ! Lsvdiota = 0; jderiv = 0; 02 Sep 14;

       MM = 1

       call dgesvx( 'N', 'N', NN, MM, dmatrix(1:NN,1:NN,0), NN, FAA(1:NN,1:NN), NN, ipiv(1:NN),  &
                 equed, Rdgesvx(1:NN), Cdgesvx(1:NN), drhs(1:NN,0:0), NN, dlambda(1:NN,0:0),    &
		 NN, rcond, ferr, berr, work4(1:4*NN), iwork4(1:NN), idgesvx )

       ;                 ldiota(innout,    0) = dlambda(1,  0) ! return intent out; 21 Apr 13;
       ;                 dlambdaout(1:NN, lvol, innout) = dlambda(1:NN,0) 

      case( 1 ) ! Lsvdiota = 0; jderiv = 1; 02 Sep 14;

       MM = 2
       if( iflag.eq.-1 ) then ; drhs(1:NN, 1) = drhs(1:NN,-1)
        ;                     ; drhs(1:NN, 2) = zero
       endif

       call DCOPY(NN*NN, omatrix(1,1), 1, dmatrix(1,1,0), 1) ! BLAS version 21 Jul 19

       call dgesvx( 'N', 'N', NN, MM, dmatrix(1:NN,1:NN,0), NN, FAA(1:NN,1:NN), NN, ipiv(1:NN),    &
                   equed, Rdgesvx(1:NN), Cdgesvx(1:NN), drhs(1:NN,1:MM), NN, dlambda(1:NN,1:MM),    &
	           NN, rcond, ferr2(1:MM), berr2(1:MM), work4(1:4*NN), iwork4(1:NN), idgesvx )

       if( iflag.eq. 2 ) ldiota(innout, 1:2) = dlambda(1,1:2) ! return intent out; 21 Apr 13;
       if( iflag.eq.-1 ) ldiota(innout,-1  ) = dlambda(1,  1) ! return intent out; 21 Apr 13;

      case default

       FATAL( tr00ab, .true., invalid jderiv )

      end select ! end of select case jderiv; 02 Sep 14;

      cput = GETTIME

      select case( idgesvx )                                                                                           !12345678901234567
      case( 0   )    ; if( Wtr00ab ) write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput-lcpu, "solved Fourier ; ", dlambda(1,0)
      case( 1:  )    ;               write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput-lcpu, "singular ;       "
      case( :-1 )    ;               write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgesvx", idgesvx, cput-lcpu, "input error ;    "
      case default ;               FATAL( tr00ab, .true., illegal ifail returned by dgesvx )
      end select

      FATAL( tr00ab, idgesvx.ne.0, failed to construct straight-fieldline angle using dgesvx )

     case( 1 ) ! Lsvdiota = 1; use least-squares to invert linear equations that define the straight fieldline angle; 01 Jul 14;


      nlvl   = max(0, int(log( real(NN)/26 )/log(2.0D0))+1)
      Lwork  = (63+8*nlvl)*NN+676
      Liwork = max(1,11*NN+3*nlvl*NN)

      SALLOCATE( work, (1:Lwork), zero )
      if (allocated(iwork)) then
        DALLOCATE(iwork)
      endif
      SALLOCATE( iwork, (1:Liwork), zero )

      select case( jderiv )

      case( 0 ) ! Lsvdiota = 1; jderiv = 0; 02 Sep 14;

       dlambda(1:NN,0) = drhs(1:NN,0) ! on entry, rhs; on exit, solution; 20 Jun 14;

       call dgelsd( NN, NN, 1, dmatrix(1:NN,1:NN,0), NN, dlambda(1:NN,0), NN, sval(1:NN), rcond, Irank, &
                    work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd )

       ldiota(innout,0) = dlambda(1,0)
       dlambdaout(1:NN, lvol, innout) = dlambda(1:NN,0) 

      case( 1 ) ! Lsvdiota = 1; jderiv = 1; 02 Sep 14;

       if(     iflag.eq. 2 ) then
        do imupf = 1, 2
         dmatrix(1:NN,1:NN,0) = omatrix(1:NN,1:NN) ; dlambda(1:NN,imupf) = drhs(1:NN,imupf)

         call dgelsd( NN, NN, 1, dmatrix(1:NN,1:NN,0), NN, dlambda(1:NN,imupf), NN, sval(1:NN), rcond, Irank, &
                      work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd )

         ldiota(innout,imupf) = dlambda(1,imupf)
        enddo
       elseif( iflag.eq.-1 ) then
        do imupf = -1, -1
         dmatrix(1:NN,1:NN,0) = omatrix(1:NN,1:NN) ; dlambda(1:NN,imupf) = drhs(1:NN,imupf)

         call dgelsd( NN, NN, 1, dmatrix(1:NN,1:NN,0), NN, dlambda(1:NN,imupf), NN, sval(1:NN), rcond, Irank, &
                      work(1:Lwork), Lwork, iwork(1:Liwork), idgelsd )

         ldiota(innout,imupf) = dlambda(1,imupf)
        enddo
       else
        FATAL( tr00ab, .true., invalid iflag )
       endif

      case default

       FATAL( tr00ab, .true., invalid jderiv )

      end select ! end of select case( jderiv) ; 02 Sep 14;

      DALLOCATE(work)

      cput = GETTIME

      select case( idgelsd )                                                                                           !12345678901234567
      case( 0   )    ; if( Wtr00ab)  write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput-lcpu, "solved Fourier ; ", dlambda(1,0)
      case( :-1 )    ;               write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput-lcpu, "input error ;    "
      case( 1:  )    ;               write(ounit,1030) cput-cpus, myid, lvol, innout, id, "idgelsd", idgelsd, cput-lcpu, "QR failed ;      "
      case default ;               FATAL( tr00ab, .true., illegal ifail returned by f04arf )
      end select

      FATAL( tr00ab, idgelsd.ne.0, failed to construct straight-fieldline angle using dgelsd )

      dmatrix(1:NN,1:NN, 0) = omatrix(1:NN,1:NN) ! original "unperturbed" matrix; 30 Jan 13;

      case default

       FATAL( tr00ab, .true., illegal Lsvdiota )

      end select ! end of select case( Lsvdiota ) ; 02 Sep 14;

1030 format("tr00ab : ",f10.2," ; myid=",i3," ; lvol=",i3," ; innout="i2" ; jderiv="i2" ; "a7"="i2" ; time="f10.4" ; "a17,:" [d]iota="es17.09" ;")



    enddo ! end of do jderiv; 31 Jan 13;

   endif ! end of if( Lsparse.eq.0 .or. Lsparse.eq.3 );


   if( Lsparse.eq.0 .or. Lsparse.eq.3 ) then
     DALLOCATE( dmatrix )
     DALLOCATE( omatrix )
     DALLOCATE( FAA )
   endif




  enddo ! end of do innout; 29 Jan 13;







end subroutine tr00ab


