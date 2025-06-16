
subroutine casing( teta, zeta, gBn, icasing )



  use constants, only : zero, pi, pi2

  use numerical, only :

  use fileunits, only : ounit, vunit

  use inputlist, only : Wmacros, Wcasing, vcasingtol, vcasingits, vcasingper

  use cputiming, only : Tcasing

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, globaljk, Dxyz, Nxyz, ext



  LOCALS

  REAL, intent(in)     :: teta, zeta ! arbitrary location; Cartesian;
  REAL, intent(out)    :: gBn ! magnetic field; Cartesian;
  INTEGER, intent(out) :: icasing

  INTEGER, parameter   :: Ndim = 2, Nfun = 1

  INTEGER              :: ldim, lfun, minpts, maxpts, Lrwk, idcuhre, jk, irestart, funcls, key, num, maxsub
  REAL                 :: integrals(1:Nfun), low(1:Ndim), upp(1:Ndim), labs, lrel, absest(1:Nfun)
  REAL, allocatable    :: rwk(:)

  external             :: dvcfield

  BEGIN(casing)



  jk = globaljk ! shorthand; globaljk is a "global" variable which must be passed through to subroutine dvcfield;



  gBn = zero ! initialize intent(out);



  ldim = Ndim ! number of dimensions of the integral; surface integral => number of dimensions = 2;

  low(1) = teta - vcasingper * pi ; upp(1) = teta + vcasingper * pi
  low(2) = zeta - vcasingper * pi ; upp(2) = zeta + vcasingper * pi

  key = 0

  minpts = vcasingits ! minimum number of function evaluations; provided on input;

  num = 65 ! see documentation for dcuhre;


  maxpts = 16777216

  lfun = Nfun ! number of functions to be integrated; require three components of magnetic field, Bx, By and Bz; and their derivatives wrt x,y,z;

  labs = vcasingtol ; lrel = vcasingtol ! absolute and relative accuracy requested; vcasingtol is an input parameter;

  maxsub = ( maxpts - num ) / ( 2 * num ) + 1

  Lrwk = maxsub * ( 2 * Ndim + 2 * Nfun + 2 ) + 17 * Nfun + 1

  SALLOCATE( rwk, (1:Lrwk), zero )

  irestart = 0 ; funcls = 0

  do ! will continually call until satisfactory accuracy has been achieved;

   idcuhre = 1

   call DCUHRE( ldim, lfun, low(1:Ndim), upp(1:Ndim), minpts, maxpts, dvcfield, labs, lrel, key, &
                Lrwk, irestart, integrals(1:lfun), absest(1:lfun), funcls, idcuhre, rwk(1:Lrwk) )

   gBn = integrals(1)

   cput = GETTIME
   select case( idcuhre ) !                                                                                      "123456789012345678901234"
   case(0)      ;
    ;           ; exit
   case(1)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "maxpts too smal;        "
    ;           ;!exit
   case(2)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "illegal key;            "
    ;           ; exit
   case(3)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "illegal Ndim;           "
    ;           ; exit
   case(4)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "key.eq.1 & Ndim.ne.2;   "
    ;           ; exit
   case(5)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "key.eq.2 & Ndim.ne.3;   "
    ;           ; exit
   case(6)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "numfun < 1;             "
    ;           ; exit
   case(7)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "volume is zero;         "
    ;           ; exit
   case(8)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "maxpts < 3*NUM;         "
    ;           ; exit
   case(9)      ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "maxpts < minpts;        "
    ;           ; exit
   case(10)     ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "epsabs < 0 & epsrel < 0;"
    ;           ; exit
   case(11)     ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "NW is too small;        "
    ;           ; exit
   case(12)     ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "illegal irestart;       "
    ;           ; exit
   case default ; write(ounit,1001) cput-cpus, myid, Dxyz(1:3,jk), gBn, absest(1:Nfun), idcuhre, minpts, maxpts, "tryin' to kill me?      "
    ;           ; exit
   end select

   maxpts = 2 * maxpts ; minpts = funcls ; irestart = 0

   maxsub = ( maxpts - num ) / ( 2 * num ) + 1

   Lrwk = maxsub * ( 2 * Ndim + 2 * Nfun + 2 ) + 17 * Nfun + 1

   DALLOCATE(rwk)

   SALLOCATE(rwk, (1:Lrwk), zero)


  enddo ! end of virtual casing accuracy infinite-do-loop; 10 Apr 13;


1001 format("casing : ",f10.2," : myid=",i3," ; [x,y,z]=["es10.2" ,"es10.2" ,"es10.2" ]; gBn="es12.4" , ",&
            "err="es8.0" ; ifail="i3" ; min/max calls="2i12" ; "a24)



  icasing = idcuhre ! this is an error flag returned by casing;



  DALLOCATE(rwk)



  RETURN(casing)



end subroutine casing



subroutine dvcfield( Ndim, tz, Nfun, vcintegrand ) ! differential virtual-casing field; format is fixed by NAG requirements;



  use constants, only : zero, half, one, three, four

  use numerical, only : small

  use fileunits, only : ounit, vunit

  use inputlist, only : Wcasing, Nvol, Igeometry, Lrad, vcasingeps

  use cputiming, only :

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, &
                        pi2nfp, &
                        Mvol, &
                        mn, im, in, &
                        iRbc, iZbs, iRbs, iZbc, &
                        Ate, Aze, Ato, Azo, &
                        TT, &
                        YESstellsym, NOTstellsym, &
                        globaljk, Dxyz, Nxyz, &
                        first_free_bound



  LOCALS

  INTEGER , intent(in)  :: Ndim, Nfun
  REAL    , intent(in)  :: tz(1:Ndim)
  REAL    , intent(out) :: vcintegrand(1:Nfun) ! integrand; components of magnetic field due to plasma currents in Cartesian coordinates;

  INTEGER               :: ii, mi, ni, ll, ideriv, jk
  REAL                  :: dR(0:3), dZ(0:3), gBut, gBuz, gtt, gtz, gzz, sqrtg, Blt, Blz, czeta, szeta, arg, carg, sarg, XX, YY, ZZ, teta, zeta
  REAL                  :: jj(1:3), rr(1:3), distance(1:3), firstorderfactor

  REAL                  :: XXt, XXz, YYt, YYz, ZZt, ZZz, ds, Bxyz(1:3)


  dR(0:3) = zero ; ideriv = 0 ; gBut = zero ; gBuz = zero ! initialize summation of coordinates and tangential field;
  dZ(0:3) = zero

  teta = tz(1) ; zeta = tz(2) ! shorthand; 09 Mar 17;

  jk = globaljk



  select case( Igeometry )

  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;

   if( YESstellsym ) then

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;

     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)

     dR(0) = dR(0) +          (                   iRbc(ii,Nvol) ) * carg
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg                                            ) * half
     dR(2) = dR(2) + mi * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )
     dR(3) = dR(3) - ni * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )


     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1) ! contravariant; Jacobian comes later;
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1)
     enddo

    enddo ! end of do ii = 1, mn ;

   else ! NOTstellsym ; 08 Feb 16;

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;

     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)

     dR(0) = dR(0) +                              iRbc(ii,Nvol)   * carg +                   iRbs(ii,Nvol)   * sarg
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg + ( iRbs(ii,Mvol) - iRbs(ii,Nvol) ) * sarg ) * half
     dR(2) = dR(2) + mi * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )
     dR(3) = dR(3) - ni * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )


     do ll = 0, Lrad(Mvol)
      gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg + Azo(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1) ! contravariant; Jacobian comes later;
      gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg + Ato(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1)
     enddo

    enddo ! end of do ii = 1, mn ;

   endif ! end of if( YESstellsym ) ; 08 Feb 16;

  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;

   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )

  case( 3 ) ! Igeometry = 3 ; 09 Mar 17;

   if( YESstellsym ) then

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;

     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     dR(0) = dR(0) +          (                   iRbc(ii,Nvol) ) * carg
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg                                            ) * half
     dR(2) = dR(2) + mi * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )
     dR(3) = dR(3) - ni * ( - (                   iRbc(ii,Nvol) ) * sarg                                            )

     dZ(0) = dZ(0) +                                                       (                 iZbs(ii,Nvol) ) * sarg
     dZ(1) = dZ(1) +      (                                                ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
     dZ(2) = dZ(2) + mi * (                                                (                 iZbs(ii,Nvol) ) * carg )
     dZ(3) = dZ(3) - ni * (                                                (                 iZbs(ii,Nvol) ) * carg )

     if (first_free_bound) then
        do ll = 0, Lrad(Nvol)  ! 1 is for outside thr volume
           gBut = gBut - ( Aze(Nvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,1,1) ! contravariant; Jacobian comes later;
           gBuz = gBuz + ( Ate(Nvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,1,1)
        enddo
     else
        do ll = 0, Lrad(Mvol)
           gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1) ! contravariant; Jacobian comes later;
           gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg                                    ) * TT(ll,0,1)
        enddo
     endif

    enddo ! end of do ii = 1, mn ;

   else ! NOTstellsym ; 08 Feb 16;

    do ii = 1, mn ; mi = im(ii) ; ni = in(ii) ! loop over Fourier modes; construct surface current; slow transform required as position is arbitrary;

     arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)

     dR(0) = dR(0) +                              iRbc(ii,Nvol)   * carg +                   iRbs(ii,Nvol)   * sarg
     dR(1) = dR(1) +      (   (   iRbc(ii,Mvol) - iRbc(ii,Nvol) ) * carg + ( iRbs(ii,Mvol) - iRbs(ii,Nvol) ) * sarg ) * half
     dR(2) = dR(2) + mi * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )
     dR(3) = dR(3) - ni * ( -                     iRbc(ii,Nvol)   * sarg +                   iRbs(ii,Nvol)   * carg )

     dZ(0) = dZ(0) +                              iZbc(ii,Nvol)   * carg +                   iZbs(ii,Nvol)   * sarg
     dZ(1) = dZ(1) +      (   (   iZbc(ii,Mvol) - iZbc(ii,Nvol) ) * carg + ( iZbs(ii,Mvol) - iZbs(ii,Nvol) ) * sarg ) * half
     dZ(2) = dZ(2) + mi * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )
     dZ(3) = dZ(3) - ni * ( -                     iZbc(ii,Nvol)   * sarg +                   iZbs(ii,Nvol)   * carg )

     if (first_free_bound) then
        do ll = 0, Lrad(Mvol) ! omit the possible current sheet due to a jump in tangential field at the plasma boundary; Zhu 20190603;
           gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg + Azo(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,1,1) ! contravariant; Jacobian comes later;
           gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg + Ato(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,1,1)
        enddo
     else
        do ll = 0, Lrad(Mvol)
           gBut = gBut - ( Aze(Mvol,ideriv,ii)%s(ll) * carg + Azo(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1) ! contravariant; Jacobian comes later;
           gBuz = gBuz + ( Ate(Mvol,ideriv,ii)%s(ll) * carg + Ato(Mvol,ideriv,ii)%s(ll) * sarg ) * TT(ll,0,1)
        enddo
     endif

    enddo ! end of do ii = 1, mn ;

   endif ! end of if( YESstellsym ) ; 08 Feb 16;

  end select ! end of select case( Igeometry ) ; 09 Mar 17;



  select case( Igeometry )

  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;

   gtt = one + dR(2)*dR(2)
   gtz =       dR(2)*dR(3)
   gzz = one + dR(3)*dR(3)

   sqrtg = dR(1)

  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;

   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )

  case( 3 ) ! Igeometry = 3 ; 09 Mar 17;

   gtt = dR(2)*dR(2) + dZ(2)*dZ(2)
   gtz = dR(2)*dR(3) + dZ(2)*dZ(3)
   gzz = dR(3)*dR(3) + dZ(3)*dZ(3) + dR(0)*dR(0)

   sqrtg = dR(0) * ( dZ(1) * dR(2) - dR(1) * dZ(2) )

  end select



  Blt = ( gBut * gtt + gBuz * gtz ) / sqrtg
  Blz = ( gBut * gtz + gBuz * gzz ) / sqrtg



  select case( Igeometry )

  case( 1 ) ! Igeometry = 1 ; 09 Mar 17;

   XX =          teta ; XXt =           one ; XXz =          zero
   YY =          zeta ; YYt =          zero ; YYz =           one
   ZZ = dR(0)         ; ZZt = dR(2)         ; ZZz = dR(3)

  case( 2 ) ! Igeometry = 2 ; 09 Mar 17;

   FATAL( casing, .true., virtual casing under construction for cylindrical geometry )

  case( 3 ) ! Igeometry = 3 ; toroidal geometry;

   czeta = cos( zeta ) ; szeta = sin( zeta )

   XX = dR(0) * czeta ; XXt = dR(2) * czeta ; XXz = dR(3) * czeta - dR(0) * szeta ! 10 Apr 13;
   YY = dR(0) * szeta ; YYt = dR(2) * szeta ; YYz = dR(3) * szeta + dR(0) * czeta
   ZZ = dZ(0)         ; ZZt = dZ(2)         ; ZZz = dZ(3)

  end select



  rr(1:3) = (/ Dxyz(1,jk) - XX, &
               Dxyz(2,jk) - YY, &
               Dxyz(3,jk) - ZZ /)

  jj(1:3) = (/ Blz * XXt - Blt * XXz, &
               Blz * YYt - Blt * YYz, &
               Blz * ZZt - Blt * ZZz /)



  distance(2) = sum( rr(1:3) * rr(1:3) ) + vcasingeps**2 ! 04 May 17;

  distance(1) = sqrt( distance(2) ) ; distance(3) = distance(1) * distance(2) ! powers of distance; 24 Nov 16;

  firstorderfactor = ( one + three * vcasingeps**2 / distance(2) ) / distance(3) ! 04 May 17;



  Bxyz(1:3) = (/ jj(2) * rr(3) - jj(3) * rr(2), &
                 jj(3) * rr(1) - jj(1) * rr(3), &
                 jj(1) * rr(2) - jj(2) * rr(1)  /)

  vcintegrand(1) = sum( Bxyz(1:3) * Nxyz(1:3,jk) ) * firstorderfactor






  return



end subroutine dvcfield


