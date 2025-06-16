

subroutine volume( lvol, vflag )



  use constants, only : zero, half, one, two, four, third, quart, pi2

  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wvolume, Igeometry, Nvol, pscale

  use cputiming

  use allglobal, only : myid, cpus, MPI_COMM_SPEC, &
                        YESstellsym, Mvol, &
                        Ntz, mn, im, in, iRbc, iZbs, iRbs, iZbc, &
                        djkp, djkm, &
                        vvolume, dvolume, &
                        Rij, Zij, cosi, sini, &
                        dBdX, &
                        pi2nfp, pi2pi2nfp, pi2pi2nfpquart



  use mpi
  implicit none
  INTEGER   :: ierr, astat, ios, nthreads, ithread
  real(8)      :: cput, cpui, cpuo=0

  INTEGER, intent(in) :: lvol

  INTEGER             :: vflag, Lcurvature

  INTEGER             :: jvol, ii, jj, kk, mi, ni, mj, nj, mk, nk, innout

  real(8)                :: vol(0:1), vint(1:Ntz)

  real(8)                :: Rei, Roi, Zei, Zoi, Rej, Roj, Zej, Zoj, Rek, Rok, Zek, Zok

  real(8)                :: AA, BB, CC, DD, lss

  






  if( lvol.gt.Nvol ) then ; vvolume(lvol) = one ; dvolume = zero ; return ! this can only be the vacuum region; provide default value; 13 Sep 13;
  endif



  ;  vol(0:1) = zero ! initialize;
  ; dvolume   = zero ! derivatives of volume wrt interface inner/outer, R/Z harmonic;



  do innout = 0, 1 ! will subtract inner volume from outer volume to obtain enclosed volume; 26 Feb 13; ! this does seem a little wasteful; 13 Sep 13;

   jvol = lvol - 1 + innout ! labels inner or outer interface; 13 Sep 13;




   select case( Igeometry )



   case( 1 ) !> <li> \c Igeometry.eq.1 : Cartesian : \f$\sqrt g = R_s\f$


    vol(innout) = iRbc(1,jvol) ! 20 Jun 14;



    if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.1 ) then ! compute derivative of volume;
     if( dBdX%issym.eq.0 ) dvolume = one ! note that the sign factor for the lower interface is included below; 20 Jun 14;
    endif



   case( 2 ) !> <li> \c Igeometry.eq.2 : cylindrical : \f$\sqrt g = R R_s = \frac{1}{2}\partial_s (R^2)\f$




    if( YESstellsym ) then

     do ii = 1, mn ; mi = im(ii) ; ni = in(ii)

      do jj = 1, mn ; mj = im(jj) ; nj = in(jj)

       vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) * ( djkp(ii,jj) + djkm(ii,jj) )


       if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.ii ) then ! compute derivative of volume;
        dvolume = dvolume + iRbc(jj,jvol) * ( djkp(jj,ii) + djkm(jj,ii) + djkp(ii,jj) + djkm(ii,jj) )
       endif

      enddo ! end of do jj; 02 Sep 14;

     enddo ! end of do ii; 02 Sep 14;

    else ! NOTstellsym;

     do ii = 1, mn ; mi = im(ii) ; ni = in(ii)
      do jj = 1, mn ; mj = im(jj) ; nj = in(jj)

       vol(innout) = vol(innout) + iRbc(ii,jvol) * iRbc(jj,jvol) * ( djkp(ii,jj) + djkm(ii,jj) ) &
                                 + iRbs(ii,jvol) * iRbs(jj,jvol) * ( djkp(ii,jj) - djkm(ii,jj) )


       if( dBdX%L .and. dBdX%innout.eq.innout .and. dBdX%ii.eq.ii ) then ! compute derivative of volume;
        if( dBdX%issym.eq.0 ) then !     stellarator-symmetric harmonic; dV/dRei ; 13 Sep 13;
        dvolume = dvolume + iRbc(jj,jvol) * ( djkp(jj,ii) + djkm(jj,ii) + djkp(ii,jj) + djkm(ii,jj) )
        else
if( .true. ) then
     write(6,'("volume :      fatal : myid=",i3," ; .true. ; derivatives of volume under construction;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "volume : .true. : derivatives of volume under construction ;"
   endif
        dvolume = dvolume + iRbs(jj,jvol) * ( djkp(jj,ii) - djkm(jj,ii) + djkp(ii,jj) - djkm(ii,jj) ) ! needs to be checked; 02 Sep 14;
        endif
       endif

      enddo
     enddo

    endif ! end of if( YESstellsym ) ; 11 Aug 14;




   case( 3 ) !> <li> \c Igeometry.eq.3 : toroidal : \f${\bf x}\cdot {\bf e}_\theta \times {\bf e}_\zeta  = R ( Z R_\theta - R Z_\theta ) \f$


    if (lvol.eq.1 .and. innout.eq.0) then
      vol(1) = zero
      dvolume = zero
    else

      Lcurvature = 1

      lss = innout * two - one
      call coords( lvol, lss, Lcurvature, Ntz, mn  )

      vint = Rij(1:Ntz,0,0) * (Zij(1:Ntz,0,0)*Rij(1:Ntz,2,0) - Zij(1:Ntz,2,0)*Rij(1:Ntz,0,0))
      vol(innout) = four * sum(vint) / float(Ntz)

      if( dBdX%L .and. dBdX%innout.eq.innout  ) then ! compute derivative of volume;

        ii = dBdX%ii ! shorthand

        if( dBdX%irz.eq.0 ) then ! compute derivatives wrt R;

          if( dBdX%issym.eq.0 ) then
            vint = cosi(1:Ntz,ii) * (Zij(1:Ntz,0,0)*Rij(1:Ntz,2,0) - Zij(1:Ntz,2,0)*Rij(1:Ntz,0,0)) &
                + Rij(1:Ntz,0,0) * (-im(ii) * Zij(1:Ntz,0,0)*sini(1:Ntz,ii)) &
                + Rij(1:Ntz,0,0) * (-Zij(1:Ntz,2,0)*cosi(1:Ntz,ii))
          else
            vint = sini(1:Ntz,ii) * (Zij(1:Ntz,0,0)*Rij(1:Ntz,2,0) - Zij(1:Ntz,2,0)*Rij(1:Ntz,0,0)) &
                + Rij(1:Ntz,0,0) * (+im(ii) * Zij(1:Ntz,0,0)*cosi(1:Ntz,ii)) &
                + Rij(1:Ntz,0,0) * (-Zij(1:Ntz,2,0)*sini(1:Ntz,ii))
          endif

        else ! matches if( dBdX%irz.eq.0 ) then; compute derivative wrt Z;

          if( dBdX%issym.eq.0 ) then
            vint = Rij(1:Ntz,0,0) * (sini(1:Ntz,ii)*Rij(1:Ntz,2,0)) &
                + Rij(1:Ntz,0,0) * (-im(ii)*cosi(1:Ntz,ii)*Rij(1:Ntz,0,0))
          else
            vint = Rij(1:Ntz,0,0) * (cosi(1:Ntz,ii)*Rij(1:Ntz,2,0)) &
                + Rij(1:Ntz,0,0) * (+im(ii)*sini(1:Ntz,ii)*Rij(1:Ntz,0,0))
          endif

        endif ! end of if( dBdX%irz.eq.0 )

        dvolume = four * sum(vint) / float(Ntz)

      else

        dvolume = zero

      endif  ! end of if( dBdX%L .and. dBdX%innout.eq.innout )
    endif  ! lvol.eq.1 .and. innout.eq.0


   end select




  enddo ! end of innout loop; 26 Feb 13;




  select case( Igeometry)
  case( 1 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfp              ; dvolume = dvolume * pi2pi2nfp                ! 20 Jun 14;
  case( 2 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfpquart         ; dvolume = dvolume * pi2pi2nfpquart
  case( 3 ) ; vvolume(lvol) = ( vol(1) - vol(0) ) * pi2pi2nfpquart * third ; dvolume = dvolume * pi2pi2nfpquart * third
  case( 4 ) ; vvolume(lvol) = one                                          ; dvolume = zero ! this is under construction; 04 Dec 14;
if( abs(pscale).gt.vsmall ) then
     write(6,'("volume :      fatal : myid=",i3," ; abs(pscale).gt.vsmall ; need to compute volume;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "volume : abs(pscale).gt.vsmall : need to compute volume ;"
   endif
  case default
if( .true. ) then
     write(6,'("volume :      fatal : myid=",i3," ; .true. ; invalid Igeometry;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "volume : .true. : invalid Igeometry ;"
   endif
  end select

  if( dBdX%innout.eq.0 ) dvolume = - dvolume



  if( Wvolume ) then
   cput = MPI_WTIME()
   write(ounit,'("volume : ",f10.2," : myid=",i3," ; Igeometry=",i2," ; vvolume(",i3," ) =",es23.15" ;")') cput-cpus, myid, Igeometry, lvol, vvolume(lvol)
  endif



if( vflag.eq.0 .and. vvolume(lvol).lt.small ) then
     write(6,'("volume :      fatal : myid=",i3," ; vflag.eq.0 .and. vvolume(lvol).lt.small ; volume cannot be zero or negative;")') myid
     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )
     stop "volume : vflag.eq.0 .and. vvolume(lvol).lt.small : volume cannot be zero or negative ;"
   endif

  if( vvolume(lvol).lt.small ) then
   write(ounit,'("volume : ", 10x ," : myid=",i3," ; lvol=",i3," ; vvolume=",es13.5," ; volume cannot be zero or negative ;")') myid, lvol, vvolume(lvol)
   vvolume(lvol) = +9.9E+09
   vflag = 1
  else
   vflag = 0
  endif







end subroutine volume



