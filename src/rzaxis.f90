

subroutine rzaxis( Mvol, mn, inRbc, inZbs, inRbs, inZbc, ivol, LcomputeDerivatives )



  use constants, only : zero, one, half, two

  use numerical, only : vsmall

  use fileunits, only : ounit

  use inputlist, only : Wrzaxis, Igeometry, Ntor, Lcheck, Wmacros, Lreflect, Ntoraxis, Lrzaxis

  use cputiming, only : Trzaxis

  use allglobal, only : ncpu, myid, cpus, im, in, MPI_COMM_SPEC, &
                        ajk, Nt, Nz, Ntz, &
                        Rij, Zij, sg, cosi, sini, &
                        ijreal, ijimag, jireal, jiimag, jkreal, jkimag, kjreal, kjimag, &
                        efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, cosi, sini, &
                        YESstellsym, NOTstellsym, Lcoordinatesingularity, &
                        dRodR, dRodZ, dZodR, dZodZ, &
                        dRadR, dRadZ, dZadR, dZadZ, &
                        iRbc, iZbs, iRbs, iZbc, &
                        dBdX



  LOCALS

  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;

  INTEGER, intent(in)    :: Mvol, mn, ivol
  REAL                   :: inRbc(1:mn,0:Mvol), inZbs(1:mn,0:Mvol), inRbs(1:mn,0:Mvol), inZbc(1:mn,0:Mvol)
  REAL                   :: jRbc(1:mn,0:Mvol), jZbs(1:mn,0:Mvol), jRbs(1:mn,0:Mvol), jZbc(1:mn,0:Mvol)
  REAL                   :: tmpRbc(1:mn,0:Mvol), tmpZbs(1:mn,0:Mvol), tmpRbs(1:mn,0:Mvol), tmpZbc(1:mn,0:Mvol) ! use as temp matrices to store iRbc etc

  REAL                   :: jacbase(1:Ntz), jacbasec(1:mn), jacbases(1:mn) ! the 2D Jacobian and its Fourier
  REAL                   :: junkc(1:mn), junks(1:mn) ! these are junk matrices used for fft

  INTEGER                :: jvol, ii, ifail, jj, id, issym, irz, imn
  INTEGER                :: idJc, idJs, idRc, idRs, idZc, idZs

  INTEGER                :: Lcurvature

  INTEGER                :: Njac, idgetrf, idgetrs ! internal variables used in Jacobian method
  REAL, allocatable      :: jacrhs(:), djacrhs(:), jacmat(:,:), djacmat(:,:), solution(:), LU(:,:) ! internal matrices used in Jacobian method
  INTEGER, allocatable   :: ipiv(:)   ! internal matrices used in  Jacobian method



  

  write(*,*) "Calling rzaxis"


  jvol = 0 ! this identifies the "surface" in which the poloidal averaged harmonics will be placed; 19 Jul 16;

  Ntoraxis = min(Ntor,Ntoraxis)



  select case( Igeometry )

  case( 1:2 )

   inRbc(1:mn,jvol) = zero
   inRbs(1:mn,jvol) = zero

   if ( Igeometry.eq.1 .and. Lreflect.eq.1) then ! reflect upper and lower bound in slab, each take half the amplitude
    inRbc(2:mn,0) = -inRbc(2:mn,Mvol)
   if( NOTstellsym ) then
    inRbs(2:mn,0) = -inRbs(2:mn,Mvol)
    endif
   endif

  case(   3 )

   if (Lrzaxis .eq. 1) then ! use centroid method

    call invfft( mn, im(1:mn), in(1:mn), im(1:mn) * inRbs(1:mn,ivol), - im(1:mn) * inRbc(1:mn,ivol), &
                                          im(1:mn) * inZbs(1:mn,ivol), - im(1:mn) * inZbc(1:mn,ivol), &
                  Nt, Nz, jkreal(1:Ntz), jkimag(1:Ntz) ) ! R_\t, Z_\t; 03 Nov 16;

    ijreal(1:Ntz) = sqrt( jkreal(1:Ntz)**2 + jkimag(1:Ntz)**2 ) ! dl ; 11 Aug 14;
    ijimag(1:Ntz) = zero

    jireal(1:Ntz) = ijreal(1:Ntz) ! dl ; 19 Sep 16;

    ifail = 0
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail ) ! Fourier harmonics of differential poloidal length; 11 Mar 16;

    efmn(1:mn) = efmn(1:mn) * ajk(1:mn) ! poloidal integration of length; only take m=0 harmonics; 11 Aug 14;
    ofmn(1:mn) = ofmn(1:mn) * ajk(1:mn)
    cfmn(1:mn) = zero
    sfmn(1:mn) = zero

    call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), & ! map length = "integrated dl" back to real space; 19 Sep 16;
                  Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz) )

    jiimag(1:Ntz) = ijreal(1:Ntz) !  L ; 19 Sep 16;


    call invfft( mn, im(1:mn), in(1:mn),            inRbc(1:mn,ivol),              inRbs(1:mn,ivol), &
                                                    inZbc(1:mn,ivol),              inZbs(1:mn,ivol), &
                  Nt, Nz, kjreal(1:Ntz), kjimag(1:Ntz) ) ! R, Z; 03 Nov 16;

    ijreal(1:Ntz) = kjreal(1:Ntz) * jireal(1:Ntz) ! R dl;
    ijimag(1:Ntz) = kjimag(1:Ntz) * jireal(1:Ntz) ! Z dl;

    ifail = 0
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! Fourier harmonics of weighted R & Z; 11 Mar 16;

    evmn(1:mn) = evmn(1:mn) * ajk(1:mn) ! poloidal integration of R dl; 19 Sep 16;
    odmn(1:mn) = odmn(1:mn) * ajk(1:mn)
    comn(1:mn) = comn(1:mn) * ajk(1:mn) ! poloidal integration of Z dl; 19 Sep 16;
    simn(1:mn) = simn(1:mn) * ajk(1:mn)

    call invfft( mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), &
                  Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz) )

    ijreal(1:Ntz) = ijreal(1:Ntz) / jiimag(1:Ntz) ! Ro; 19 Sep 16;
    ijimag(1:Ntz) = ijimag(1:Ntz) / jiimag(1:Ntz) ! Zo; 19 Sep 16;

    kjreal(1:Ntz) = kjreal(1:Ntz) - ijreal(1:Ntz) ! \Delta R = R_1 - R_0 ; 03 Nov 16;
    kjimag(1:Ntz) = kjimag(1:Ntz) - ijimag(1:Ntz) ! \Delta R = Z_1 - Z_0 ; 03 Nov 16;

    ifail = 0
    call tfft( Nt, Nz, ijreal(1:Ntz), ijimag(1:Ntz), &
                mn, im(1:mn), in(1:mn), inRbc(1:mn,jvol), inRbs(1:mn,jvol), inZbc(1:mn,jvol), inZbs(1:mn,jvol), ifail )

   else if ( Lrzaxis .eq. 2) then ! use Jacobian m=1 harmonic elimination method

    tmpRbc = iRbc
    tmpZbs = iZbs
    tmpRbs = iRbs
    tmpZbc = iZbc

    jRbc = inRbc
    jZbs = inZbs
    jRbs = inRbs
    jZbc = inZbc

    if( YESstellsym ) then
      Njac = 2 * Ntoraxis + 1
    else
      Njac = 2 * (2 * Ntoraxis + 1)
    end if

    SALLOCATE( jacrhs, (1:Njac), zero )
    SALLOCATE( jacmat, (1:Njac, 1:Njac), zero )
    SALLOCATE( LU, (1:Njac, 1:Njac), zero )
    SALLOCATE( solution, (1:Njac), zero )
    SALLOCATE( ipiv, (1:Njac), 0)

    iRbc(1:mn,1) = jRbc(1:mn, ivol)
    iZbs(1:mn,1) = jZbs(1:mn, ivol)
    iRbs(1:mn,1) = jRbs(1:mn, ivol)
    iZbc(1:mn,1) = jZbc(1:mn, ivol)

    iRbc(1:mn,0) = zero
    iZbs(1:mn,0) = zero
    iRbs(1:mn,0) = zero
    iZbc(1:mn,0) = zero

    iRbc(1:Ntor+1,0) = jRbc(1:Ntor+1, ivol)
    iZbs(1:Ntor+1,0) = jZbs(1:Ntor+1, ivol)
    iRbs(1:Ntor+1,0) = jRbs(1:Ntor+1, ivol)
    iZbc(1:Ntor+1,0) = jZbc(1:Ntor+1, ivol)

    Lcoordinatesingularity = .true.
    Lcurvature = 1
    dBdX%innout = 1

    idJc = Ntoraxis+1                 ! rhs index for J cos n=0 term
    idJs = Ntoraxis+1 + 2*Ntoraxis+1  ! rhs index for J sin n=0 term
    idRc = 1
    idZs = Ntoraxis + 1
    idRs = 2 * Ntoraxis + 1
    idZc = 3 * Ntoraxis + 2

    WCALL( rzaxis, coords, (1, one, Lcurvature, Ntz, mn ))

    jacbase = sg(1:Ntz,0) / Rij(1:Ntz,0,0)  ! extract the baseline 2D jacobian, note the definition here does not have the R factor

    call tfft( Nt, Nz, jacbase, Rij, &
               mn, im(1:mn), in(1:mn), jacbasec(1:mn), jacbases(1:mn), junkc(1:mn), junks(1:mn), ifail )

    if (YESstellsym) then
      jacrhs = -jacbasec(2*(Ntor+1)-Ntoraxis:2*(Ntor+1)+Ntoraxis)
    else
      jacrhs(1:2*Ntoraxis+1) = -jacbasec(2*(Ntor+1)-Ntoraxis:2*(Ntor+1)+Ntoraxis)
      jacrhs(2*Ntoraxis+2:Njac) = -jacbases(2*(Ntor+1)-Ntoraxis:2*(Ntor+1)+Ntoraxis)
    end if !if (YESstellsym)

    if (YESstellsym) then

      do ii = -Ntoraxis, Ntoraxis
        do jj = 1, Ntoraxis

          if (ii-jj .ge. -Ntor) then
            id = 2 * (Ntor + 1) + ii - jj
            jacmat(ii+Ntoraxis+1, jj+1) = jacmat(ii+Ntoraxis+1, jj+1) - jZbs(id,ivol)
            jacmat(ii+Ntoraxis+1, Ntoraxis+1+jj) = jacmat(ii+Ntoraxis+1, Ntoraxis+1+jj) + jRbc(id,ivol)
          end if ! if (ii-jj .ge. -Ntor)

          if (ii+jj .le. Ntor) then
            id = 2 * (Ntor + 1) + ii + jj
            jacmat(ii+Ntoraxis+1, jj+1) = jacmat(ii+Ntoraxis+1, jj+1) - jZbs(id,ivol)
            jacmat(ii+Ntoraxis+1, Ntoraxis+1+jj) = jacmat(ii+Ntoraxis+1, Ntoraxis+1+jj) - jRbc(id,ivol)
          end if ! if (ii+jj .le. Ntor)

        end do ! jj

        id = 2 * (Ntor + 1) + ii
        jacmat(ii+Ntoraxis+1, 1) = - two * jZbs(id,ivol)

      end do ! ii

    else ! for NOTstellsym

      do ii = -Ntoraxis, Ntoraxis
        do jj = 1, Ntoraxis

          if (ii-jj .ge. -Ntor) then
            id = 2 * (Ntor + 1) + ii - jj
            jacmat(ii+idJc, jj+idRc) = jacmat(ii+idJc, jj+idRc) - jZbs(id,ivol)
            jacmat(ii+idJc, jj+idZs) = jacmat(ii+idJc, jj+idZs) + jRbc(id,ivol)
            jacmat(ii+idJc, jj+idRs) = jacmat(ii+idJc, jj+idRs) - jZbc(id,ivol)
            jacmat(ii+idJc, jj+idZc) = jacmat(ii+idJc, jj+idZc) + jRbs(id,ivol)

            jacmat(ii+idJs, jj+idRc) = jacmat(ii+idJs, jj+idRc) + jZbc(id,ivol)
            jacmat(ii+idJs, jj+idZs) = jacmat(ii+idJs, jj+idZs) + jRbs(id,ivol)
            jacmat(ii+idJs, jj+idRs) = jacmat(ii+idJs, jj+idRs) - jZbs(id,ivol)
            jacmat(ii+idJs, jj+idZc) = jacmat(ii+idJs, jj+idZc) - jRbc(id,ivol)

          end if ! if (ii-jj .ge. -Ntor)

          if (ii+jj .le. Ntor) then
            id = 2 * (Ntor + 1) + ii + jj

            jacmat(ii+idJc, jj+idRc) = jacmat(ii+idJc, jj+idRc) - jZbs(id,ivol)
            jacmat(ii+idJc, jj+idZs) = jacmat(ii+idJc, jj+idZs) - jRbc(id,ivol)
            jacmat(ii+idJc, jj+idRs) = jacmat(ii+idJc, jj+idRs) + jZbc(id,ivol)
            jacmat(ii+idJc, jj+idZc) = jacmat(ii+idJc, jj+idZc) + jRbs(id,ivol)

            jacmat(ii+idJs, jj+idRc) = jacmat(ii+idJs, jj+idRc) + jZbc(id,ivol)
            jacmat(ii+idJs, jj+idZs) = jacmat(ii+idJs, jj+idZs) - jRbs(id,ivol)
            jacmat(ii+idJs, jj+idRs) = jacmat(ii+idJs, jj+idRs) + jZbs(id,ivol)
            jacmat(ii+idJs, jj+idZc) = jacmat(ii+idJs, jj+idZc) - jRbc(id,ivol)
          end if ! if (ii+jj .le. Ntor)

        end do ! jj

        id = 2 * (Ntor + 1) + ii
        jacmat(ii+idJc, idRc) = - two * jZbs(id,ivol)
        jacmat(ii+idJc, idZc) = + two * jRbs(id,ivol)
        jacmat(ii+idJs, idRc) = + two * jZbc(id,ivol)
        jacmat(ii+idJs, idZc) = - two * jRbc(id,ivol)

      end do ! ii

    endif ! if (YESstellsym)

    jacmat = jacmat * half ! because we are using (1+s)/2 instead of s

    LU = jacmat
    call DGETRF( Njac, Njac, LU, Njac, ipiv, idgetrf ) ! LU factorization
    solution = jacrhs
    call DGETRS('N', Njac, 1, LU, Njac, ipiv, solution, Njac, idgetrs ) ! sovle linear equation

    if( idgetrf .lt. 0 .or. idgetrs .lt. 0 ) then
    ;             write(ounit,1010) cput-cpus, myid, ivol, idgetrf, idgetrs, "input error ;     "
    elseif( idgetrf .gt. 0 ) then
    ;             write(ounit,1010) cput-cpus, myid, ivol, idgetrf, idgetrs,  "singular ;        "
    endif

1010 format("rzaxis : ",f10.2," : myid=",i3," ; ivol=",i3," idgetrf idgetrs=",i3,' ',i3," ; "a34)

    iRbc = tmpRbc
    iZbs = tmpZbs
    iRbs = tmpRbs
    iZbc = tmpZbc


    inRbc(:,jvol) = zero
    inZbs(:,jvol) = zero
    inRbs(:,jvol) = zero
    inZbc(:,jvol) = zero

    inRbc(1:Ntoraxis+1,jvol) = inRbc(1:Ntoraxis+1,ivol) - solution(idRc:idRc+Ntoraxis)
    inZbs(2:Ntoraxis+1 ,jvol) = inZbs(2:Ntoraxis+1,ivol) - solution(idZs+1:idZs+Ntoraxis)
    if (YESstellsym) then
      inRbs(1:Ntoraxis+1,jvol) = zero
      inZbc(2:Ntoraxis+1,jvol) = zero
    else
      inRbs(2:Ntoraxis+1,jvol) = inRbs(2:Ntoraxis+1,ivol) - solution(idRs+1:idRs+Ntoraxis)
      inZbc(1:Ntoraxis+1,jvol) = inZbc(1:Ntoraxis+1,ivol) - solution(idZc:idZc+Ntoraxis)
    endif ! YESstellsym

    DALLOCATE( jacrhs )
    DALLOCATE( jacmat )
    DALLOCATE( LU )
    DALLOCATE( solution )
    DALLOCATE( ipiv )



   end if ! end of forking based on Lrzaxis ; 10 Jan 20



  end select ! end of select case( Igeometry ) ; 08 Feb 16;





end subroutine rzaxis

