
subroutine mp00ac( Ndof, Xdof, Fdof, Ddof, Ldfjac, iflag ) ! argument list is fixed by NAG; ma02aa calls mp00ac through C05PCF;




  use constants, only : zero, half, one

  use numerical, only : small, machprec

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wmp00ac, Wtr00ab, Wcurent, Wma02aa, &
                        mu, helicity, iota, oita, curtor, curpol, Lrad, Ntor,&
                        Lconstraint, mupftol, &
                        Lmatsolver, NiterGMRES, epsGMRES, LGMRESprec, epsILU

  use cputiming, only : Tmp00ac

  use allglobal, only : myid, ncpu, cpus, ivol, MPI_COMM_SPEC, &
                        YESstellsym, NOTstellsym, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        mn, im, in, mns, &
                        Nt, Nz, & ! only required to pass through as arguments to tr00ab;
                        NAdof, &
                        dMA, dMB,      dMD,           dMG, &
                        Adotx, Ddotx,&
                        NdMASmax, NdMAS, dMAS, dMDS, idMAS, jdMAS, & ! preconditioning matrix
                        solution, GMRESlastsolution, &
                        dtflux, dpflux, &
                        diotadxup, dItGpdxtp, &
                        lBBintegral, lABintegral, &
                        xoffset, &
                        ImagneticOK, &
                        Ate, Aze, Ato, Azo, Mvol, Iquad, &
                        LILUprecond, GMRESlastsolution, ext



  LOCALS

  INTEGER, intent(in)  :: Ndof, Ldfjac
  REAL   , intent(in)  :: Xdof(1:Ndof)
  REAL                 :: Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof)
  INTEGER              :: iflag


  INTEGER, parameter   :: NB = 4 ! optimal workspace block size for LAPACK:DGECON;

  INTEGER              :: lvol, NN, MM, ideriv, lmns, ii, jj, nnz, Lwork

  INTEGER              :: idgetrf(0:1), idgetrs(0:1), idgerfs(0:1), idgecon(0:1)

  REAL                 :: lmu, dpf, dtf, dpsi(1:2), tpsi(1:2), ppsi(1:2), lcpu, test(2,2)

  REAL                 :: anorm, rcond, ferr(2), berr(2), signfactor

  CHARACTER            :: packorunpack

  INTEGER, allocatable :: ipiv(:), Iwork(:)

  REAL   , allocatable :: matrix(:,:), rhs(:,:), LU(:,:)

  REAL   , allocatable :: RW(:), RD(:,:)

  REAL   , allocatable :: matrixC(:,:)

  INTEGER, parameter   :: nrestart = 5 ! do GMRES restart after nrestart iterations
  INTEGER              :: maxfil   ! bandwidth for ILU subroutines, will be estimated

  INTEGER              :: NS, itercount, Nbilut

  REAL   , allocatable :: matrixS(:), bilut(:)
  INTEGER, allocatable :: ibilut(:),  jbilut(:)

  INTEGER, parameter   :: ipar_SIZE = 128
  INTEGER              :: ipar(ipar_SIZE), iluierr, RCI_REQUEST, nw, t1, t2, t3
  REAL                 :: fpar(ipar_SIZE), v1
  REAL, allocatable    :: wk(:)
  INTEGER,allocatable  :: jw(:), iperm(:)

  BEGIN(mp00ac)



  lvol = ivol ! recall that ivol is global;




  if( Lplasmaregion ) then

   ;                    ; lmu  = Xdof(1) - xoffset
   ;                    ; dtf  = dtflux(lvol)
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif

  else ! Lvacuumregion;


   ;                    ; lmu  = zero               ! restrict attention to strict vacuum field;

   ;                    ; dtf  = Xdof(1) - xoffset
   if( Ndof.eq.2 ) then ; dpf  = Xdof(2) - xoffset
   else                 ; dpf  = dpflux(lvol)
   endif

  endif ! end of if( Lplasmaregion ) ;

  dpsi(1:2) = (/  dtf,  dpf /) ! enclosed poloidal fluxes and their derivatives;
  tpsi(1:2) = (/  one, zero /) ! enclosed toroidal fluxes and their derivatives;
  ppsi(1:2) = (/ zero,  one /) ! enclosed toroidal fluxes and their derivatives;



  diotadxup(0:1,-1:2,lvol) = zero ! rotational-transform, and its derivatives with respect to lmu and dpf, or toroidal current, on the inner/outer interface;
  dItGpdxtp(0:1,-1:2,lvol) = zero ! plasma and linking currents;



  NN = NAdof(lvol) ! shorthand;

  SALLOCATE( rhs   , (1:NN,0:2 ), zero )
    SALLOCATE( matrix, (1:NN,1:NN), zero )


  solution(1:NN,-1:2) = zero ! this is a global array allocated in dforce;

  select case (Lmatsolver)
  case (1) ! direct matrix solver
    Lwork = NB*NN

    SALLOCATE( RW,    (1:Lwork ),  zero )
    SALLOCATE( RD,    (1:NN,0:2),  zero )
    SALLOCATE( LU,    (1:NN,1:NN), zero )
    SALLOCATE( ipiv,  (1:NN),         0 )
    SALLOCATE( Iwork, (1:NN),         0 )
  case (2:3) ! GMRES
    if (LILUprecond) then
      NS = NdMAS(lvol) ! shorthand
      SALLOCATE( matrixS, (1:NS), zero )

      if (Lcoordinatesingularity) then
        maxfil = Lrad(lvol) + 10
        if (NOTstellsym) maxfil = maxfil + Lrad(lvol) + 10
      else
        maxfil = 2 * Lrad(lvol) + 10
        if (NOTstellsym) maxfil = maxfil + 2 * Lrad(lvol) + 10
      end if

      Nbilut = (2*maxfil+2)*NN
      SALLOCATE( bilut, (1:Nbilut), zero)
      SALLOCATE( jbilut, (1:Nbilut), 0)
      SALLOCATE( ibilut, (1:NN+1), 0)

    endif
    nw = (NN+3)*(nrestart+2) + (nrestart+1)*nrestart
    SALLOCATE( wk, (1:nw), zero)
    SALLOCATE( jw, (1:2*NN), 0)
    SALLOCATE( iperm, (1:2*NN), 0)
  end select


  idgetrf(0:1) = 0 ! error flags;
  idgetrs(0:1) = 0 ! error flags;
  idgerfs(0:1) = 0 ! error flags;
  idgecon(0:1) = 0 ! error flags;




  do ideriv = 0, 1 ! loop over derivatives;

   if( iflag.eq.1 .and. ideriv.eq.1 ) cycle ! only need to return function; recall the derivative estimate requires function evaluation;

   if( Lcoordinatesingularity ) then

      ;matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)



    ;select case( ideriv )
    ;case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2)                       , dpsi(1:2) )
    ;case( 1 )    ! construct dMD*solution
    ; ;           ; rhs(1:NN,1) =                                                              - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )
    ; ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2)                       , ppsi(1:2) )
    ;end select

   else ! .not.Lcoordinatesingularity;

    if( Lplasmaregion ) then

       matrix(1:NN,1:NN) = dMA(1:NN,1:NN) - lmu * dMD(1:NN,1:NN)


    ;select case( ideriv )
    ;case( 0 )    ; rhs(1:NN,0) = - matmul(  dMB(1:NN,1:2)                       , dpsi(1:2) )
    ;case( 1 )    ! construct dMD*solution
    ; ;           ; rhs(1:NN,1) =                                                              - matmul( - one  * dMD(1:NN,1:NN), solution(1:NN,0) )

    ; ;           ; rhs(1:NN,2) = - matmul(  dMB(1:NN,1:2)                       , ppsi(1:2) )
    ;end select

    else ! Lvacuumregion ;

      matrix(1:NN,1:NN) = dMA(1:NN,1:NN) ! - lmu * dMD(1:NN,1:NN) ;

     select case( ideriv )
     case( 0 )    ; rhs(1:NN,0) = - dMG(1:NN) - matmul( dMB(1:NN,1:2), dpsi(1:2) ) ! perhaps there is an lmu term missing here;
     case( 1 )    ; rhs(1:NN,1) =             - matmul( dMB(1:NN,1:2), tpsi(1:2) ) ! perhaps there is an lmu term missing here;
      ;           ; rhs(1:NN,2) =             - matmul( dMB(1:NN,1:2), ppsi(1:2) ) ! perhaps there is an lmu term missing here;
     end select

    endif ! end of if( Lplasmaregion ) ;

   endif ! end of if( Lcoordinatesingularity ) ;



   select case( Lmatsolver )

   case(1) ! Using direct matrix solver (LU factorization), must not be matrix free

    select case( ideriv )

    case( 0 ) ! ideriv=0;


      

      

      MM = 1
      call DCOPY(NN*NN, matrix, 1, LU, 1) ! BLAS version
      solution(1:NN,0   ) = rhs(:,0   )
      call DGETRF(NN, NN, LU, NN, ipiv, idgetrf(ideriv) ) ! LU factorization

      anorm=maxval(sum(abs(matrix),1))
      call DGECON('I', NN, LU, NN, anorm, rcond, RW, Iwork, idgecon(ideriv)) ! estimate the condition number

      call DGETRS('N', NN, MM, LU, NN, ipiv, solution(1:NN,0   ), NN, idgetrs(ideriv) ) ! sovle linear equation
      call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1,0), NN, solution(1,0), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv)) ! refine the solution

      
    case( 1 ) ! ideriv=1;

      MM = 2
      solution(1:NN,1:MM) = rhs(:,1:MM)
      call DGETRS( 'N', NN, MM, LU, NN, ipiv, solution(1:NN,1:MM), NN, idgetrs(ideriv) )
      call DGERFS('N', NN, MM, matrix, NN, LU, NN, ipiv, rhs(1,1), NN, solution(1,1), NN, FERR, BERR, RW, Iwork, idgerfs(ideriv))

    end select ! ideriv;

    cput = GETTIME

    if(     idgetrf(ideriv) .eq. 0 .and. idgetrs(ideriv) .eq. 0 .and. idgerfs(ideriv) .eq. 0 .and. rcond .ge. machprec) then
      if( Wmp00ac ) write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "success ;         ", cput-lcpu
    elseif( idgetrf(ideriv) .lt. 0 .or. idgetrs(ideriv) .lt. 0 .or. idgerfs(ideriv) .lt. 0   ) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "input error ;     "
    elseif( idgetrf(ideriv) .gt. 0 ) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "singular ;        "
    elseif( rcond .le. machprec) then
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "ill conditioned ; "
    else
      ;             write(ounit,1010) cput-cpus, myid, lvol, ideriv, "idgetrf idgetrs idgerfs", idgetrf(ideriv), idgetrs(ideriv), idgetrf(ideriv), "invalid error ; "
    endif

  1010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; ideriv="i2" ; "a23"=",i3,' ',i3,' ',i3," ; "a34,:" time=",f10.2," ;")

   end select ! Lmatsolver

  enddo ! end of do ideriv;


  do ideriv = 0, 2

   if( iflag.eq.1 .and. ideriv.gt.0 ) cycle

   packorunpack = 'U'
   WCALL( mp00ac, packab, ( packorunpack, lvol, NN, solution(1:NN,ideriv), ideriv ) ) ! unpacking; this assigns oAt, oAz through common;


  enddo ! do ideriv = 0, 2;

   lBBintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMA(1:NN,1:NN), solution(1:NN,0) ) ) &
                     +        sum( solution(1:NN,0) * matmul( dMB(1:NN,1: 2),     dpsi(1: 2  ) ) ) !

   lABintegral(lvol) = half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) !

  DALLOCATE( matrix )
  DALLOCATE( rhs    )

  select case (Lmatsolver)
  case (1) ! LU
      DALLOCATE( RW )
      DALLOCATE( RD )
      DALLOCATE( LU )
      DALLOCATE( ipiv )
      DALLOCATE( Iwork )
  case (2:3) ! GMRES
    if (LILUprecond) then
      DALLOCATE( matrixS )
      DALLOCATE( bilut )
      DALLOCATE( jbilut )
      DALLOCATE( ibilut )
    endif
    DALLOCATE( wk )
    DALLOCATE( jw )
    DALLOCATE( iperm )
  end select


  idgetrf(0:1) = abs(idgetrf(0:1)) + abs(idgetrs(0:1)) + abs(idgerfs(0:1)) + abs(idgecon(0:1))
  if( idgetrf(0).ne.0 .or. idgetrf(1).ne.0 ) then ! failed to construct Beltrami/vacuum field and/or derivatives;

   ImagneticOK(lvol) = .false. ! set error flag;

   if( iflag.eq.1 ) Fdof(1:Ndof       ) = zero ! provide dummy intent out;
   if( iflag.eq.2 ) Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   iflag = -1 ! this value will be returned by C05PCF to ma02aa;

   goto 9999

  else

   ImagneticOK(lvol) = .true. ! set error flag; used in dforce;

  endif



  if( YESstellsym ) then ; lmns = 1 + (mns-1)           ! number of independent degrees of freedom in angle transformation;
  else                   ; lmns = 1 + (mns-1) + (mns-1) ! only required for dense, Fourier angle transformation;
  endif






  select case( Lconstraint )

  case( -1 ) ! Lconstraint=-1;

   if( Lplasmaregion ) then

    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif

    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint=-1 indicates no iterations over mu   , dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   else ! Lvacuumregion

    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif

    if( Wcurent ) then ! compute enclosed currents    only for diagnostic purposes;
     WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )
     curtor = dItGpdxtp(0,0,lvol) ! icurrent(0) ! update input variables;
     curpol = dItGpdxtp(1,0,lvol) ! gcurrent(0)
    endif

    Fdof(1:Ndof       ) = zero ! provide dummy intent out;Lconstraint=-1 indicates no iterations over dtflux, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   endif ! end of if( Lplasmaregion) ;

  case(  0 ) ! Lconstraint= 0;

   if( Lplasmaregion ) then

    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif

    Fdof(1:Ndof       ) = zero ! provide dummy intent out; Lconstraint= 0 indicates no iterations over mu, dpflux are required;
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   else ! Lvacuumregion

    WCALL( mp00ac, curent,( lvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )

    if( iflag.eq.1 ) Fdof(1:2  ) = (/ dItGpdxtp(0,0,lvol) - curtor, dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ dItGpdxtp(0,1,lvol)         , dItGpdxtp(1,1,lvol)          /)
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ dItGpdxtp(0,2,lvol)         , dItGpdxtp(1,2,lvol)          /)

   endif ! end of if( Lplasmaregion) ;

  case(  1 ) ! Lconstraint= 1;

   WCALL( mp00ac, tr00ab,( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) ) ! required for both plasma and vacuum region;

   if( Lplasmaregion ) then

    if( Lcoordinatesingularity ) then ! Ndof = 1;
     if( iflag.eq.1 ) Fdof(1    ) = diotadxup(1,  0,lvol) - iota(        lvol )
     if( iflag.eq.2 ) Ddof(1  ,1) = diotadxup(1,  1,lvol)                        ! derivative of outer rotational-transform wrt helicity multiplier
    endif

    if( Ndof.eq.2 ) then
     if( iflag.eq.1 ) Fdof(1:2  ) = diotadxup(0:1,0,lvol) - (/ oita(lvol-1), iota(lvol) /)
     if( iflag.eq.2 ) Ddof(1:2,1) = diotadxup(0:1,1,lvol)
     if( iflag.eq.2 ) Ddof(1:2,2) = diotadxup(0:1,2,lvol)
    endif

   else ! Lvacuumregion

    WCALL( mp00ac, curent, ( lvol, mn,     Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,lvol) ) )

    curtor = dItGpdxtp(0,0,lvol) ! update input variables; 08 Jun 16;

    if( iflag.eq.1 ) Fdof(1:2  ) = (/ diotadxup(0,0,lvol) - oita(lvol-1), dItGpdxtp(1,0,lvol) - curpol /)
    if( iflag.eq.2 ) Ddof(1:2,1) = (/ diotadxup(0,1,lvol)               , dItGpdxtp(1,1,lvol)          /)
    if( iflag.eq.2 ) Ddof(1:2,2) = (/ diotadxup(0,2,lvol)               , dItGpdxtp(1,2,lvol)          /)

   endif ! end of if( Lplasmaregion) ;

  case(  2 )

   if ( iflag.eq.1 ) Fdof(1     ) = lABintegral(lvol) - helicity(lvol)


    if ( iflag.eq.2 ) Ddof(1   ,1) = half * sum( solution(1:NN,1) * matmul( dMD(1:NN,1:NN), solution(1:NN,0) ) ) &
                                    + half * sum( solution(1:NN,0) * matmul( dMD(1:NN,1:NN), solution(1:NN,1) ) )

  case(  3 )

   if( Lplasmaregion ) then

    if( Wtr00ab ) then ! compute rotational transform only for diagnostic purposes;
     WCALL( mp00ac, tr00ab, ( lvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,lvol) ) )
    endif

    Fdof(1:Ndof       ) = zero ! provide dummy intent out; no iteration other mu and psip locally
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   else ! Lvacuumregion



    Fdof(1:Ndof       ) = zero ! provide dummy intent out; no iteration other mu and psip locally
    Ddof(1:Ndof,1:Ndof) = zero ! provide dummy intent out;

   endif ! end of if( Lplasmaregion) ;

  end select ! end of select case( Lconstraint ) ;



  if( Wmp00ac .or. Wma02aa ) then ! the following is screen output;

   cput = GETTIME

   if( Lplasmaregion ) then
    select case( iflag )
    case( 0 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3000) cput-cpus, myid, lvol, lmu, dpf, iflag, Fdof(1:Ndof)
    case( 2 )    ; write(ounit,3010) cput-cpus, myid, lvol, lmu, dpf, iflag, Ddof(1:Ndof,1:Ndof)
    case default ; FATAL( mp00ac, .true., illegal iflag on entry )
    end select
   else ! Lvacuumregion
    select case( iflag )
    case( 0 )    ; write(ounit,3001) cput-cpus, myid, lvol, dtf, dpf, iflag                         ! this is impossible by above logic;
    case( 1 )    ; write(ounit,3001) cput-cpus, myid, lvol, dtf, dpf, iflag, Fdof(1:Ndof)
    case( 2 )    ; write(ounit,3011) cput-cpus, myid, lvol, dtf, dpf, iflag, Ddof(1:Ndof,1:Ndof)
    case default ; FATAL( mp00ac, .true., illegal iflag on entry )
    end select
   endif

3000 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3010 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (mu,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")
3001 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" F="2es13.05" ;")
3011 format("mp00ac : ",f10.2," : myid=",i3," ; lvol=",i3," ; (dt,dp)=("es23.15" ,"es23.15" ) ; iflag="i2" ;":" D="4es13.05" ;")

  endif ! end of if( Wmp00ac .or. Wma02aa ) ;



  if( iflag.eq.1 ) then ! only in this case is Fdof defined;

   if( sum( abs( Fdof(1:Ndof) ) ) / Ndof .lt. mupftol ) then ! satisfactory;

    if ( Lplasmaregion ) then ; mu(lvol) = lmu  ;                    ; dpflux(lvol) = dpf

    else                      ; mu(lvol) = zero ; dtflux(lvol) = dtf ; dpflux(lvol) = dpf

    endif

    iflag = -2 ! return "acceptance" flag through to ma02aa via ifail; early termination;

   endif ! end of if( sum(Fdof) ) ;

  endif ! end of if( iflag.eq.1 ) ;



  RETURN(mp00ac)



end subroutine mp00ac
