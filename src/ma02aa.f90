
subroutine ma02aa( lvol, NN )



  use constants, only : zero, half, one, ten

  use numerical, only : vsmall, small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wma02aa, &
                        Lconstraint, mu, helicity, &
                        mupftol, mupfits, Lrad, Lcheck

  use cputiming

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        Mvol, mn, im, in, &
                        LBlinear, LBnewton, LBsequad, &
                        dMA, dMB,      dMD,           solution, &
                        MBpsi,  Ate,                          &
                        ImagneticOK, &
                        lBBintegral, lABintegral, &
                        ivol, Nfielddof, &
                        dtflux, dpflux, &
                        xoffset, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint



  LOCALS

  INTEGER, intent(in)  :: lvol, NN


  INTEGER              :: ideriv
  REAL                 :: tol, dpsi(1:2), lastcpu
  CHARACTER            :: packorunpack

  INTEGER              :: Nxdof, Ndof, Ldfjac, iflag, maxfev, mode, LRR, nfev, njev, nprint, ihybrj
  REAL                 :: Xdof(1:2), Fdof(1:2), Ddof(1:2,1:2), oDdof(1:2,1:2)
  REAL                 :: factor, diag(1:2), RR(1:2*(2+1)/2), QTF(1:2), wk(1:2,1:4)

  INTEGER              :: irevcm

  INTEGER              :: pNN

  REAL                 :: xi(0:NN), Fxi(0:NN), xo(0:NN), Mxi(1:NN)

  external             :: mp00ac

  INTEGER              :: ihybrj1, Ldfmuaa, lengthwork
  REAL                 :: NewtonError
  REAL   , allocatable :: DFxi(:,:), work(:)
  external             :: df00ab

  INTEGER              :: NLinearConstraints, NNonLinearConstraints, LDA, LDCJ, LDR, iterations, LIWk, LRWk, ie04uff
  INTEGER, allocatable :: Istate(:), NEEDC(:), IWk(:)
  REAL                 :: objectivefunction
  REAL   , allocatable :: LinearConstraintMatrix(:,:), LowerBound(:), UpperBound(:)
  REAL   , allocatable :: constraintfunction(:), constraintgradient(:,:), multipliers(:), objectivegradient(:), RS(:,:), RWk(:)
  CHARACTER            :: optionalparameter*33

  


  ivol = lvol ! various subroutines (e.g. mp00ac, df00ab) that may be called below require volume identification, but the argument list is fixed by NAG;




  if( LBsequad ) then ! sequential quadratic programming (SQP); construct minimum energy with constrained helicity;
   lastcpu = GETTIME

   NLinearConstraints = 0 ! no linear constraints;

   NNonLinearConstraints = 1 ! single non-linear constraint = conserved helicity;

   LDA = max(1,NLinearConstraints)

   LDCJ = max(1,NNonLinearConstraints)

   LDR = NN

   SALLOCATE( LinearConstraintMatrix, (1:LDA,1:1), zero ) ! linear constraint matrix;

   SALLOCATE( LowerBound, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! lower bounds on variables, linear constraints and non-linear constraints;
   SALLOCATE( UpperBound, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! upper bounds on variables, linear constraints and non-linear constraints;

   LowerBound(                       1 : NN                                          ) = -1.0E+21       !   variable constraints; no constraint;
   UpperBound(                       1 : NN                                          ) = +1.0E+21       !
   LowerBound( NN+                   1 : NN+NLinearConstraints                       ) = -1.0E+21       !     linear constraints; no constraint;
   UpperBound( NN+                   1 : NN+NLinearConstraints                       ) = +1.0E+21       !
   LowerBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) ! non-linear constraints; enforce helicity constraint;
   UpperBound( NN+NLinearConstraints+1 : NN+NLinearConstraints+NNonLinearConstraints ) = helicity(lvol) !

   iterations = 0 ! iteration counter;

   SALLOCATE( Istate, (1:NN+NLinearConstraints+NNonLinearConstraints), 0 )

   SALLOCATE( constraintfunction, (1:NNonLinearConstraints), zero ) ! constraint functions;

   SALLOCATE( constraintgradient, (1:LDCJ,1:NN), zero ) ! derivatives of constraint functions;

   SALLOCATE( multipliers, (1:NN+NLinearConstraints+NNonLinearConstraints), zero ) ! Lagrange multipliers ?;

   objectivefunction = zero ! objective function;

   SALLOCATE( objectivegradient, (1:NN), zero ) ! derivatives of objective function;

   SALLOCATE( RS, (1:LDR,1:NN), zero )
   ideriv = 0 ; dpsi(1:2) = (/ dtflux(lvol), dpflux(lvol) /) ! these are also used below;

   packorunpack = 'P'

   CALL( ma02aa, packab, ( packorunpack, lvol, NN, xi(1:NN), ideriv ) )

   SALLOCATE( NEEDC, (1:NNonLinearConstraints), 0 )

   LIWk = 3*NN + NLinearConstraints + 2*NNonLinearConstraints ! workspace;
   SALLOCATE( IWk, (1:LIWk), 0 )       ! workspace;

   LRWk = 2*NN**2 + NN * NLinearConstraints + 2 * NN * NNonLinearConstraints + 21 * NN + 11 * NLinearConstraints + 22 * NNonLinearConstraints + 1 ! workspace;
   SALLOCATE( RWk, (1:LRWk), zero )                                                              ! workspace;

   irevcm = 0 ; ie04uff = 1 ! reverse communication loop control; ifail error flag;




   MBpsi(1:NN) =                         matmul( dMB(1:NN,1: 2), dpsi(1:2) )






   DALLOCATE(RWk)
   DALLOCATE(IWk)
   DALLOCATE(NEEDC)
   DALLOCATE(RS)
   DALLOCATE(objectivegradient)
   DALLOCATE(multipliers)
   DALLOCATE(constraintgradient)
   DALLOCATE(constraintfunction)
   DALLOCATE(Istate)
   DALLOCATE(LowerBound)
   DALLOCATE(UpperBound)
   DALLOCATE(LinearConstraintMatrix)


   packorunpack = 'U'
   CALL( ma02aa, packab ( packorunpack, lvol, NN, xi(1:NN), ideriv ) )

   lBBintegral(lvol) = half * sum( xi(1:NN) * matmul( dMA(1:NN,1:NN), xi(1:NN) ) ) + sum( xi(1:NN) * MBpsi(1:NN) ) ! + psiMCpsi
   lABintegral(lvol) = half * sum( xi(1:NN) * matmul( dMD(1:NN,1:NN), xi(1:NN) ) ) ! + sum( xi(1:NN) * MEpsi(1:NN) ) ! + psiMFpsi

   solution(1:NN,0) = xi(1:NN)


  endif ! end of if( LBsequad ) then;

  if( LBlinear ) then ! assume Beltrami field is parameterized by helicity multiplier (and poloidal flux);

   lastcpu = GETTIME

   if( Lplasmaregion ) then

    Xdof(1:2) = xoffset + (/     mu(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

    select case( Lconstraint )
    case( -1 )    ;                                   ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
    ;             ; iflag = 1 !(we don't need derivatives)
    case(  0 )    ;                                   ; Nxdof = 0 ! multiplier & poloidal flux NOT varied
    ;             ; iflag = 1 !(we don't need derivatives)                            ;
    case(  1 )    ; if( Lcoordinatesingularity ) then ; Nxdof = 1 ! multiplier                 IS  varied to match       outer transform;
     ;              else                              ; Nxdof = 2 ! multiplier & poloidal flux ARE varied to match inner/outer transform;
     ;              endif
    case(  2 )    ;                                     Nxdof = 1 ! multiplier                 IS  varied to match             helicity ;
    case(  3 )    ; if( Lcoordinatesingularity ) then ; Nxdof = 0 ! multiplier & poloidal flux NOT varied                               ;
     ;              else                              ; Nxdof = 0 ! Global constraint, no dof locally
     ;              endif
     ;            ; iflag = 2 !(we still need derivatives)
    end select

   else ! Lvacuumregion ;

    Xdof(1:2) = xoffset + (/ dtflux(lvol), dpflux(lvol) /) ! initial guess for degrees of freedom; offset from zero so that relative error is small;

    select case( Lconstraint )
    case( -1 )    ;                                   ; Nxdof = 0 ! poloidal   & toroidal flux NOT varied                                                  ;
    case(  0 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current      ;
    case(  1 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and transform-constraint;
    case(  2 )    ;                                   ; Nxdof = 2 ! poloidal   & toroidal flux ARE varied to match linking current and plasma current      ;
    case(  3 )    ;                                   ; Nxdof = 0 ! Fluxes are determined in dforce via a linear system
                                                      ; iflag = 2
    end select

   endif ! end of if( Lplasmaregion) ;

   select case( Nxdof )

   case( 0   ) ! need only call mp00ac once, to calculate Beltrami field for given helicity multiplier and enclosed fluxes;

    ;         ; Ndof = 1     ; Ldfjac = Ndof ; nfev = 1 ; njev = 0 ; ihybrj = 1;  ! provide dummy values for consistency;

    WCALL( ma02aa, mp00ac, ( Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) )

    helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;


   case( 1:2 ) ! will iteratively call mp00ac, to calculate Beltrami field that satisfies constraints;

    ;         ; Ndof = Nxdof ; Ldfjac = Ndof ; nfev = 0 ; njev = 0 ; ihybrj = 0;

    tol = mupftol ; LRR = Ndof * ( Ndof+1 ) / 2 ; mode = 0 ; diag(1:2) = zero ; factor = one ; maxfev = mupfits ; nprint = 0

    FATAL( ma02aa, Ndof.gt.2, illegal )

    WCALL( ma02aa, hybrj2, ( mp00ac, Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, tol, &
                             maxfev, diag(1:Ndof), mode, factor, nprint, ihybrj, nfev, njev, RR(1:LRR), LRR, QTF(1:Ndof), &
                 WK(1:Ndof,1), WK(1:Ndof,2), WK(1:Ndof,3), WK(1:Ndof,4) ) )

    if( Lplasmaregion ) then

     select case( ihybrj )
     case( 0: ) ;     mu(lvol) = Xdof(1)      - xoffset
      ;         ; dpflux(lvol) = Xdof(2)      - xoffset
     case( :-1) ;      Xdof(1) = mu(lvol)     + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
      ;         ;      Xdof(2) = dpflux(lvol) + xoffset ! mu    and dpflux have been updated in mp00ac; early termination;
     end select

    else ! Lvacuumregion;

     select case( ihybrj )
     case( 0: ) ; dtflux(lvol) = Xdof(1)      - xoffset
      ;         ; dpflux(lvol) = Xdof(2)      - xoffset
     case( :-1) ; Xdof(1)      = dtflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
      ;         ; Xdof(2)      = dpflux(lvol) + xoffset ! dtflux and dpflux have been updated in mp00ac; early termination;
     end select

    endif ! end of if( Lplasmaregion ) ;

    if (Lconstraint .ne. 2) helicity(lvol) = lABintegral(lvol) ! this was computed in mp00ac;



    if( Lconstraint.eq.1 .or. Lconstraint.eq.3 .or. ( Lvacuumregion .and. Lconstraint.eq.0 ) ) then

     iflag = 2 ; Ldfjac = Ndof ! call mp00ac: tr00ab/curent to ensure the derivatives of B, transform, currents, wrt mu/dtflux & dpflux are calculated;

     WCALL( ma02aa, mp00ac, ( Ndof, Xdof(1:Ndof), Fdof(1:Ndof), Ddof(1:Ldfjac,1:Ndof), Ldfjac, iflag ) )

    endif ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.eq.0 ) ) ;



   end select ! end of select case( Nxdof ) ;


   cput = GETTIME

   select case(ihybrj) ! this screen output may not be correct for Lvacuumregion;
   case(    1   )
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "success         ", Fdof(1:Ndof)
   case(   -2   )
    if( Wma02aa ) write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "|F| < mupftol   ", Fdof(1:Ndof)
   case(   -1   )
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "Beltrami fail   ", Fdof(1:Ndof)
   case(    0   )
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "input error     ", Fdof(1:Ndof)
   case(    2   )
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "consider restart", Fdof(1:Ndof)
   case(    3   )
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "xtol too small  ", Fdof(1:Ndof)
   case(    4:5 )
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "bad progress    ", Fdof(1:Ndof)
   case default
    ;             write(ounit,1040) cput-cpus, myid, lvol, ihybrj, helicity(lvol), mu(lvol), dpflux(lvol), cput-lastcpu, "illegal ifail   ", Fdof(1:Ndof)
    FATAL( ma02aa, .true., illegal ifail returned by hybrj )
   end select

  endif ! end of if( LBlinear ) then;





1010 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; SQP    : ie04uff=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ":,a36)
1020 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Newton : ihybrj1=",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  "error="es7.0" ; ":,a18)
1040 format("ma02aa : ",f10.2," : myid=",i3," ; lvol=",i3," ; Linear : ihybrj =",i3," hel="es12.4" mu="es12.4" dpflux="es12.4" time="f9.1" ; "&
  :,a16" ; F="2es08.0)



end subroutine ma02aa


