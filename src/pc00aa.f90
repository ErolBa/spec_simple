
subroutine pc00aa( NGdof, position, Nvol, mn, ie04dgf ) ! argument list is optional;



  use constants, only : zero, ten

  use numerical, only :

  use fileunits, only : ounit

  use inputlist, only : Wpc00aa, verify, maxstep, maxiter, forcetol

  use cputiming, only : Tpc00aa

  use allglobal, only : ncpu, myid, cpus, Energy, ForceErr



  LOCALS

  INTEGER, intent(in)    :: Nvol, mn, NGdof
  REAL   , intent(inout) :: position(0:NGdof)
  INTEGER                :: ie04dgf

  LOGICAL                :: LComputeDerivatives!, Lexit = .true.
  INTEGER                :: niterations, Iwork(1:NGdof+1), iuser(1:2)
  REAL                   :: lEnergy, Gradient(0:NGdof), work(1:13*NGdof), ruser(1:1)
  CHARACTER              :: smaxstep*34

  external               :: pc00ab

  BEGIN(pc00aa)



  if( myid.eq.0 ) then
   cput = GETTIME
   write(ounit,'("pc00aa : ", 10x ," : ")')
   write(ounit,1000) cput-cpus, myid, NGdof, maxstep, maxiter, verify
  endif

1000 format("pc00aa : ",f10.2," : myid=",i3," ; calling E04DGF : NGdof="i6" ; maxstep="es10.2" ; maxiter="i9" ; verify=",i3," ;")



  iuser(1:2) = (/ 0, 0 /) ! this is only used to pass information through to pc00ab; some resolution settings & iteration counter;
  ruser(1:1) = zero       ! this will be assigned in first call to pc00ab;






  select case( verify ) ! supply optional parameters to E04DGF;

  case( -1 ) ! no checks; no screen output;
   call E04DKF('Nolist')
   call E04DKF('Print Level = 0')
   call E04DKF('Verify = -1')
  case(  0 ) ! simple check;
   call E04DKF('Verify =  0') ! simple check
  case(  1 ) ! extensive test;
   call E04DKF('Verify =  1') ! extensive test;
  case default
   FATAL(pc00aa, .true., invalid verify supplied on input)
  end select

  call E04DKF('Iteration Limit = 99999999')
  write(smaxstep,'("Maximum Step Length ="es13.5)')maxstep
  call E04DKF(smaxstep)



  ie04dgf = 1

  call E04DGF( NGdof, pc00ab, niterations, lEnergy, Gradient(1:NGdof), position(1:NGdof), &
               Iwork(1:NGdof+1), work(1:13*NGdof), iuser(1:2), ruser(1:1), ie04dgf )

  cput = GETTIME

  select case( ie04dgf )
  case(:-1)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : user requested termination      ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  0)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : success                         ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  3)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : iteration limit reached         ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  4)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : step length too small           ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  6)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : not all minimum conditions met  ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  7)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : user supplied derivatives error ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  8)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : initial gradient too small      ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case(  9)    ; if( myid.eq.0 ) write(ounit,'("pc00aa : ",f10.2," : input error                     ; ie04dgf=",i3," ;")')cput-cpus,ie04dgf
  case default
   FATAL(pc00aa, .true., E04DGF ifail error)
  end select



  RETURN(pc00aa)



end subroutine pc00aa


