









subroutine dfp100(Ndofgl, x, Fvec, LComputeDerivatives)

  use constants, only : zero, half, one, two, pi2, pi, mu0

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp100, Igeometry, Nvol, Lrad, Isurf, &
                        Lconstraint, curpol

  use cputiming, only : Tdfp100

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        ImagneticOK, NAdof, mn, &
                        Mvol, Iquad, &
                        dBdX, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, Localconstraint, &
                        IPDt, IPDtdPf, xoffset, dpflux, &
                        IsMyVolume, IsMyVolumeValue, WhichCpuID, &
                        IconstraintOK, &
                        DToocc, DToocs, DToosc, DTooss, &
                        TTsscc, TTsscs, TTsssc, TTssss, &
                        TDstcc, TDstcs, TDstsc, TDstss, &
                        TDszcc, TDszcs, TDszsc, TDszss, &
                        DDttcc, DDttcs, DDttsc, DDttss, &
                        DDtzcc, DDtzcs, DDtzsc, DDtzss, &
                        DDzzcc, DDzzcs, DDzzsc, DDzzss, &
                        dMA, dMB, dMD, dMG, MBpsi, solution, &
                        Nt, Nz, LILUprecond, Lsavedguvij, guvijsave, izbs

  LOCALS

  INTEGER              :: vvol, Ndofgl, iflag, cpu_send_one, cpu_send_two, ll, NN, ideriv, iocons
  INTEGER              :: status(MPI_STATUS_SIZE), request1, request2
  REAL                 :: Fvec(1:Ndofgl), x(1:Mvol-1), Bt00(1:Mvol, 0:1, -1:2), ldItGp(0:1, -1:2)
  LOGICAL              :: LComputeDerivatives
  INTEGER              :: deriv, Lcurvature



  BEGIN(dfp100)

  dpflux(2:Mvol) = x - xoffset

  do vvol = 1, Mvol

    LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
    ImagneticOK(vvol) = .false.
    

    call IsMyVolume(vvol)

    if( IsMyVolumeValue .EQ. 0 ) then
      cycle
    else if( IsMyVolumeValue .EQ. -1) then
      FATAL(dfp100, .true., Unassociated volume)
    endif

    NN = NAdof(vvol)
    ll = Lrad(vvol)

    call allocate_geometry_matrices(vvol, LcomputeDerivatives)
    call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
    call intghs_workspace_init(vvol)

    dBdX%L = .false. ! No need to take derivatives of matrices w.r.t geometry.

    ideriv = 0 ; Lcurvature = 1
    WCALL( dfp100, compute_guvijsave, (Iquad(vvol), vvol, ideriv, Lcurvature) )
    Lsavedguvij = .true.

    if (LILUprecond) then
      WCALL( dfp100, spsint, ( Iquad(vvol), mn, vvol, ll ) )
      WCALL( dfp100, spsmat, ( vvol, mn, ll) )
    endif

    WCALL( dfp100, ma00aa, ( Iquad(vvol), mn, vvol, ll ) ) ! compute volume integrals of metric elements;
    WCALL( dfp100, matrix, ( vvol, mn, ll ) )

    WCALL( dfp100, ma02aa, ( vvol, NN ) )

    Lsavedguvij = .false.



    call intghs_workspace_destroy()
    call deallocate_Beltrami_matrices(LcomputeDerivatives)
    call deallocate_geometry_matrices(LcomputeDerivatives)

    if( Lconstraint.EQ.3 ) then
      do iocons=0,1
        if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle

        ideriv = 0
        WCALL( dfp100, lbpol, (vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons) ) !Compute field at interface for global constraint

        ideriv = 2
        WCALL( dfp100, lbpol, (vvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons) ) !Compute field at interface for global constraint, d w.r.t. pflux
      enddo

      if( Lvacuumregion ) then



        ideriv=1 ! derivatives of Btheta w.r.t tflux
        iocons=0 ! Only need inner side of volume derivatives
        WCALL( dfp100, lbpol, (Mvol, Bt00(1:Mvol, 0:1, -1:2), ideriv, iocons) )

        iflag=2 ! derivatives of poloidal linking current w.r.t geometry not required
        WCALL( dfp100, curent, (Mvol, mn, Nt, Nz, iflag, ldItGp(0:1,-1:2) ) )
      endif
    endif
  enddo


  if( .not.LocalConstraint ) then

    select case (Lconstraint)

      case( 3 )

        do vvol = 1, Mvol-1

          call WhichCpuID(vvol  , cpu_send_one)
          call WhichCpuID(vvol+1, cpu_send_two)

          RlBCAST(Bt00(vvol  , 1, 0), 1, cpu_send_one)
          RlBCAST(Bt00(vvol+1, 0, 0), 1, cpu_send_two)
          RlBCAST(Bt00(vvol  , 1, 2), 1, cpu_send_one)
          RlBCAST(Bt00(vvol+1, 0, 2), 1, cpu_send_two)

          IPDt(vvol) = pi2 * (Bt00(vvol+1, 0, 0) - Bt00(vvol, 1, 0))

          IPDtdPf(vvol,vvol) = pi2 * Bt00(vvol+1, 0, 2)
          if (vvol .ne. 1) IPDtdPf(vvol,vvol-1) = -pi2 * Bt00(vvol, 1, 2)
        enddo

        if( myid.EQ.0 ) then
            Fvec(1:Mvol-1) = IPDt(1:Mvol-1) - Isurf(1:Mvol-1)
        endif

      case default
        FATAL(dfp100, .true., Unaccepted value for Lconstraint)
    end select
  endif

  RETURN(dfp100)

end subroutine dfp100
