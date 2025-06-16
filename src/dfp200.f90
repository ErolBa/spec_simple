









subroutine dfp200( LcomputeDerivatives, vvol)

  use constants, only : zero, half, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp200, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        epsilon, &
                        Lfindzero, &
                        Lconstraint, Lcheck, LHmatrix, &
                        Lextrap

  use cputiming, only : Tdfp200

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC,&
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                        Mvol, &
                        Iquad, &                  ! convenience; provided to ma00aa as argument to avoid allocations;
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NAdof, &
                        YESstellsym, NOTstellsym, &
                        mn, im, in, mns, Ntz, &
                        Ate, Aze, Ato, Azo, & ! only required for debugging;
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace;
                        dBdX, &
                        dtflux, dpflux, sweight, &
                        mmpp, &
                        dMA, dMB, dMD, dMG, &
                        Bemn, Bomn, Iomn, Iemn, Somn, Semn, &
                        LGdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        diotadxup, dItGpdxtp, &
                        dFFdRZ, dBBdmp, dmupfdx, denergydrr,denergydzr, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        lmns, &
                        mn, mne, &
                        dRodR, dRodZ, dZodR, dZodZ, &
                        LocalConstraint, solution, &
                        IsMyVolume, IsMyVolumeValue, Btemn, WhichCpuID

  use typedefns



  LOCALS

  LOGICAL, intent(in)  :: LComputeDerivatives ! indicates whether derivatives are to be calculated;
  LOGICAL              :: LInnerVolume

  INTEGER              :: NN, IA, ifail, if01adf, vflag, MM, idgetrf, idgetri, Lwork, lvol, pvol
  INTEGER              :: vvol, innout, ii, jj, irz, issym, iocons, idoc, idof, imn, ll
  INTEGER              :: Lcurvature, ideriv, id
  INTEGER              :: iflag, cpu_id, cpu_id1, even_or_odd, vol_parity
  INTEGER              :: stat(MPI_STATUS_SIZE), tag, tag2, req1, req2, req3, req4

  REAL                 :: lastcpu, lss, lfactor, DDl, MMl
  REAL                 :: det
  REAL   , allocatable :: XX(:), YY(:), dBB(:,:), dII(:), dLL(:), dPP(:), length(:), dRR(:,:), dZZ(:,:), constraint(:)
  REAL   , allocatable :: ddFcol1(:), ddFcol2(:), ddFcol3(:), ddFcol4(:)


  CHARACTER            :: packorunpack

  type(MatrixLU)  :: oBI(1:Mvol)

  


  SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
  SALLOCATE(  XX       , (1:Ntz     ), zero )
  SALLOCATE(  YY       , (1:Ntz     ), zero )
  SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;
  if( LocalConstraint ) then

  do vvol = 1, Mvol

    WCALL(dfp200, IsMyVolume, (vvol))

    if( IsMyVolumeValue .EQ. 0 ) then
      cycle
    else if( IsMyVolumeValue .EQ. -1) then
      FATAL(dfp200, .true., Unassociated volume)
    endif


    LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

    dBdX%vol = vvol  ! Label
    ll = Lrad(vvol)  ! Shorthand
    NN = NAdof(vvol) ! shorthand;



    vflag = 1
    WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;



    do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;

      if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
      if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;

      ideriv = 0 ; id = ideriv
      iflag = 0 ! XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
      WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dBB(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

    enddo ! end of do iocons = 0, 1;

  enddo ! matches do vvol = 1, Mvol

else ! CASE SEMI GLOBAL CONSTRAINT

     do vvol = 1, Mvol
         WCALL(dfp200, IsMyVolume, (vvol))

         if( IsMyVolumeValue .EQ. 0 ) then
             cycle
         else if( IsMyVolumeValue .EQ. -1) then
             FATAL(dfp200, .true., Unassociated volume)
         endif


         LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
         ll = Lrad(vvol)  ! Shorthand
         NN = NAdof(vvol) ! shorthand;



         vflag = 1
         WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute volume;

         do iocons = 0, 1 ! construct field magnitude on inner and outer interfaces; inside do vvol;

             if( vvol.eq.1    .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis);
             if( vvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; there are no constraints at outer boundary;

             ideriv = 0 ; id = ideriv
             iflag = 0 ! dAt, dAz, XX & YY are returned by lforce; Bemn(1:mn,vvol,iocons), Iomn(1:mn,vvol) etc. are returned through global;
             WCALL( dfp200, lforce, ( vvol, iocons, ideriv, Ntz, dBB(1:Ntz,id), XX(1:Ntz), YY(1:Ntz), length(1:Ntz), DDl, MMl, iflag ) )

         enddo ! end of do iocons = 0, 1;
     enddo

  endif ! End of if( LocalConstraint )

  DALLOCATE(dBB)
  DALLOCATE( XX) ! spectral constraints; not used;
  DALLOCATE( YY)
  DALLOCATE(length)


2000 continue



end subroutine dfp200