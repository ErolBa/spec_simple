









subroutine dfp200( LcomputeDerivatives, vvol)

  use constants, only : zero, half, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp200, Nvol, Mpol, Ntor, Lrad, tflux, Igeometry, &
                        gamma, adiabatic, pscale, mu, &
                        epsilon, &
                        Lfindzero, &
                        Lconstraint, Lcheck, LHmatrix, &
                        Lextrap, &
                        Lfreebound

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
                        dFFdRZ, dBBdmp, dmupfdx, hessian, dessian, Lhessianallocated,Lhessian2Dallocated,Lhessian3Dallocated,denergydrr,denergydzr, &
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

  BEGIN(dfp200)


  SALLOCATE( dBB       , (1:Ntz,-1:2), zero ) ! magnetic field strength (on interfaces) in real space and derivatives;
  SALLOCATE(  XX       , (1:Ntz     ), zero )
  SALLOCATE(  YY       , (1:Ntz     ), zero )
  SALLOCATE( length    , (1:Ntz     ), zero ) ! this is calculated in lforce;

  if( LComputeDerivatives ) then
    SALLOCATE( dRR       , (1:Ntz,-1:1), zero )
    SALLOCATE( dZZ       , (1:Ntz,-1:1), zero )
    SALLOCATE( dII       , (1:Ntz     ), zero ) ! spectral constraint;
    SALLOCATE( dLL       , (1:Ntz     ), zero ) ! length   constraint;
    SALLOCATE( dPP       , (1:Ntz     ), zero ) ! poloidal constraint;
    SALLOCATE( constraint, (1:Ntz     ), zero )
    SALLOCATE( ddFcol1   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol2   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol3   , (1:Ntz     ), zero )
    SALLOCATE( ddFcol4   , (1:Ntz     ), zero )
  endif

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

    if( LcomputeDerivatives ) then ! compute inverse of Beltrami matrices;

      SALLOCATE( oBI(vvol)%mat, (0:NN,0:NN ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
      SALLOCATE( oBI(vvol)%ipivot, (0:NN), zero)

      call allocate_geometry_matrices(vvol, LcomputeDerivatives)
      call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)
      call intghs_workspace_init(vvol)

      packorunpack = 'P'
      WCALL( dfp200, packab, ( packorunpack, vvol, NN, solution(1:NN,0), 0 ) ) ! packing, put the solutions back to the solution matrix;



      call get_LU_Beltrami_matrices(vvol, oBI(vvol), NN)

      dBdX%L = .true. ! will need derivatives;
      idof = 0 ! labels degree of freedom; local to interface;

      do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;

        dBdX%ii = ii ! controls construction of derivatives in subroutines called below;
        do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;

          if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
          dBdX%irz = irz ! controls construction of derivatives;

          do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms;

            if( issym.eq.1 .and. YESstellsym               ) cycle ! no dependence on non-stellarator symmetric harmonics;
            if( ii.eq.1    .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
            if( ii.eq.1    .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};

            dBdX%issym = issym ! controls construction of derivatives;
            idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; this is outside do innout;


            do innout = 0, 1 ! loop over deformations to inner and outer interface; inside do vvol; inside do ii; inside do irz;

              if( vvol.eq.1    .and. innout.eq.0 ) cycle ! no degrees of freedom at coordinate axis / fixed inner boundary;
              if( vvol.eq.Mvol .and. innout.eq.1 ) cycle ! no degress of freedom                      fixed outer boundary; for linearized displacement;

              dBdX%innout = innout

              call get_perturbed_solution(vvol, oBI(vvol), NN)


              call evaluate_dmupfdx(innout, idof, ii, issym, irz)


              if(Lhessianallocated) then
                call evaluate_dBB(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz, &
                                  LComputeDerivatives)
              endif

                if (Lhessian3Dallocated .and. Igeometry.ge.3) then
                    call hessian3D_dFFdRZ(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
                endif


            enddo ! matches do innout;
          enddo ! matches do issym;
        enddo ! matches do irz;
      enddo ! matches do ii;

      dBdX%L = .false. ! probably not needed, but included anyway;

      call intghs_workspace_destroy()
      call deallocate_Beltrami_matrices(LcomputeDerivatives)
      call deallocate_geometry_matrices(LcomputeDerivatives)

      DALLOCATE(oBI(vvol)%mat)
      DALLOCATE(oBI(vvol)%ipivot)

    endif ! end of if( LComputeDerivatives ) ;

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



    if( LComputeDerivatives ) then
        dBdX%L = .true. ! will need derivatives;

        if( Lfreebound.eq.1 ) then
            do ideriv = 0, 2
                do ii = 1, mn

                    call WhichCpuID(Mvol, cpu_id)

                    RlBCAST( Ate(Mvol,ideriv,ii)%s(0:Lrad(Mvol)), Lrad(Mvol)+1, cpu_id )
                    RlBCAST( Aze(Mvol,ideriv,ii)%s(0:Lrad(Mvol)), Lrad(Mvol)+1, cpu_id )
                enddo
             enddo
        endif






        do vvol = 1, Mvol
            LREGION(vvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ; TODO: maybe not necessary, remove
            NN = NAdof(vvol) ! shorthand;
            ll = Lrad(vvol)

            SALLOCATE( oBI(vvol)%mat, (0:NN,0:NN ), zero ) ! inverse of ``original'', i.e. unperturbed, Beltrami matrix;
            SALLOCATE( oBI(vvol)%ipivot, (0:NN), zero)

            WCALL(dfp200, IsMyVolume, (vvol))
            if( IsMyVolumeValue .EQ. 0 ) then
                cycle
            else if( IsMyVolumeValue .EQ. -1) then
                FATAL(dfp200, .true., Unassociated volume)
            endif



            call allocate_geometry_matrices(vvol, LcomputeDerivatives)
            call allocate_Beltrami_matrices(vvol, LcomputeDerivatives)

            call get_LU_beltrami_matrices(vvol, oBI(vvol), NN)

            call deallocate_Beltrami_matrices(LcomputeDerivatives)
            call deallocate_geometry_matrices(LcomputeDerivatives)
        enddo ! end of vvol = 1, Mvol

        if( ncpu .gt. 1 ) then
            do vvol = 1, Mvol
                NN = NAdof(vvol)
                call WhichCpuID(vvol, cpu_id)
                call MPI_BCAST( oBI(vvol)%mat(0:NN,0:NN), (NN+1)**2, MPI_DOUBLE_PRECISION, cpu_id , MPI_COMM_SPEC , ierr)
                call MPI_BCAST( oBI(vvol)%ipivot(0:NN)  ,  NN+1    , MPI_INTEGER         , cpu_id , MPI_COMM_SPEC , ierr)
            enddo
        endif

        do even_or_odd = 0, 1 ! First loop on even interfaces perturbation, then on odd interfaces. This allow efficient parallelization

            do vvol = 1, Mvol-1 !labels which interface is perturbed


            vol_parity = MODULO(vvol,2)
            if( (vol_parity.eq.0 ) .and. (even_or_odd.eq.1) ) cycle ! even_or_odd=1 thus perturb only odd  interfaces
            if( (vol_parity.eq.1 ) .and. (even_or_odd.eq.0) ) cycle ! even_or_odd=0 thus perturb only even interfaces

            WCALL(dfp200, IsMyVolume, (vvol))

            if( IsMyVolumeValue.EQ.0 ) then ! This CPU does not deal with interface's inner volume

                WCALL(dfp200, IsMyVolume, (vvol+1))

                if( IsMyVolumeValue.eq.0 ) then ! This CPU does not deal with interface's outer volume either - cycle
                    cycle
                else if( IsMyVolumeValue.eq.-1 ) then
                    FATAL(dfp200, .true., Unassociated volume)
                else
                    LInnerVolume = .false.
                endif

            else if( IsMyVolumeValue.EQ.-1 ) then
                FATAL(dfp200, .true., Unassociated volume)
            else
                LinnerVolume = .true.
            endif


            dBdX%vol = vvol     ! Perturbed interface
            idof = 0 ! labels degree of freedom of interface vvol

            do ii = 1, mn ! loop over deformations in Fourier harmonics; inside do vvol;
            dBdX%ii = ii ! controls construction of derivatives in subroutines called below;

            do irz = 0, 1 ! loop over deformations in R and Z; inside do vvol; inside do ii;
                if( irz.eq.1 .and. Igeometry.lt.3 ) cycle ! no dependence on Z;
                dBdX%irz = irz ! controls construction of derivatives;

                do issym = 0, 1 ! loop over stellarator and non-stellarator symmetric terms;
                    if( issym.eq.1 .and. YESstellsym               ) cycle ! no dependence on non-stellarator symmetric harmonics;
                    if( ii.eq.1    .and. irz.eq.1 .and. issym.eq.0 ) cycle ! no dependence on Zbs_{m=0,n=0};
                    if( ii.eq.1    .and. irz.eq.0 .and. issym.eq.1 ) cycle ! no dependence on Rbs_{m=0,n=0};
                    dBdX%issym = issym ! controls construction of derivatives;

                    idof = idof + 1 ! this labels the degree-of-freedom that the derivative is taken with respect to; this is outside do innout;



                    do lvol = vvol, vvol+1

                        if( ncpu.gt.1 ) then
                            if(      LinnerVolume .and. (lvol.eq.vvol+1) ) cycle
                            if( .not.LinnerVolume .and. (lvol.eq.vvol  ) ) cycle
                        endif

                        if( lvol.eq.vvol   ) innout=1 ! Perturb w.r.t outer interface
                        if( lvol.eq.vvol+1 ) innout=0 ! Perturb w.r.t inner interface

                        dBdX%innout = innout
                        dBdX%L      = .true.

                        LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;
                        NN = NAdof(lvol)

                        call allocate_geometry_matrices(lvol, LcomputeDerivatives)
                        call allocate_Beltrami_matrices(lvol, LcomputeDerivatives)

                        call intghs_workspace_init(lvol)
                        call get_perturbed_solution(lvol, oBI(lvol), NN)
                        call intghs_workspace_destroy()

                        call deallocate_Beltrami_matrices(LcomputeDerivatives)
                        call deallocate_geometry_matrices(LcomputeDerivatives)

                    enddo ! end of do lvol = vvol, vvol+1


                    call WhichCpuID(vvol  , cpu_id )
                    call WhichCpuID(vvol+1, cpu_id1)


                    if( ncpu.gt. 1) then
                        if( LinnerVolume ) then
                            do jj = 1, mn
                                tag  = 1 ! Tags for MPI communications
                                tag2 = 2

                                call MPI_RECV(Ate(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag , MPI_COMM_SPEC, stat, ierr)
                                call MPI_RECV(Aze(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag2, MPI_COMM_SPEC, stat, ierr)
                            enddo

                            if( NOTstellsym ) then
                                do jj = 1, mn
                                    tag  = 3
                                    tag2 = 4

                                    call MPI_RECV(Ato(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag , MPI_COMM_SPEC, stat, ierr)
                                    call MPI_RECV(Azo(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id1, tag2, MPI_COMM_SPEC, stat, ierr)
                                enddo
                            endif

                        else
                            do jj = 1, mn
                                tag  = 1
                                tag2 = 2

                                call MPI_SEND(Ate(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag , MPI_COMM_SPEC, ierr)
                                call MPI_SEND(Aze(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag2, MPI_COMM_SPEC, ierr)
                            enddo

                            if( NOTstellsym ) then
                                do jj = 1, mn
                                    tag  = 3
                                    tag2 = 4

                                    call MPI_SEND(Ato(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag , MPI_COMM_SPEC, ierr)
                                    call MPI_SEND(Azo(vvol+1,-1,jj)%s(0:Lrad(vvol+1)), Lrad(vvol+1)+1, MPI_DOUBLE_PRECISION, cpu_id , tag2, MPI_COMM_SPEC, ierr)
                                enddo
                            endif
                        endif
                    endif

                    if( LinnerVolume) then
                        call evaluate_dmupfdx(1, idof, ii, issym, irz)
                    else
                        dmupfdx(1:Mvol, vvol, 1:2, idof, 1) = zero
                    endif

                    do lvol = vvol, vvol+1
                        WCALL(dfp200, IsMyVolume, (lvol))

                        if( IsMyVolumeValue .EQ. 0 ) then
                            cycle
                        else if( IsMyVolumeValue .EQ. -1) then
                            FATAL(dfp200, .true., Unassociated volume)
                        endif

                        if( lvol.eq.vvol ) then
                            innout      = 1
                        else
                            innout      = 0
                        endif
                        dBdX%innout = innout

                        LREGION(lvol) ! assigns Lcoordinatesingularity, Lplasmaregion, etc. ;

                        call evaluate_dBB(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz, &
                                          LComputeDerivatives)

                        if (Lhessian3Dallocated .and. Igeometry.ge.3) then
                            call hessian3D_dFFdRZ(vvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)
                        endif

                    enddo     ! matches do lvol = vvol, vvol+1

                enddo ! matches do issym;
            enddo ! matches do irz;
        enddo ! matches do ii;
      enddo ! matches do vvol;
    enddo ! matches do even_or_odd;

    do vvol = 1, Mvol
        DALLOCATE(oBi(vvol)%mat)
        DALLOCATE(oBI(vvol)%ipivot)
    enddo

    dBdX%L = .false. ! probably not needed, but included anyway;

    do vvol = 1, Mvol
        call WhichCpuID(vvol, cpu_id)

        if( vvol.ne.Mvol ) then
            call MPI_BCAST( dmupfdx(1:Mvol, vvol ,1:2, 1:LGdof,   1), Mvol*2*LGdof  , MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr )
        endif

        call MPI_BCAST( dFFdRZ(1:LGdof, 0:1, 1:LGdof, 0:1, vvol), 2*2*(LGdof**2), MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr )
        call MPI_BCAST( dBBdmp(1:LGdof, vvol ,0:1, 1:2         ), 2*2*LGdof, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_SPEC, ierr )
    enddo


    if (Lhessian3Dallocated) then
        do vvol = 1, Mvol
            call WhichCpuID(vvol, cpu_id)

            call MPI_BCAST( denergydrr(1:LGdof, vvol, 0:1, 1:LGdof, 0:1), 2*2*(LGdof**2), MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
            call MPI_BCAST( denergydzr(1:LGdof, vvol, 0:1, 1:LGdof, 0:1), 2*2*(LGdof**2), MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
            call MPI_BCAST( dBBdmp(1:LGdof, vvol ,0:1, 1:2         ), 2*2*LGdof, MPI_DOUBLE_PRECISION, cpu_id, MPI_COMM_WORLD, ierr )
        enddo
    endif

    endif ! End of if( LComputeDerivatives )



  endif ! End of if( LocalConstraint )

  if( LcomputeDerivatives ) then
      DALLOCATE(constraint)
      DALLOCATE(dPP)
      DALLOCATE(dLL)
      DALLOCATE(dII)
      DALLOCATE(dZZ)
      DALLOCATE(dRR)
      DALLOCATE(ddFcol1)
  DALLOCATE(ddFcol2)
  DALLOCATE(ddFcol3)
  DALLOCATE(ddFcol4)
  endif

  DALLOCATE(dBB)
  DALLOCATE( XX) ! spectral constraints; not used;
  DALLOCATE( YY)
  DALLOCATE(length)


2000 continue

  RETURN(dfp200)

end subroutine dfp200




subroutine get_LU_beltrami_matrices(vvol, oBI, NN)





  use constants, only :   zero, half, one, two

  use fileunits, only :   ounit

  use cputiming, only :   Tdfp200

  use inputlist, only :   Wmacros, Wdfp200, Wdforce, Lrad, mu, Lconstraint

  use allglobal, only :   ncpu, myid, cpus, MPI_COMM_SPEC, &
                          Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                          Nt, Nz, &
                          dBdX, &
                          dMA, dMB, dMD, dMG, &
                          mn, mne, Iquad, solution, Lsavedguvij

  use typedefns

  LOCALS
  TYPE(MatrixLU), intent(inout) :: oBI
  INTEGER, intent(in)    :: vvol, NN

  INTEGER             :: IA, MM, LDA, Lwork, ll
  INTEGER             :: idgetrf, idgetri
  REAL                :: lastcpu
  REAL, allocatable   :: work(:)


  lastcpu = GETTIME

  ll = Lrad(vvol)

  dBdX%L = .false.

  Lsavedguvij = .false.

  WCALL( dfp200, ma00aa, (Iquad(vvol), mn, vvol, ll) )
  WCALL( dfp200, matrix, (vvol, mn, ll) )

  lastcpu = GETTIME

  if(Lconstraint .eq. 2) then  ! for helicity constraint

    call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019
    dMA(0,0)       = zero

    dMA(1:NN,0)    = -matmul(dMD(1:NN,1:NN),solution(1:NN,0))
    dMA(0,1:NN)    = dMA(1:NN,0) ! This is the transpose of the above

    IA = NN + 1 + 1
    MM = NN ; LDA = NN+1 ; Lwork = NN+1
    idgetrf = 1 ; call DGETRF( MM+1, NN+1, dMA(0:LDA-1,0:NN), LDA, oBI%ipivot, idgetrf )

    cput = GETTIME
    select case( idgetrf ) !                                                                     0123456789012345678
    case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
    case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
    case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
    case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
    end select

    call DCOPY((1+NN)*(1+NN), dMA, 1, oBI%mat, 1)            ! BLAS version; 24 Jul 2019

  else ! for other constraints

    call DAXPY((NN+1)*(NN+1), -mu(vvol), dMD, 1, dMA, 1) ! BLAS version; 24 Jul 2019

    call DCOPY((NN+1)*(NN+1), dMA, 1, dMD, 1)            ! BLAS version; 24 Jul 2019

    IA = NN + 1

    MM = NN ; LDA = NN + 1 ; Lwork = NN

    idgetrf = 1 ; call DGETRF( MM, NN, dMA(1,1), LDA, oBI%ipivot(1:NN), idgetrf )

    cput = GETTIME
    select case( idgetrf ) !                                                                     0123456789012345678
    case(  :-1 ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "input error;      "
    case(  0   ) ; if( Wdforce ) write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "success;          "
    case( 1:   ) ;               write(ounit,1010) cput-cpus, myid, vvol, cput-lastcpu, idgetrf, "singular;         "
    case default ;               FATAL( dforce, .true., illegal ifail returned from F07ADF )
    end select

    oBI%mat(1:NN,1:NN) = dMA(1:NN,1:NN)

  endif

1011 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRI ; time=",f10.2,"s ; inverse of Beltrami matrix; idgetrf=",i2," ; ",a18)
1010 format("dforce : ",f10.2," : myid=",i3," ; vvol=",i3," ; called DGETRF ; time=",f10.2,"s ; LU factorization of matrix; idgetrf=",i2," ; ",a18)

end subroutine get_LU_beltrami_matrices



subroutine get_perturbed_solution(vvol, oBI, NN)


  use constants, only :   zero, half, one, two

  use fileunits, only :   ounit

  use cputiming, only :   Tdfp200

  use inputlist, only :   Wmacros, Wdfp200, Lrad, mu, Lconstraint

  use allglobal, only :   ncpu, myid, cpus, MPI_COMM_SPEC, &
                          mn, Iquad, NAdof, &
                          dMA, dMB, dMD, dMG, solution, &
                          dtflux, dpflux, dBdX

  use typedefns

 LOCALS

  INTEGER, intent(in)     :: vvol, NN
  TYPE(MatrixLU),intent(inout) :: oBI

  INTEGER                 :: ideriv, ll, idgetrs
  REAL                    :: dpsi(1:2), work(1:NN+1), rhs(0:NN), dVA(0:NN), dVD(0:NN)
  CHARACTER               :: packorunpack


  ll = Lrad(vvol)  ! Shorthand
  dBdX%L = .true.



  WCALL( dfp200, intghs, ( Iquad(vvol), mn, vvol, ll, 0 ) )


  rhs(0)    = zero
  rhs(1:NN) = -dVA(1:NN)

  if (Lconstraint .eq. 2) then
    work(1:NN+1) = rhs(0:NN)
    call DGETRS('N',NN+1,1,oBI%mat(0,0),NN+1,oBI%ipivot(0:NN),work(1),NN+1,idgetrs) ! Change to DGETRS; 22 Jul 19
    solution(1:NN,-1) = work(2:NN+1)
  else
    solution(1:NN,-1) = rhs(1:NN)
    call DGETRS('N',NN,1,oBI%mat(1,1),NN+1,oBI%ipivot(1:NN),solution(1,-1),NN,idgetrs) ! Change to DGETRS; 22 Jul 19
  endif

  packorunpack = 'U'
  WCALL( dfp200, packab,( packorunpack, vvol, NN,  solution(1:NN,-1), -1 ) ) ! derivatives placed in Ate(lvol,ideriv,1:mn)%s(0:Lrad),

end subroutine get_perturbed_solution







subroutine evaluate_dmupfdx(innout, idof, ii, issym, irz)


    use constants, only :   zero, half, one, two

    use fileunits, only :   ounit

    use numerical, only :   small

    use cputiming, only :   Tdfp200

    use inputlist, only :   Wmacros, Wdfp200, Lrad, mu, Lconstraint, Lcheck, Nvol, Igeometry, mupftol, Lfreebound, dRZ

    use allglobal, only :   ncpu, myid, cpus, MPI_COMM_SPEC, &
                            Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, &
                            Mvol, Iquad, NGdof, &
                            iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                            NAdof, &
                            mn, im, in, mns, &
                            Ate, Aze, Ato, Azo, &     ! only required for debugging;
                            Nt, Nz, dBdX, &
                            dtflux, dpflux, sweight, &
                            Rij, Zij, &
                            diotadxup, dItGpdxtp, dmupfdx, &
                            psifactor, &
                            lmns, &
                            mn, mne, &
                            LocalConstraint, &
                            vvolume, dvolume, &
                            IsMyVolume, IsMyVolumeValue, &
                            Btemn, xoffset, &
                            IPdtdPf, &
                            dMA, dMB, dMG, dMD


  LOCALS:

    INTEGER             ::  vvol, innout, idof, iflag, ii, issym, irz, ll, NN, ifail
    INTEGER             ::  vflag, N, iwork(1:Nvol-1), idgesvx, pvol, order, IDGESV
    INTEGER             ::  iocons
    INTEGER, allocatable::  IPIV(:)
    REAL                ::  det, lfactor, Bt00(1:Mvol, 0:1, -1:2)
    REAL                ::  R(1:Nvol-1), C(1:Nvol-1), work(1:4*Nvol-4), ferr, berr, rcond, tmp(2:Nvol)
    LOGICAL             ::  Lonlysolution, LcomputeDerivatives, dfp100_logical
    REAL, allocatable   ::  dBdmpf(:,:), dBdx2(:)





    vvol = dBdX%vol     ! shorthand

    ll = Lrad(vvol)        ! shorthand
    NN = NAdof(vvol)     ! shorthand;





    lfactor = psifactor(ii,vvol-1+innout)     ! this "pre-conditions" the geometrical degrees-of-freedom;

    if( Lconstraint.eq.1 .or. Lconstraint.eq.3 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) ) then ! will need to accommodate constraints

        if (Lconstraint.eq.3) then

            if( Lfreebound.eq.1 ) then
                order = Mvol
            else
                order = Mvol-1
            endif

            SALLOCATE(dBdmpf , ( 1:order, 1:order ), zero)
            SALLOCATE(dBdx2  , ( 1:order ), zero)
            SALLOCATE(IPIV   , ( 1:order ), zero)


            dmupfdx(1:Nvol,vvol,1,idof,1) = zero ! The helicity multiplier is not varied (constrained). However dtflux is varied in vacuum region, if free boundary

            do pvol = 1, Mvol
                LREGION(pvol)

                do iocons = 0, 1
            if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle
                    WCALL(dfp200, lbpol, (pvol, Bt00(1:Mvol, 0:1, -1:2), 2, iocons)) ! Stores derivative in global variable Btemn
                enddo
            enddo

            dBdmpf(1:order,1:order) = zero ! Initialize. TODO: useless?
            do pvol=1, Mvol-2
                dBdmpf(pvol,   pvol) =  Bt00(pvol+1, 0, 2)
                dBdmpf(pvol+1, pvol) = -Bt00(pvol+1, 1, 2)
            enddo
            dBdmpf(Mvol-1,Mvol-1) = Bt00(Mvol, 0, 2)

            do pvol = vvol, vvol+1
                LREGION(pvol)
                if( pvol.eq.vvol ) then
                    dBdX%innout = 1 ! take derivative w.r.t outer interface
                else !pvol.eq.vvol+1
                    dBdX%innout = 0 ! w.r.t inner interface
                endif

                do iocons = 0, 1
           if( ( Lcoordinatesingularity .and. iocons.eq.0 ) .or. ( Lvacuumregion .and. iocons.eq.1 ) ) cycle
                    WCALL(dfp200, lbpol, (pvol, Bt00(1:Mvol, 0:1, -1:2), -1, iocons)) ! derivate w.r.t geometry
                enddo


            enddo

            dBdx2(1:Mvol-1) = zero
            if( vvol.gt.1 ) then
                dBdx2(vvol-1)   =                 - Bt00(vvol,   0, -1)
            endif
            ;   dBdx2(vvol  )   = Bt00(vvol  , 1, -1) - Bt00(vvol+1, 0, -1)
            if (vvol.lt.Mvol-1) then
            ; dBdx2(vvol+1) = Bt00(vvol+1, 1, -1)
            endif


            if( Lfreebound.eq.1 ) then ! Need to modify last two equations

                dBdmpf(1:Mvol, Mvol  ) = zero

                iocons = 0
                WCALL(dfp200, lbpol, (Mvol, Bt00(1:Mvol, 0:1, -1:2), 1, iocons))

                iflag =  2 ; WCALL( dfp200, curent, ( Mvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,Mvol) ) )

                dBdmpf(Mvol-1, Mvol  ) =  Bt00(Mvol, 0, 1)         !dBdpsit
                dBdmpf(Mvol  , Mvol-1) =  dItGpdxtp( 1, 2, Mvol)    !dIpdpsip
                dBdmpf(Mvol  , Mvol  ) =  dItGpdxtp( 1, 1, Mvol)    !dIpdpsit

                if( vvol.eq.Mvol-1 ) then ! Plasma interface is perturbed
                    iflag = -1 ; WCALL( dfp200, curent, ( Mvol, mn, Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,Mvol) ) )
                    dBdx2( Mvol ) = -dItGpdxtp( 1,-1, Mvol)    !-dIpdxj
                else ! Inner interface is perturbed
                    dBdx2( Mvol ) = zero
                endif
            endif



            call DGESV(order, 1   , dBdmpf(1:order,1:order), order, IPIV, dBdx2(1:order), order, IDGESV )

            dmupfdx(1     , vvol, 2, idof, 1) = zero ! First poloidal flux is always zero
            dmupfdx(2:Mvol, vvol, 2, idof, 1) = lfactor * dBdx2(1:Mvol-1) ! These are the derivatives of pflux

            if( Lfreebound.eq.1 ) then
                dmupfdx(Mvol, vvol, 1, idof, 1) = lfactor * dBdx2(Mvol) ! This is the derivative of tflux
            endif

            DALLOCATE( dBdmpf )
            DALLOCATE( dBdx2  )
            DALLOCATE( IPIV   )


    else ! LocalConstraint

        if( Lconstraint.eq.1 ) then
            iflag = -1 ; WCALL( dfp200, tr00ab, ( vvol, mn, lmns, Nt, Nz, iflag, diotadxup(0:1,-1:2,vvol) ) ) ! compute d(transform)/dx;
        endif

        if( Lvacuumregion .and. Lconstraint.ge.0 ) then
            iflag = -1 ; WCALL( dfp200, curent, ( vvol, mn,       Nt, Nz, iflag, dItGpdxtp(0:1,-1:2,vvol) ) ) ! compute d(Itor,Gpol)/dx;
        endif

        dmupfdx(vvol,1,1,idof,innout) = zero        ! Prepare array
        dmupfdx(vvol,1,2,idof,innout) = zero        ! Prepare array

        if( Lplasmaregion ) then

            if( Lconstraint.eq.1) then

                if( Lcoordinatesingularity ) then ! solution does not depend on dpflux, and only outer transform is a constraint;

                    det = diotadxup(1,1,vvol)
                    FATAL( dfp200, abs(det).lt.small, error computing derivatives of mu          wrt geometry at fixed transform )

                    dmupfdx(vvol,1,1,idof,innout) = - lfactor * (                                                                        diotadxup(1,-1,vvol) ) / det
                    dmupfdx(vvol,1,2,idof,innout) =   zero

                else ! if( .not.Lcoordinatesingularity ) ;

                    det = diotadxup(0,1,vvol) * diotadxup(1,2,vvol) - diotadxup(0,2,vvol) * diotadxup(1,1,vvol)
                    FATAL( dfp200, abs(det).lt.small, error computing derivatives of mu & dpflux wrt geometry at fixed transform )

                    dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + diotadxup(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * diotadxup(1,-1,vvol) ) / det
                    dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - diotadxup(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * diotadxup(1,-1,vvol) ) / det

                endif ! end of if( Lcoordinatesingularity ) ;

            endif ! end of if( Lconstraint.eq.1 ) ;

        else ! Vacuum region

            if    ( Lconstraint.eq.0 ) then ! THIS NEEDS ATTENTION;

                det = dItGpdxtp(0,1,vvol) * dItGpdxtp(1,2,vvol) - dItGpdxtp(0,2,vvol) * dItGpdxtp(1,1,vvol)
                FATAL( dfp200, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )

                dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * dItGpdxtp(0,-1,vvol) - dItGpdxtp(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
                dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * dItGpdxtp(0,-1,vvol) + dItGpdxtp(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det

            else if( Lconstraint.eq.1 ) then

                det = diotadxup(0,1,vvol) * dItGpdxtp(1,2,vvol) - diotadxup(0,2,vvol) * dItGpdxtp(1,1,vvol)
                FATAL( dfp200, abs(det).lt.small, error computing derivatives of dtflux & dpflux wrt geometry at fixed Itor and Gpol )

                dmupfdx(vvol,1,1,idof,innout) = - lfactor * ( + dItGpdxtp(1, 2,vvol) * diotadxup(0,-1,vvol) - diotadxup(0, 2,vvol) * dItGpdxtp(1,-1,vvol) ) / det
                dmupfdx(vvol,1,2,idof,innout) = - lfactor * ( - dItGpdxtp(1, 1,vvol) * diotadxup(0,-1,vvol) + diotadxup(0, 1,vvol) * dItGpdxtp(1,-1,vvol) ) / det

            endif

        endif ! end of if( Lplasmaregion ) ;

    endif ! end of if Lconstraint.eq.3

    endif ! end of if( Lconstraint.eq.1 .or. ( Lvacuumregion .and. Lconstraint.ge.0 ) then;

    vflag = 1 ! this flag instructs volume to continue even if the volume is invalid;
    WCALL( dfp200, volume, ( vvol, vflag ) ) ! compute derivative of volume; wrt to harmonic described by dBdX structure;

end subroutine evaluate_dmupfdx




subroutine evaluate_dBB(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz, LcomputeDerivatives)



  use constants, only : zero, half, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp200, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad

  use cputiming, only : Tdfp200

  use allglobal, only : ncpu, myid, cpus, MPI_COMM_SPEC, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, &
                        Mvol, dpflux, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, &
                        mmpp, &
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ,HdFFdRZ, dBBdmp, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS

LOGICAL, intent(in)     :: LComputeDerivatives
INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)


do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;



    if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) then ! first, determine how B^2 varies with mu and dpflux;

        iflag = 1
        do ideriv=1, 2
            WCALL(dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )
        enddo

        call tfft(    Nt, Nz, dBB(1:Ntz,1), dBB(1:Ntz,2), & ! derivatives of B^2 wrt mu and dpflux;
                    mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )
        idoc = 0
        dBBdmp(idoc+1:idoc+mn  ,lvol,iocons,1) = efmn(1:mn) * BBweight(1:mn) ! pressure;
        dBBdmp(idoc+1:idoc+mn  ,lvol,iocons,2) = cfmn(1:mn) * BBweight(1:mn) ! pressure;
        idoc = idoc + mn   ! even;


        if( NOTstellsym ) then
            dBBdmp(idoc+1:idoc+mn-1,lvol,iocons,1) = ofmn(2:mn) * BBweight(2:mn) ! pressure;
            dBBdmp(idoc+1:idoc+mn-1,lvol,iocons,2) = sfmn(2:mn) * BBweight(2:mn) ! pressure;
            idoc = idoc + mn-1 ! oddd;


        endif ! end of if( NOTstellsym) ;

    endif ! end of if( Lconstraint.eq.1 .OR. Lconstraint.eq.3 ) ;



    ideriv = 0; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )



    ideriv = -1; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)




    dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

    if( Igeometry.ge.3 ) then ! spectral constraints are only required in toroidal or extended-cylindrical geometry;

        if( innout.eq.1 .and. iocons.eq.1 ) then ! include derivatives of spectral constraints;

            if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                if( irz.eq.0 ) then ! take derivative wrt Rbc;
                    dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iRbc(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
                else ! take derivative wrt Zbs;
                    dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iZbs(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
                endif ! end of if( irz.eq.0 ) ;
            else                  ! take derivatives wrt Rbs and Zbc;
                if( irz.eq.0 ) then
                    dII(1:Ntz) = + im(ii) * cosi(1:Ntz,ii) * ( XX(1:Ntz) - MMl * iRij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iRbs(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Rij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * sini(1:Ntz,ii)
                else
                    dII(1:Ntz) = - im(ii) * sini(1:Ntz,ii) * ( YY(1:Ntz) - MMl * iZij(1:Ntz,lvol) ) &
                    - two * ( mmpp(ii) - MMl ) * iZbc(ii,lvol) * ( Rij(1:Ntz,2,0) * iRij(1:Ntz,lvol) + Zij(1:Ntz,2,0) * iZij(1:Ntz,lvol) ) / DDl &
                    + Zij(1:Ntz,2,0) * ( mmpp(ii) - MMl ) * cosi(1:Ntz,ii)
                endif !matches ( irz.eq.0 )
            endif! matches ( issym.eq.0 )

        else

        dII(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;

        endif ! end of if( innout.eq.1 .and. iocons.eq.1 ) ;

        constraint(1:Ntz) = + ( dRij(1:Ntz,lvol) * tRij(1:Ntz,lvol-1+iocons) + dZij(1:Ntz,lvol) * tZij(1:Ntz,lvol-1+iocons) ) / length(1:Ntz)

        if( iocons.eq.0 ) then ! take derivatives of constraints at inner boundary;

            if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,lvol-1) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,lvol-1) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,lvol-1) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,lvol-1) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) &
                                                                        - constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) &
                                                                        - constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else ! if issym.eq.1 ; take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) &
                                                                        - constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,lvol-1)                                              ) / length(1:Ntz) &
                                                                        - constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            endif

        else ! if iocons.eq.1 ; take derivatives of constraints at outer boundary;

            if( innout.eq.0 ) then ! derivative wrt inner boundary coefficient;
                if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tRij(1:Ntz,lvol  )                                              ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tZij(1:Ntz,lvol  )                                              ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                else                  ! take derivatives wrt Rbs and Zbc;
                    if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( - sini(1:Ntz,ii) * tRij(1:Ntz,lvol  )                                              ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    else                ; dLL(1:Ntz) = + ( - cosi(1:Ntz,ii) * tZij(1:Ntz,lvol  )                                              ) / length(1:Ntz) &
                                                                        + constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                    endif
                endif
            else ! if innout.eq.1 ; derivative wrt outer boundary coefficient;
                if( Igeometry.eq.3 .and. lvol.eq.1 ) then ! need to accomodate derivatives of coordinate axis;
                    if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                        if( irz.eq.0 ) then ; dLL(1:Ntz) = ( &   ! d/dRbc ;
                                                            + ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) ) * tRij(1:Ntz,lvol) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) &
                                                            + (                - dZodR(1:Ntz,0,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * ( cosi(1:Ntz,ii) - dRodR(1:Ntz,0,ii) )   &
                                                            + dZij(1:Ntz,lvol) * (                - dZodR(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = ( &   ! d/dZbs ;
                                                            + (                - dRodZ(1:Ntz,1,ii) ) * tRij(1:Ntz,lvol) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) &
                                                            + ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * (                - dRodZ(1:Ntz,1,ii) )   &
                                                            + dZij(1:Ntz,lvol) * ( sini(1:Ntz,ii) - dZodZ(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    else
                        if( irz.eq.0 ) then ; dLL(1:Ntz) =    ( &   ! d/dRbs ;
                                                            + ( sini(1:Ntz,ii) - dRodR(1:Ntz,1,ii) ) * tRij(1:Ntz,lvol) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii)     &
                                                            + (                - dZodR(1:Ntz,1,ii) ) * tZij(1:Ntz,lvol)                                                 &
                                                            - constraint(1:Ntz)                                                                                         &
                                                            * (   dRij(1:Ntz,lvol) * ( sini(1:Ntz,ii)   - dRodR(1:Ntz,1,ii) )                                              &
                                                                + dZij(1:Ntz,lvol) * (                  - dZodR(1:Ntz,1,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        else                ; dLL(1:Ntz) =     ( &   ! d/dZbs ;
                                                            + (                - dRodZ(1:Ntz,0,ii) ) * tRij(1:Ntz,lvol) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) &
                                                            + ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) * tZij(1:Ntz,lvol) &
                                                            - constraint(1:Ntz) &
                                                            * ( dRij(1:Ntz,lvol) * (                - dRodZ(1:Ntz,0,ii) )   &
                                                            + dZij(1:Ntz,lvol) * ( cosi(1:Ntz,ii) - dZodZ(1:Ntz,0,ii) ) ) / length(1:Ntz) ) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    endif
                else
                    if( issym.eq.0 ) then ! take derivatives wrt Rbc and Zbs;
                        if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tRij(1:Ntz,lvol  ) - dRij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) &
                                                                            - constraint(1:Ntz) * dRij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tZij(1:Ntz,lvol  ) + dZij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) &
                                                                            - constraint(1:Ntz) * dZij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        endif
                    else                  ! take derivatives wrt Rbs and Zbc;
                        if( irz.eq.0 ) then ; dLL(1:Ntz) = + ( + sini(1:Ntz,ii) * tRij(1:Ntz,lvol  ) + dRij(1:Ntz,lvol) * im(ii) * cosi(1:Ntz,ii) ) / length(1:Ntz) &
                                                                            - constraint(1:Ntz) * dRij(1:Ntz,lvol) * sini(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        else                ; dLL(1:Ntz) = + ( + cosi(1:Ntz,ii) * tZij(1:Ntz,lvol  ) - dZij(1:Ntz,lvol) * im(ii) * sini(1:Ntz,ii) ) / length(1:Ntz) &
                                                                            - constraint(1:Ntz) * dZij(1:Ntz,lvol) * cosi(1:Ntz,ii) / length(1:Ntz) / length(1:Ntz)
                        endif ! end of if( irz.eq.0 ) ;
                    endif ! end of if( issym.eq.0 ) ;
                endif  ! end of if( Igeometry.eq.3 .and. lvol.eq.1 ) ;
            endif ! end of if( innout.eq.0 ) ;

        endif ! end of if( iocons.eq.0 ) ;

    endif ! end of if( Igeometry.ge.3 ) ;

    call tfft(  Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

    call tfft(  Nt, Nz, dPP(1:Ntz)   , dLL(1:Ntz), &
                mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! evmn and odmn are available as workspace;


    FATAL( dfp200, lvol-1+innout.gt.Mvol, psifactor needs attention )


    idoc = 0

    ;   dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = + efmn(1:mn) * psifactor(ii,lvol-1+innout) * BBweight(1:mn)

    idoc = idoc + mn   ! even;
    if( Igeometry.ge.3 ) then ! Add spectral constraints;
        dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = - sfmn(2:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
                                                             - simn(2:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
      idoc = idoc + mn-1 ! odd;
    endif ! end of if( Igeometry.ge.3) ;

    if( NOTstellsym ) then ! Construct non-stellarator symmetric terms

    ;       dFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = + ofmn(2:mn) * psifactor(ii,lvol-1+innout) * BBweight(2:mn)

        idoc = idoc + mn-1 ! odd;
        if( Igeometry.ge.3 ) then ! Add spectral constraints;
            dFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = - cfmn(1:mn) * psifactor(ii,lvol-1+innout) * epsilon       & ! spectral condensation;
                                                                 - comn(1:mn) * psifactor(ii,lvol-1+innout) * sweight(lvol)   ! poloidal length constraint;
            idoc = idoc + mn   ! even;
        endif ! end of if( Igeometry.ge.3) ;

    endif ! end of if( NOTstellsym) ;

enddo ! end of do iocons;

end subroutine evaluate_dBB


subroutine hessian_dFFdRZ(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)



  use constants, only : zero, half, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp200, Ntor, Igeometry, Lconstraint,&
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad

  use cputiming, only : Tdfp200

  use allglobal, only : ncpu, myid, cpus, &
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, MPI_COMM_SPEC, &
                        Mvol, dpflux, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, &
                        mmpp, &
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ, dBBdmp, HdFFdRZ, &
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS

INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)



do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;


    ideriv = 0; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )


    ideriv = -1; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)*Rij(1:Ntz,0,0)


    dLL(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero ! either no spectral constraint, or not the appropriate interface;
    dII(1:Ntz) = zero


    call tfft(  Nt, Nz, ijreal(1:Ntz), dII(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )


    FATAL( dfp200, lvol-1+innout.gt.Mvol, psifactor needs attention )

    idoc = 0

    ;   HdFFdRZ(idoc+1:idoc+mn    ,iocons,idof,innout,lvol) = + efmn(1:mn) !* psifactor(ii,lvol-1+innout) * BBweight(1:mn)
    idoc = idoc + mn   ! even;

    if( NOTstellsym ) then ! Construct non-stellarator symmetric terms

    ;       HdFFdRZ(idoc+1:idoc+mn-1  ,iocons,idof,innout,lvol) = + ofmn(2:mn) !* psifactor(ii,lvol-1+innout) * BBweight(2:mn)

        idoc = idoc + mn-1 ! odd;

    endif ! end of if( NOTstellsym) ;


enddo ! end of do iocons;

end subroutine hessian_dFFdRZ

subroutine hessian3D_dFFdRZ(lvol, idof, innout, issym, irz, ii, dBB, XX, YY, length, dRR, dZZ, dII, dLL, dPP, Ntz)



  use constants, only : zero, half, one, two

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wmacros, Wdfp200, Ntor, Igeometry, Lconstraint, &
                        gamma, adiabatic, pscale, Lcheck, dRZ, Nvol, &
                        epsilon, Lrad

  use cputiming, only : Tdfp200

  use allglobal, only : ncpu, myid, cpus, YESstellsym, NOTstellsym, MPI_COMM_SPEC,&
                        Lcoordinatesingularity, Lplasmaregion, Lvacuumregion, LocalConstraint, &
                        Mvol, dpflux, &
                        iRbc, iZbs, iRbs, iZbc, & ! Fourier harmonics of geometry; vector of independent variables, position, is "unpacked" into iRbc,iZbs;
                        NOTstellsym, &
                        mn, im, in, mns, &
                        ijreal, &
                        efmn, ofmn, cfmn, sfmn, &
                        evmn, odmn, comn, simn, &
                        Nt, Nz, &
                        cosi, sini, & ! FFT workspace
                        sweight, Lhessian3Dallocated, &
                        mmpp, &
                        LGdof, NGdof, NAdof, &
                        vvolume, dvolume, &
                        Rij, Zij, sg, guvij, iRij, iZij, dRij, dZij, tRij, tZij, & ! Jacobian and metrics; computed in coords;
                        dFFdRZ, dBBdmp, HdFFdRZ, &
                        denergydrr, denergydrz,denergydzr,denergydzz,&
                        BBweight, & ! exponential weight on force-imbalance harmonics;
                        psifactor, &
                        mn, Iquad, &
                        dRodR, dRodZ, dZodR, dZodZ, dBdX, &
                        xoffset

 LOCALS

INTEGER                 :: iocons, lvol, ideriv, id, iflag, Lcurvature, innout, issym, irz, ii, ifail, idoc, idof, Ntz
REAL                    :: lss, DDl, MMl
REAL                    :: dAt(1:Ntz,-1:2), dAz(1:Ntz,-1:2), XX(1:Ntz), YY(1:Ntz), dBB(1:Ntz,-1:2), dII(1:Ntz), dLL(1:Ntz)
REAL                    :: dPP(1:Ntz), length(1:Ntz), dRR(1:Ntz,-1:2), dZZ(1:Ntz,-1:2), constraint(1:Ntz)
REAL                    :: ddFcol1(1:Ntz), ddFcol2(1:Ntz),ddFcol3(1:Ntz), ddFcol4(1:Ntz)



do iocons = 0, 1

    if( lvol.eq.   1 .and. iocons.eq.0 ) cycle ! fixed inner boundary (or coordinate axis); no constraints;
    if( lvol.eq.Mvol .and. iocons.eq.1 ) cycle ! fixed outer boundary                     ; no constraints;


    ideriv = 0; iflag=0

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl,iflag) )



    ideriv = -1; iflag=1

    WCALL( dfp200, lforce, (lvol, iocons, ideriv, Ntz, dBB, XX, YY, length, DDl, MMl, iflag) )

    FATAL( dfp200, vvolume(lvol).lt.small, shall divide by vvolume(lvol)**(gamma+one) )

    ijreal(1:Ntz) = - adiabatic(lvol) * pscale * gamma * dvolume / vvolume(lvol)**(gamma+one) + dBB(1:Ntz,-1)

        if(YESstellsym ) then

            if( irz.eq.0 ) then !derivative wrt Rbs of d/dr(dF)

                 if(innout.eq.iocons) then
                    ddFcol1(1:Ntz) = ijreal(1:Ntz) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) & !Frr       d/dr(dw/dr)
                                   +  (dBB(1:Ntz,0)+adiabatic(lvol) * pscale/vvolume(lvol)**gamma ) * Zij(1:Ntz,2,0) * cosi(1:Ntz,ii)

                    ddFcol2(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) &  !Fzr        d/dr(dw/dz)
                                    + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,2,0) * cosi(1:Ntz,ii) &
                                    + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,0,0) * sini(1:Ntz,ii) * (-im(ii))
                  else

                  ddFcol1(1:Ntz) = ijreal(1:Ntz) * Zij(1:Ntz,2,0) * Rij(1:Ntz,0,0) !Frr       d/dr(dw/dr)
                  ddFcol2(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0) !Fzr        d/dr(dw/dz)
                 endif
            else !derivative wrt to Zbs of d/dz(dF)

                if(innout.eq.iocons) then
                   ddFcol3(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0) &   !Frz      d/dz(dw/dr)
                                  + (dBB(1:Ntz,0) +adiabatic(lvol) * pscale/vvolume(lvol)**gamma )* Rij(1:Ntz,0,0) * im(ii) *cosi(1:Ntz,ii)
                   ddFcol4(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0)  !Fzz    d/dz(dw/dz)
                else
                   ddFcol3(1:Ntz) = ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Zij(1:Ntz,2,0)  !Frz      d/dz(dw/dr)
                    ddFcol4(1:Ntz) =  ijreal(1:Ntz) * Rij(1:Ntz,0,0) * Rij(1:Ntz,2,0)  !Fzz    d/dz(dw/dz)
                endif
           endif
        else


                 ddFcol1(1:Ntz) = zero
                 ddFcol2(1:Ntz) = zero
                 ddFcol3(1:Ntz) = zero
                 ddFcol4(1:Ntz) = zero
        endif


                   call tfft(  Nt, Nz, ddFcol1(1:Ntz),ddFcol2(1:Ntz), & ! recall that ijreal contains derivatives of pressure term;
                               mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ifail )

                   call tfft(  Nt, Nz, ddFcol3(1:Ntz),ddFcol4(1:Ntz), &
                               mn, im(1:mn), in(1:mn), evmn(1:mn), odmn(1:mn), comn(1:mn), simn(1:mn), ifail ) ! evmn and odmn are available as workspace;



           if (irz.eq.0 .and. ii.eq.1) then
               ;idoc=0
                      ;denergydrr(idoc+1:idoc+mn ,lvol,iocons,idof,innout) = -efmn(1:mn)*half !1 wrr
               ;idoc=idoc+mn
                      ;denergydrr(idoc+1:idoc+mn-1,lvol,iocons,idof,innout) = -sfmn(2:mn)*half ! 2 wrz
               ;idoc=idoc+mn-1
           endif
          if (irz.eq.0 .and. ii.gt.1) then
              ;idoc=0
                     ;denergydrr(idoc+1:idoc+mn ,lvol,iocons,idof,innout) = -efmn(1:mn) !1 wrr
              ;idoc=idoc+mn
                     ;denergydrr(idoc+1:idoc+mn-1,lvol,iocons,idof,innout) = -sfmn(2:mn) ! 2 wrz
              ;idoc=idoc+mn-1
            endif

          if(irz.eq.1 .and. ii.gt.1) then
              ;idoc=0
                     ;denergydzr(idoc+1:idoc+mn-1 ,lvol,iocons,idof,innout) =  evmn(1:mn-1) !wzr
              ;idoc=idoc+mn
                      ;denergydzr(idoc+1:idoc+mn-1 ,lvol,iocons,idof,innout) =  simn(2:mn) !wzz
              ;idoc=idoc+mn-1
          endif





    dLL(1:Ntz) = zero !  no spectral constraint, or not the appropriate interface;
    dPP(1:Ntz) = zero !  no spectral constraint, or not the appropriate interface;
    dII(1:Ntz) = zero !  no angle/spectral width constraint


enddo ! end of do iocons;

end subroutine hessian3D_dFFdRZ
