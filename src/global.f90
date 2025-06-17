
module constants

    implicit none

    real(8), parameter :: zero = 0.0
    real(8), parameter :: one = 1.0
    real(8), parameter :: two = 2.0
    real(8), parameter :: three = 3.0
    real(8), parameter :: four = 4.0
    real(8), parameter :: five = 5.0
    real(8), parameter :: six = 6.0
    real(8), parameter :: seven = 7.0
    real(8), parameter :: eight = 8.0
    real(8), parameter :: nine = 9.0
    real(8), parameter :: ten = 10.0

    real(8), parameter :: eleven = 11.0
    real(8), parameter :: twelve = 12.0

    real(8), parameter :: hundred = 100.0
    real(8), parameter :: thousand = 1000.0

    real(8), parameter :: half = one/two
    real(8), parameter :: third = one/three
    real(8), parameter :: quart = one/four
    real(8), parameter :: fifth = one/five
    real(8), parameter :: sixth = one/six

    real(8), parameter :: pi2 = 6.28318530717958623
    real(8), parameter :: pi = pi2/two
    real(8), parameter :: mu0 = 2.0e-07*pi2
    real(8), parameter :: goldenmean = 1.618033988749895

    real(8), parameter :: version = 3.23

end module constants

module numerical

    implicit none

    real(8), parameter :: machprec = 1.11e-16
    real(8), parameter :: vsmall = 100*machprec
    real(8), parameter :: small = 10000*machprec
    real(8), parameter :: sqrtmachprec = sqrt(machprec)
    real(8), parameter :: logtolerance = 1.0e-32

end module numerical

module fileunits

    implicit none

    integer :: iunit = 10
    integer :: ounit = 6
    integer :: gunit = 13

    integer :: aunit = 11
    integer :: dunit = 12
    integer :: hunit = 14
    integer :: munit = 14
    integer :: lunit = 20
    integer :: vunit = 15

contains

end module fileunits

module cputiming

    real(8) :: Tmanual = 0.0, manualT = 0.0
    real(8) :: Trzaxis = 0.0, rzaxisT = 0.0
    real(8) :: Tpackxi = 0.0, packxiT = 0.0
    real(8) :: Tvolume = 0.0, volumeT = 0.0
    real(8) :: Tcoords = 0.0, coordsT = 0.0
    real(8) :: Tbasefn = 0.0, basefnT = 0.0
    real(8) :: Tmemory = 0.0, memoryT = 0.0
    real(8) :: Tmetrix = 0.0, metrixT = 0.0
    real(8) :: Tma00aa = 0.0, ma00aaT = 0.0
    real(8) :: Tmatrix = 0.0, matrixT = 0.0
    real(8) :: Tmp00ac = 0.0, mp00acT = 0.0
    real(8) :: Tma02aa = 0.0, ma02aaT = 0.0
    real(8) :: Tpackab = 0.0, packabT = 0.0
    real(8) :: Ttr00ab = 0.0, tr00abT = 0.0
    real(8) :: Tcurent = 0.0, curentT = 0.0
    real(8) :: Tlforce = 0.0, lforceT = 0.0
    real(8) :: Tforce_real = 0.0, force_realT = 0.0
    real(8) :: Tforce_real_helper = 0.0, force_real_helperT = 0.0
    real(8) :: Tforce_real_dforce = 0.0, force_real_dforceT = 0.0
    real(8) :: Tforce_real_dfp200 = 0.0, force_real_dfp200T = 0.0
    real(8) :: Tforce_real_lforce = 0.0, force_real_lforceT = 0.0
    real(8) :: Tbcast_freal_jac = 0.0, bcast_freal_jacT = 0.0
    real(8) :: Tget_jac = 0.0, get_jacT = 0.0
    real(8) :: Tget_angle_constraint = 0.0, get_angle_constraintT = 0.0
    real(8) :: Tangconst_dfp200 = 0.0, angconst_dfp200T = 0.0
    real(8) :: Tintghs = 0.0, intghsT = 0.0
    real(8) :: Tmtrxhs = 0.0, mtrxhsT = 0.0
    real(8) :: Tlbpol = 0.0, lbpolT = 0.0
    real(8) :: Tbrcast = 0.0, brcastT = 0.0
    real(8) :: Tdfp100 = 0.0, dfp100T = 0.0
    real(8) :: Tdfp200 = 0.0, dfp200T = 0.0
    real(8) :: Tdforce = 0.0, dforceT = 0.0
    real(8) :: Tjo00aa = 0.0, jo00aaT = 0.0
    real(8) :: Tpp00aa = 0.0, pp00aaT = 0.0
    real(8) :: Tpp00ab = 0.0, pp00abT = 0.0
    real(8) :: Tbfield = 0.0, bfieldT = 0.0
    real(8) :: Tstzxyz = 0.0, stzxyzT = 0.0
    real(8) :: Thesian = 0.0, hesianT = 0.0
    real(8) :: Tra00aa = 0.0, ra00aaT = 0.0
    real(8) :: Tnumrec = 0.0, numrecT = 0.0
    real(8) :: Tdcuhre = 0.0, dcuhreT = 0.0
    real(8) :: Tminpack = 0.0, minpackT = 0.0
    real(8) :: Tiqpack = 0.0, iqpackT = 0.0
    real(8) :: Trksuite = 0.0, rksuiteT = 0.0
    real(8) :: Ti1mach = 0.0, i1machT = 0.0
    real(8) :: Td1mach = 0.0, d1machT = 0.0
    real(8) :: Tilut = 0.0, ilutT = 0.0
    real(8) :: Titers = 0.0, itersT = 0.0
    real(8) :: Tsphdf5 = 0.0, sphdf5T = 0.0
    real(8) :: Tpreset = 0.0, presetT = 0.0
    real(8) :: Tglobal = 0.0, globalT = 0.0
    real(8) :: Txspech = 0.0, xspechT = 0.0
    real(8) :: Tinputlist = 0.0, inputlistT = 0.0

    real(8) :: Treadin = 0.0
    real(8) :: Twrtend = 0.0

end module cputiming

module typedefns

    type subgrid
        real(8), allocatable :: s(:)
        integer, allocatable :: i(:)
    end type subgrid

    type MatrixLU
        real(8), allocatable :: mat(:, :)
        integer, allocatable :: ipivot(:)
    end type MatrixLU

    type derivative
        logical :: L
        integer :: vol
        integer :: innout
        integer :: ii
        integer :: irz
        integer :: issym
    end type derivative

end module typedefns

module allglobal

    use constants
    use typedefns

    implicit none

    integer :: myid
    integer :: ncpu
    integer :: IsMyVolumeValue
    real(8) :: cpus
    integer :: MPI_COMM_SPEC

    logical :: skip_write = .false.

    real(8) :: pi2nfp
    real(8) :: pi2pi2nfp
    real(8) :: pi2pi2nfphalf
    real(8) :: pi2pi2nfpquart

    character(LEN=1000) :: ext

    real(8) :: ForceErr
    real(8) :: Energy
    real(8) :: BnsErr

    real(8), allocatable :: IPDt(:), IPDtDpf(:, :)

    integer :: Mvol

    logical :: YESstellsym
    logical :: NOTstellsym

    real(8), allocatable :: cheby(:, :)
    real(8), allocatable :: zernike(:, :, :)

    real(8), allocatable :: TT(:, :, :)
    real(8), allocatable :: RTT(:, :, :, :)

    real(8), allocatable :: RTM(:, :)
    real(8), allocatable :: ZernikeDof(:)

    integer :: mne
    integer, allocatable :: ime(:)
    integer, allocatable :: ine(:)

    integer :: mns
    integer, allocatable :: ims(:)
    integer, allocatable :: ins(:)

    integer :: lMpol
    integer :: lNtor
    integer :: sMpol
    integer :: sNtor

    real(8) :: xoffset = 1.0

    logical, allocatable :: ImagneticOK(:)

    logical :: IconstraintOK

    real(8), allocatable :: beltramierror(:, :)

    integer :: mn
    integer, allocatable :: im(:)
    integer, allocatable :: in(:)

    real(8), allocatable :: halfmm(:)
    real(8), allocatable :: regumm(:)

    real(8) :: Rscale
    real(8), allocatable :: psifactor(:, :)
    real(8), allocatable :: inifactor(:, :)

    real(8), allocatable :: BBweight(:)

    real(8), allocatable :: mmpp(:)

    real(8), allocatable :: iRbc(:, :)
    real(8), allocatable :: iZbs(:, :)
    real(8), allocatable :: iRbs(:, :)
    real(8), allocatable :: iZbc(:, :)

    real(8), allocatable :: dRbc(:, :)
    real(8), allocatable :: dZbs(:, :)
    real(8), allocatable :: dRbs(:, :)
    real(8), allocatable :: dZbc(:, :)

    real(8), allocatable :: iRij(:, :)
    real(8), allocatable :: iZij(:, :)
    real(8), allocatable :: dRij(:, :)
    real(8), allocatable :: dZij(:, :)
    real(8), allocatable :: tRij(:, :)
    real(8), allocatable :: tZij(:, :)

    real(8), allocatable :: iVns(:)
    real(8), allocatable :: iBns(:)
    real(8), allocatable :: iVnc(:)
    real(8), allocatable :: iBnc(:)

    real(8), allocatable :: lRbc(:)
    real(8), allocatable :: lZbs(:)
    real(8), allocatable :: lRbs(:)
    real(8), allocatable :: lZbc(:)

    integer :: num_modes
    integer, allocatable :: mmRZRZ(:), nnRZRZ(:)
    real(8), allocatable :: allRZRZ(:, :, :)

    integer :: Nt
    integer :: Nz
    integer :: Ntz
    integer :: hNt
    integer :: hNz
    real(8) :: soNtz

    real(8), allocatable :: Rij(:, :, :)
    real(8), allocatable :: Zij(:, :, :)
    real(8), allocatable :: Xij(:, :, :)
    real(8), allocatable :: Yij(:, :, :)
    real(8), allocatable :: sg(:, :)
    real(8), allocatable :: guvij(:, :, :, :)
    real(8), allocatable :: gvuij(:, :, :)
    real(8), allocatable :: guvijsave(:, :, :, :)

    integer, allocatable :: ki(:, :)
    integer, allocatable :: kijs(:, :, :)
    integer, allocatable :: kija(:, :, :)

    integer, allocatable :: iotakkii(:)
    integer, allocatable :: iotaksub(:, :)
    integer, allocatable :: iotakadd(:, :)
    integer, allocatable :: iotaksgn(:, :)

    real(8), allocatable :: efmn(:)
    real(8), allocatable :: ofmn(:)
    real(8), allocatable :: cfmn(:)
    real(8), allocatable :: sfmn(:)
    real(8), allocatable :: evmn(:)
    real(8), allocatable :: odmn(:)
    real(8), allocatable :: comn(:)
    real(8), allocatable :: simn(:)

    real(8), allocatable :: ijreal(:)
    real(8), allocatable :: ijimag(:)
    real(8), allocatable :: jireal(:)
    real(8), allocatable :: jiimag(:)

    real(8), allocatable :: jkreal(:)
    real(8), allocatable :: jkimag(:)
    real(8), allocatable :: kjreal(:)
    real(8), allocatable :: kjimag(:)

    real(8), allocatable :: Bsupumn(:, :, :)
    real(8), allocatable :: Bsupvmn(:, :, :)

    real(8), allocatable :: goomne(:, :)
    real(8), allocatable :: goomno(:, :)
    real(8), allocatable :: gssmne(:, :)
    real(8), allocatable :: gssmno(:, :)
    real(8), allocatable :: gstmne(:, :)
    real(8), allocatable :: gstmno(:, :)
    real(8), allocatable :: gszmne(:, :)
    real(8), allocatable :: gszmno(:, :)
    real(8), allocatable :: gttmne(:, :)
    real(8), allocatable :: gttmno(:, :)
    real(8), allocatable :: gtzmne(:, :)
    real(8), allocatable :: gtzmno(:, :)
    real(8), allocatable :: gzzmne(:, :)
    real(8), allocatable :: gzzmno(:, :)

    real(8), allocatable :: DToocc(:, :, :, :)
    real(8), allocatable :: DToocs(:, :, :, :)
    real(8), allocatable :: DToosc(:, :, :, :)
    real(8), allocatable :: DTooss(:, :, :, :)
    real(8), allocatable :: TTsscc(:, :, :, :)
    real(8), allocatable :: TTsscs(:, :, :, :)
    real(8), allocatable :: TTsssc(:, :, :, :)
    real(8), allocatable :: TTssss(:, :, :, :)
    real(8), allocatable :: TDstcc(:, :, :, :)
    real(8), allocatable :: TDstcs(:, :, :, :)
    real(8), allocatable :: TDstsc(:, :, :, :)
    real(8), allocatable :: TDstss(:, :, :, :)
    real(8), allocatable :: TDszcc(:, :, :, :)
    real(8), allocatable :: TDszcs(:, :, :, :)
    real(8), allocatable :: TDszsc(:, :, :, :)
    real(8), allocatable :: TDszss(:, :, :, :)
    real(8), allocatable :: DDttcc(:, :, :, :)
    real(8), allocatable :: DDttcs(:, :, :, :)
    real(8), allocatable :: DDttsc(:, :, :, :)
    real(8), allocatable :: DDttss(:, :, :, :)
    real(8), allocatable :: DDtzcc(:, :, :, :)
    real(8), allocatable :: DDtzcs(:, :, :, :)
    real(8), allocatable :: DDtzsc(:, :, :, :)
    real(8), allocatable :: DDtzss(:, :, :, :)
    real(8), allocatable :: DDzzcc(:, :, :, :)
    real(8), allocatable :: DDzzcs(:, :, :, :)
    real(8), allocatable :: DDzzsc(:, :, :, :)
    real(8), allocatable :: DDzzss(:, :, :, :)

    real(8), allocatable :: Tsc(:, :)
    real(8), allocatable :: Tss(:, :)
    real(8), allocatable :: Dtc(:, :)
    real(8), allocatable :: Dts(:, :)
    real(8), allocatable :: Dzc(:, :)
    real(8), allocatable :: Dzs(:, :)
    real(8), allocatable :: Ttc(:, :)
    real(8), allocatable :: Tzc(:, :)
    real(8), allocatable :: Tts(:, :)
    real(8), allocatable :: Tzs(:, :)

    real(8), allocatable :: dtflux(:)
    real(8), allocatable :: dpflux(:)

    real(8), allocatable :: sweight(:)

    integer, allocatable :: NAdof(:)
    integer, allocatable :: Nfielddof(:)

    type(subgrid), allocatable :: Ate(:, :, :)
    type(subgrid), allocatable :: Aze(:, :, :)
    type(subgrid), allocatable :: Ato(:, :, :)
    type(subgrid), allocatable :: Azo(:, :, :)

    integer, allocatable :: Lma(:, :)
    integer, allocatable :: Lmb(:, :)
    integer, allocatable :: Lmc(:, :)
    integer, allocatable :: Lmd(:, :)
    integer, allocatable :: Lme(:, :)
    integer, allocatable :: Lmf(:, :)
    integer, allocatable :: Lmg(:, :)
    integer, allocatable :: Lmh(:, :)

    real(8), allocatable :: Lmavalue(:, :)
    real(8), allocatable :: Lmbvalue(:, :)
    real(8), allocatable :: Lmcvalue(:, :)
    real(8), allocatable :: Lmdvalue(:, :)
    real(8), allocatable :: Lmevalue(:, :)
    real(8), allocatable :: Lmfvalue(:, :)
    real(8), allocatable :: Lmgvalue(:, :)
    real(8), allocatable :: Lmhvalue(:, :)

    integer, allocatable :: Fso(:, :)
    integer, allocatable :: Fse(:, :)

    logical :: Lcoordinatesingularity
    logical :: Lplasmaregion = .true.
    logical :: Lvacuumregion = .false.
    logical :: Lsavedguvij
    logical :: Localconstraint

    real(8), allocatable :: Remn_ext(:, :, :)
    real(8), allocatable :: Romn_ext(:, :, :)
    real(8), allocatable :: Zemn_ext(:, :, :)
    real(8), allocatable :: Zomn_ext(:, :, :)

    logical :: use_ext_mesh = .false.

    real(8), allocatable :: dMA(:, :)
    real(8), allocatable :: dMB(:, :)
    real(8), allocatable :: dMD(:, :)

    real(8), allocatable :: dMAS(:)
    real(8), allocatable :: dMDS(:)
    integer, allocatable :: idMAS(:)
    integer, allocatable :: jdMAS(:)
    integer, allocatable :: NdMASmax(:)
    integer, allocatable :: NdMAS(:)

    real(8), allocatable :: dMG(:)

    real(8), allocatable :: AdotX(:)
    real(8), allocatable :: DdotX(:)

    real(8), allocatable :: solution(:, :)

    real(8), allocatable :: GMRESlastsolution(:, :, :)

    real(8), allocatable :: MBpsi(:)

    real(8), allocatable :: BeltramiInverse(:, :)

    real(8), allocatable :: diotadxup(:, :, :)
    real(8), allocatable :: dItGpdxtp(:, :, :)

    real(8), allocatable :: glambda(:, :, :, :)

    integer :: lmns

    real(8), allocatable :: dlambdaout(:, :, :)

    real(8), allocatable :: Bemn(:, :, :)
    real(8), allocatable :: Iomn(:, :)
    real(8), allocatable :: Somn(:, :, :)
    real(8), allocatable :: Pomn(:, :, :)

    real(8), allocatable :: Bomn(:, :, :)
    real(8), allocatable :: Iemn(:, :)
    real(8), allocatable :: Semn(:, :, :)
    real(8), allocatable :: Pemn(:, :, :)

    real(8), allocatable :: BBe(:)
    real(8), allocatable :: IIo(:)
    real(8), allocatable :: BBo(:)
    real(8), allocatable :: IIe(:)

    real(8), allocatable :: Btemn(:, :, :)
    real(8), allocatable :: Bzemn(:, :, :)
    real(8), allocatable :: Btomn(:, :, :)
    real(8), allocatable :: Bzomn(:, :, :)

    real(8), allocatable :: Bloweremn(:, :)
    real(8), allocatable :: Bloweromn(:, :)

    integer :: LGdof
    integer :: NGdof

    real(8), allocatable :: dBBdRZ(:, :, :)
    real(8), allocatable :: dIIdRZ(:, :)

    real(8), allocatable :: dFFdRZ(:, :, :, :, :)
    real(8), allocatable :: dBBdmp(:, :, :, :)

    real(8), allocatable :: freal_jac(:, :, :, :, :)
    real(8), allocatable :: freal_dBBdmp(:, :, :, :)
    integer :: dbbdmp_in_frealjac = 1

    real(8), allocatable :: HdFFdRZ(:, :, :, :, :)

    real(8), allocatable :: denergydrr(:, :, :, :, :)
    real(8), allocatable :: denergydrz(:, :, :, :, :)
    real(8), allocatable :: denergydzr(:, :, :, :, :)
    real(8), allocatable :: denergydzz(:, :, :, :, :)

    real(8), allocatable :: dmupfdx(:, :, :, :, :)

    real(8), allocatable :: freal_jac_full(:, :)
    real(8), allocatable :: dessian(:, :)

    real(8), allocatable :: force_final(:)

    real(8), allocatable :: cosi(:, :)
    real(8), allocatable :: sini(:, :)
    real(8), allocatable :: gteta(:)
    real(8), allocatable :: gzeta(:)

    real(8), allocatable :: ajk(:)

    real(8), allocatable :: dRadR(:, :, :, :)
    real(8), allocatable :: dRadZ(:, :, :, :)
    real(8), allocatable :: dZadR(:, :, :, :)
    real(8), allocatable :: dZadZ(:, :, :, :)

    real(8), allocatable :: dRodR(:, :, :)
    real(8), allocatable :: dRodZ(:, :, :)
    real(8), allocatable :: dZodR(:, :, :)
    real(8), allocatable :: dZodZ(:, :, :)

    integer, allocatable :: djkp(:, :)
    integer, allocatable :: djkm(:, :)

    real(8), allocatable :: lBBintegral(:)
    real(8), allocatable :: lABintegral(:)

    real(8), allocatable :: vvolume(:)
    real(8) :: dvolume

    integer :: ivol

    real(8) :: gBzeta

    integer, allocatable :: Iquad(:)

    real(8), allocatable :: gaussianweight(:, :)
    real(8), allocatable :: gaussianabscissae(:, :)

    logical :: LBlinear
    logical :: LBnewton
    logical :: LBsequad

    real(8) :: oRZp(1:3)

    type(derivative) :: dBdX

    integer :: globaljk
    real(8), allocatable :: Dxyz(:, :)
    real(8), allocatable :: Nxyz(:, :)
    real(8), allocatable :: Jxyz(:, :)

    real(8) :: tetazeta(1:2)

    real(8) :: virtualcasingfactor = -one/(four*pi)

    integer :: IBerror

    integer :: nfreeboundaryiterations

    integer, parameter :: Node = 2

    logical :: first_free_bound = .false.

contains

    subroutine build_vector_potential(lvol, iocons, aderiv, tderiv)

        use constants, only: zero, half

        use fileunits, only: ounit

        use inputlist, only: Lrad, Wbuild_vector_potential, Wmacros

        use cputiming

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        integer :: aderiv
        integer :: tderiv
        integer :: ii, &
                   ll, &
                   mi, &
                   lvol, &
                   iocons
        real(8) :: mfactor

        efmn(1:mn) = zero; sfmn(1:mn) = zero; cfmn(1:mn) = zero; ofmn(1:mn) = zero

        do ii = 1, mn

            if (Lcoordinatesingularity) then
                mi = im(ii)
                do ll = mi, Lrad(lvol), 2
                    ; ; efmn(ii) = efmn(ii) + Ate(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    ; ; cfmn(ii) = cfmn(ii) + Aze(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    if (NOTstellsym) then; ofmn(ii) = ofmn(ii) + Ato(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                        ; ; sfmn(ii) = sfmn(ii) + Azo(lvol, aderiv, ii)%s(ll)*RTT(ll, mi, iocons, 1)*half
                    end if
                end do
            else
                do ll = 0, Lrad(lvol)
                    ; ; efmn(ii) = efmn(ii) + Ate(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                    ; ; cfmn(ii) = cfmn(ii) + Aze(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                    if (NOTstellsym) then; ofmn(ii) = ofmn(ii) + Ato(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                        ; ; sfmn(ii) = sfmn(ii) + Azo(lvol, aderiv, ii)%s(ll)*TT(ll, iocons, 1)
                    end if
                end do
            end if
        end do

    end subroutine build_vector_potential

    subroutine set_mpi_comm(comm)

        implicit none

        integer, intent(in) :: comm
        integer :: ierr

        MPI_COMM_SPEC = comm

        myid = 0; ncpu = 1

        ierr = 0
        call MPI_COMM_RANK(MPI_COMM_SPEC, myid, ierr)
        if (ierr /= 0) write (*, *) "error in call to MPI_COMM_RANK"

        ierr = 0
        call MPI_COMM_SIZE(MPI_COMM_SPEC, ncpu, ierr)
        if (ierr /= 0) write (*, *) "error in call to MPI_COMM_SIZE"

    end subroutine

    subroutine read_inputlists_from_file()

        use constants
        use fileunits
        use inputlist

#ifdef IFORT
        use ifport
#endif

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        logical :: Lspexist
        integer :: filepos, seek_status, cpfile, instat, idx_mode

        character(len=1000) :: line

        integer :: mm, nn, MNMAX
        real(8), allocatable :: RZRZ(:, :)

        inquire (file=trim(ext)//".sp", exist=Lspexist)
        if (.not. Lspexist) then
            write (6, '("readin :      fatal : myid=",i3," ; .not.Lspexist ; the input file does not exist;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : .not.Lspexist : the input file does not exist ;"
        end if

        open (iunit, file=trim(ext)//".sp", status="old")

        instat = 0

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput - cpus
        end if

        read (iunit, physicslist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in physicslist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput - cpus
        end if

        read (iunit, numericlist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in numericlist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput - cpus
        end if

        read (iunit, locallist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in locallist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput - cpus
        end if

        read (iunit, globallist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in globallist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput - cpus
        end if

        read (iunit, diagnosticslist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in diagnosticslist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput - cpus
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput - cpus
        end if

        read (iunit, screenlist, iostat=instat)
        if (instat /= 0) then
            backspace (iunit)
            read (iunit, fmt='(A)') line
            write (*, '(A)') 'Invalid line in screenlist: '//trim(line)
        end if

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput - cpus
        end if

        instat = 0

        num_modes = 0

        MNMAX = MNtor + 1 + MMpol*(2*MNtor + 1)
        if (allocated(mmRZRZ)) deallocate (mmRZRZ, nnRZRZ, allRZRZ)
        allocate (mmRZRZ(1:MNMAX), nnRZRZ(1:MNMAX), allRZRZ(1:4, 1:Nvol, 1:MNMAX))

        if (Linitialize <= 0) then

            if (Nvol < 1 .or. Nvol > MNvol) then
                write (6, '("readin :      fatal : myid=",i3," ; Nvol.lt.1 .or. Nvol.gt.MNvol ; invalid Nvol: may need to recompile with higher MNvol;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "readin : Nvol.lt.1 .or. Nvol.gt.MNvol : invalid Nvol: may need to recompile with higher MNvol ;"
            end if

            if (allocated(RZRZ)) deallocate (RZRZ)
            allocate (RZRZ(1:4, 1:Nvol), stat=astat)
            RZRZ(1:4, 1:Nvol) = zero

#ifdef IFORT
            filepos = ftell(iunit) + 1
#else
            call ftell(iunit, filepos)
#endif
            do
                read (iunit, *, iostat=instat) mm, nn, RZRZ(1:4, 1:Nvol)
                if (instat /= 0) exit

                num_modes = num_modes + 1
            end do

            rewind (iunit)

#ifdef IFORT
            seek_status = fseek(iunit, filepos, 0)
#else
            call fseek(iunit, filepos, 0, seek_status)
#endif
            if (seek_status /= 0) then
                write (6, '("inplst :      fatal : myid=",i3," ; seek_status.ne.0 ; failed to seek to end of input namelists;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "inplst : seek_status.ne.0 : failed to seek to end of input namelists ;"
            end if

            do idx_mode = 1, num_modes
                read (iunit, *, iostat=instat) mmRZRZ(idx_mode), nnRZRZ(idx_mode), allRZRZ(1:4, 1:Nvol, idx_mode)
            end do

            deallocate (RZRZ, stat=astat)

        end if

        close (iunit)

    end subroutine

    subroutine write_spec_namelist()
        use constants
        use fileunits
        use inputlist

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        logical :: exist
        character(LEN=100), parameter :: example = 'example.sp'

        if (myid == 0) then
            inquire (file=trim(example), EXIST=exist)
            if (exist) then
                write (6, '("global :      fatal : myid=",i3," ; exist ; example input file example.sp already existed;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "global : exist : example input file example.sp already existed ;"
            end if
            open (iunit, file=trim(example), status='unknown', action='write')
            write (iunit, physicslist)
            write (iunit, numericlist)
            write (iunit, locallist)
            write (iunit, globallist)
            write (iunit, diagnosticslist)
            write (iunit, screenlist)
            close (iunit)
        end if

        return
    end subroutine

    subroutine check_inputs()

        use numerical
        use constants
        use fileunits
        use inputlist
        use cputiming, only: Treadin

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        integer :: vvol
        real(8) :: xx, toroidalflux, toroidalcurrent

        Mvol = Nvol + Lfreebound

        write (ounit, '("readin : ", 10x ," : ")')

        cput = MPI_WTIME()

        write (ounit, 1010) cput - cpus, Igeometry, Istellsym, Lreflect
        write (ounit, 1011) Lfreebound, phiedge, curtor, curpol
        write (ounit, 1012) gamma
        write (ounit, 1013) Nfp, Nvol, Mvol, Mpol, Ntor
        write (ounit, 1014) pscale, Ladiabatic, Lconstraint, mupftol, mupfits
        write (ounit, 1015) Lrad(1:min(Mvol, 32))

1010    format("readin : ", f10.2, " : Igeometry=", i3, " ; Istellsym=", i3, " ; Lreflect="i3" ;")
1011    format("readin : ", 10x, " : Lfreebound=", i3, " ; phiedge="es23.15" ; curtor="es23.15" ; curpol="es23.15" ;")
1012    format("readin : ", 10x, " : gamma="es23.15" ;")
1013    format("readin : ", 10x, " : Nfp=", i3, " ; Nvol=", i3, " ; Mvol=", i3, " ; Mpol=", i3, " ; Ntor=", i3, " ;")
1014    format("readin : ", 10x, " : pscale="es13.5" ; Ladiabatic="i2" ; Lconstraint="i3" ; mupf: tol,its="es10.2" ,"i4" ;")
1015    format("readin : ", 10x, " : Lrad = "257(i2, ",", :))

        if (Igeometry < 1 .or. Igeometry > 3) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.lt.1 .or. Igeometry.gt.3 ; invalid geometry;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.lt.1 .or. Igeometry.gt.3 : invalid geometry ;"
        end if
        if (Nfp <= 0) then
            write (6, '("readin :      fatal : myid=",i3," ; Nfp.le.0 ; invalid Nfp;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Nfp.le.0 : invalid Nfp ;"
        end if
        if (Mpol < 0 .or. Mpol > MMpol) then
            write (6, '("readin :      fatal : myid=",i3," ; Mpol.lt.0 .or. Mpol.gt.MMpol ; invalid poloidal resolution: may need to recompile with higher MMpol;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Mpol.lt.0 .or. Mpol.gt.MMpol : invalid poloidal resolution: may need to recompile with higher MMpol ;"
        end if
        if (Ntor < 0 .or. Ntor > MNtor) then
            write (6, '("readin :      fatal : myid=",i3," ; Ntor.lt.0 .or. Ntor.gt.MNtor ; invalid toroidal resolution: may need to recompile with higher MNtor;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Ntor.lt.0 .or. Ntor.gt.MNtor : invalid toroidal resolution: may need to recompile with higher MNtor ;"
        end if
        if (Nvol < 1 .or. Nvol > MNvol) then
            write (6, '("readin :      fatal : myid=",i3," ; Nvol.lt.1 .or. Nvol.gt.MNvol ; invalid Nvol: may need to recompile with higher MNvol;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Nvol.lt.1 .or. Nvol.gt.MNvol : invalid Nvol: may need to recompile with higher MNvol ;"
        end if
        if (mupftol <= zero) then
            write (6, '("readin :      fatal : myid=",i3," ; mupftol.le.zero ; mupftol is too small;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : mupftol.le.zero : mupftol is too small ;"
        end if
        if (abs(one + gamma) < vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(one+gamma).lt.vsmall ; 1+gamma appears in denominator in dforce;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(one+gamma).lt.vsmall : 1+gamma appears in denominator in dforce ;"
        end if
        if (abs(one - gamma) < vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(one-gamma).lt.vsmall ; 1-gamma appears in denominator in fu00aa;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(one-gamma).lt.vsmall : 1-gamma appears in denominator in fu00aa ;"
        end if
        if (Lconstraint < -1 .or. Lconstraint > 3) then
            write (6, '("readin :      fatal : myid=",i3," ; Lconstraint.lt.-1 .or. Lconstraint.gt.3 ; illegal Lconstraint;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Lconstraint.lt.-1 .or. Lconstraint.gt.3 : illegal Lconstraint ;"
        end if
        if (Igeometry == 1 .and. rpol < vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.1 .and. rpol.lt.vsmall ; poloidal extent of slab too small or negative;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.1 .and. rpol.lt.vsmall : poloidal extent of slab too small or negative ;"
        end if
        if (Igeometry == 1 .and. rtor < vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.1 .and. rtor.lt.vsmall ; toroidal extent of slab too small or negative;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.1 .and. rtor.lt.vsmall : toroidal extent of slab too small or negative ;"
        end if

        if (Istellsym == 1) then
            Rbs(-MNtor:MNtor, -MMpol:MMpol) = zero
            Zbc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Rws(-MNtor:MNtor, -MMpol:MMpol) = zero
            Zwc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Vnc(-MNtor:MNtor, -MMpol:MMpol) = zero
            Bnc(-MNtor:MNtor, -MMpol:MMpol) = zero
        end if

        if (abs(tflux(Nvol)) < vsmall) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(tflux(Nvol)).lt. vsmall ; enclosed toroidal flux cannot be zero;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(tflux(Nvol)).lt. vsmall : enclosed toroidal flux cannot be zero ;"
        end if

        toroidalflux = tflux(Nvol)

        tflux(1:Mvol) = tflux(1:Mvol)/toroidalflux
        pflux(1:Mvol) = pflux(1:Mvol)/toroidalflux

        if (tflux(1) < zero) then
            write (6, '("readin :      fatal : myid=",i3," ; tflux(1).lt.zero ; enclosed toroidal flux cannot be zero;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : tflux(1).lt.zero : enclosed toroidal flux cannot be zero ;"
        end if
        do vvol = 2, Mvol
        end do

        do vvol = 1, Mvol
            if (Lrad(vvol) < 2) then
                write (6, '("readin :      fatal : myid=",i3," ; Lrad(vvol ).lt.2 ; require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied;")') myid
                call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
                stop "readin : Lrad(vvol ).lt.2 : require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied ;"
            end if
        end do

        if (Igeometry >= 2 .and. Lrad(1) < Mpol) then
            write (ounit, '("readin : ",f10.2," : Minimum Lrad(1) is Mpol, automatically adjusted it to Mpol+4")') cput - cpus
            Lrad(1) = Mpol + 4
        end if
        if (mupfits <= 0) then
            write (6, '("readin :      fatal : myid=",i3," ; mupfits.le.0 ; must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : mupfits.le.0 : must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input ;"
        end if

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1020) cput - cpus, Linitialize, LautoinitBn, Lzerovac, Ndiscrete
        write (ounit, 1021) Nquad, iMpol, iNtor
        write (ounit, 1022) Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
        write (ounit, 1023) Lextrap, Mregular, Lrzaxis, Ntoraxis

1020    format("readin : ", f10.2, " : Linitialize=", i3, " ;LautoinitBn=", i3, " ; Lzerovac=", i2, " ; Ndiscrete="i2" ;")
1021    format("readin : ", 10x, " : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022    format("readin : ", 10x, " : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023    format("readin : ", 10x, " : Lextrap="i2" ; Mregular="i3" ; Lrzaxis="i2" ; Ntoraxis="i2" ;")

        if (Ndiscrete <= 0) then
            write (6, '("readin :      fatal : myid=",i3," ; Ndiscrete.le.0 ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Ndiscrete.le.0 : error ;"
        end if

        if (iotatol > one) then
            write (6, '("readin :      fatal : myid=",i3," ; iotatol.gt.one ; illegal value for sparse tolerance;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : iotatol.gt.one : illegal value for sparse tolerance ;"
        end if

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1030) cput - cpus, LBeltrami, Linitgues, Lmatsolver, LGMRESprec, NiterGMRES, epsGMRES, epsILU

1030    format("readin : ", f10.2, " : LBeltrami="i2" ; Linitgues="i2" ; Lmatsolver="i2" ; LGMRESprec="i2" ; NiterGMRES="i4" ; epsGMRES="es13.5" ; epsILU="es13.5" ;")

        if (LBeltrami < 0 .or. LBeltrami > 7) then
            write (6, '("readin :      fatal : myid=",i3," ; LBeltrami.lt.0 .or. LBeltrami.gt.7 ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : LBeltrami.lt.0 .or. LBeltrami.gt.7 : error ;"
        end if
        if (LGMRESprec < 0 .or. LGMRESprec > 1) then
            write (6, '("readin :      fatal : myid=",i3," ; LGMRESprec.lt.0 .or. LGMRESprec.gt.1 ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : LGMRESprec.lt.0 .or. LGMRESprec.gt.1 : error ;"
        end if
        if (NiterGMRES < 0) then
            write (6, '("readin :      fatal : myid=",i3," ; NiterGMRES.lt.0 ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : NiterGMRES.lt.0 : error ;"
        end if
        if (abs(epsGMRES) <= machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(epsGMRES).le.machprec ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(epsGMRES).le.machprec : error ;"
        end if
        if (abs(epsILU) <= machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(epsILU).le.machprec ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(epsILU).le.machprec : error ;"
        end if

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1040) cput - cpus, Lfindzero
        write (ounit, 1041) escale, opsilon, pcondense, epsilon, wpoloidal, upsilon
        write (ounit, 1042) forcetol, c05xmax, c05xtol, c05factor, LreadGF
        write (ounit, 1043) mfreeits, gBntol, gBnbld
        write (ounit, 1044) vcasingeps, vcasingtol, vcasingits, vcasingper

1040    format("readin : ", f10.2, " : Lfindzero="i2" ;")
1041    format("readin : ", 10x, " : escale="es13.5" ; opsilon="es13.5" ; pcondense="f7.3" ; epsilon="es13.5" ; wpoloidal="f7.4" ; upsilon="es13.5" ;")
1042    format("readin : ", 10x, " : forcetol="es13.5" ; c05xmax="es13.5" ; c05xtol="es13.5" ; c05factor="es13.5" ; LreadGF="L2" ; ")
1043    format("readin : ", 10x, " : mfreeits="i4" ; gBntol="es13.5" ; gBnbld="es13.5" ;")
1044    format("readin : ", 10x, " : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i6" ; vcasingper="i6" ;")

        if (escale < zero) then
            write (6, '("readin :      fatal : myid=",i3," ; escale      .lt.zero ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : escale      .lt.zero : error ;"
        end if
        if (pcondense < one) then
            write (6, '("readin :      fatal : myid=",i3," ; pcondense   .lt.one ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : pcondense   .lt.one : error ;"
        end if
        if (abs(c05xtol) <= machprec) then
            write (6, '("readin :      fatal : myid=",i3," ; abs(c05xtol).le.machprec ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : abs(c05xtol).le.machprec : error ;"
        end if
        if (c05factor <= zero) then
            write (6, '("readin :      fatal : myid=",i3," ; c05factor   .le.zero ; error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : c05factor   .le.zero : error ;"
        end if

        if (Igeometry == 3 .and. pcondense <= zero) then
            write (6, '("readin :      fatal : myid=",i3," ; Igeometry.eq.3 .and. pcondense.le.zero ; pcondense must be positive;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : Igeometry.eq.3 .and. pcondense.le.zero : pcondense must be positive ;"
        end if

        write (ounit, '("readin : ", 10x ," : ")')

        write (ounit, 1050) cput - cpus, odetol, nPpts
        write (ounit, 1051) LHevalues, LHevectors, LHmatrix, Lperturbed, dpp, dqq, dRZ, Lcheck, Ltiming

1050    format("readin : ", f10.2, " : odetol="es10.2" ; nPpts="i6" ;")
1051    format("readin : ", 10x, " : LHevalues="L2" ; LHevectors="L2" ; LHmatrix="L2" ; Lperturbed="i2" ; dpp="i3" ; dqq="i3" ; dRZ="es16.8" ; Lcheck="i3" ; Ltiming="L2" ;")

        if (odetol <= zero) then
            write (6, '("readin :      fatal : myid=",i3," ; odetol.le.zero ; input error;")') myid
            call MPI_ABORT(MPI_COMM_SPEC, 1, ierr)
            stop "readin : odetol.le.zero : input error ;"
        end if

        write (ounit, '("readin : ", 10x ," : ")')

    end subroutine

    subroutine broadcast_inputs

        use fileunits
        use inputlist

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        call MPI_BCAST(ext, 100, MPI_CHARACTER, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting physicslist     from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Igeometry, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Istellsym, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lfreebound, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(phiedge, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(curtor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(curpol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(gamma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Nfp, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Nvol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Mpol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ntor, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lrad, MNvol + 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(tflux, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pflux, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(helicity, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pscale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pressure, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ladiabatic, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(adiabatic, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(mu, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ivolume, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Isurf, MNvol + 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lconstraint, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pl, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(ql, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pr, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(qr, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iota, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(lp, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(lq, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(rp, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(rq, MNvol, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(oita, MNvol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(mupftol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(mupfits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lreflect, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(rpol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(rtor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting numericlist     from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Linitialize, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(LautoinitBn, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lzerovac, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ndiscrete, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Nquad, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iMpol, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iNtor, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lsparse, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lsvdiota, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(imethod, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iorder, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iprecon, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(iotatol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lextrap, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Mregular, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lrzaxis, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ntoraxis, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting globallist      from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Lfindzero, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(escale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(opsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(pcondense, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(epsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(wpoloidal, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(upsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(forcetol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(c05xmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(c05xtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(c05factor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(LreadGF, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(mfreeits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(gBntol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(gBnbld, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(vcasingeps, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(vcasingtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(vcasingits, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(vcasingper, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting locallist       from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(LBeltrami, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Linitgues, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(maxrndgues, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lmatsolver, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(NiterGMRES, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(epsGMRES, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(LGMRESprec, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(epsILU, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting diagnosticslist from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(odetol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(nPpts, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ppts, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(nPtrj, MNvol + 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)

        call MPI_BCAST(LHevalues, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(LHevectors, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ltransform, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(LHmatrix, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lperturbed, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(dpp, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(dqq, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lerrortype, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ngrid, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(dRZ, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Lcheck, 1, MPI_INTEGER, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Ltiming, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

        if (Wreadin) then; cput = MPI_WTIME(); write (ounit, '("readin : ",f10.2," : broadcasting screenlist      from ext.sp ;")') cput - cpus
        end if

        call MPI_BCAST(Wreadin, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Wwrtend, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)
        call MPI_BCAST(Wmacros, 1, MPI_LOGICAL, 0, MPI_COMM_SPEC, ierr)

    end subroutine

    subroutine wrtend

        use constants, only:

        use numerical, only: machprec

        use fileunits, only: ounit, iunit

        use cputiming, only: Twrtend

        use inputlist

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        integer :: vvol
        integer :: imn
        integer :: ii
        integer :: jj
        integer :: kk
        integer :: jk
        integer :: Lcurvature
        integer :: mm
        integer :: nn

        real(8) :: lss
        real(8) :: teta
        real(8) :: zeta
        real(8) :: st(1:Node)
        real(8) :: Bst(1:Node)
        real(8) :: BR
        real(8) :: BZ
        real(8) :: BP

        if (myid /= 0) return

        open (iunit, file=trim(ext)//".sp.end", status="unknown")

        write (iunit, '("&physicslist")')
        write (iunit, '(" Igeometry   = ",i9        )') Igeometry
        write (iunit, '(" Istellsym   = ",i9        )') Istellsym
        write (iunit, '(" Lfreebound  = ",i9        )') Lfreebound
        write (iunit, '(" phiedge     = ",es23.15   )') phiedge
        write (iunit, '(" curtor      = ",es23.15   )') curtor
        write (iunit, '(" curpol      = ",es23.15   )') curpol
        write (iunit, '(" gamma       = ",es23.15   )') gamma
        write (iunit, '(" Nfp         = ",i9        )') Nfp
        write (iunit, '(" Nvol        = ",i9        )') Nvol
        write (iunit, '(" Mpol        = ",i9        )') Mpol
        write (iunit, '(" Ntor        = ",i9        )') Ntor
        write (iunit, '(" Lrad        = ",257i23    )') Lrad(1:Mvol)
        write (iunit, '(" tflux       = ",257es23.15)') tflux(1:Mvol)
        write (iunit, '(" pflux       = ",257es23.15)') pflux(1:Mvol)
        write (iunit, '(" helicity    = ",256es23.15)') helicity(1:Mvol)
        write (iunit, '(" pscale      = ",es23.15   )') pscale
        write (iunit, '(" Ladiabatic  = ",i9        )') Ladiabatic
        write (iunit, '(" pressure    = ",257es23.15)') pressure(1:Mvol)
        write (iunit, '(" adiabatic   = ",257es23.15)') adiabatic(1:Mvol)
        write (iunit, '(" mu          = ",257es23.15)') mu(1:Mvol)
        write (iunit, '(" Ivolume     = ",257es23.15)') Ivolume(1:Mvol)
        write (iunit, '(" Isurf       = ",257es23.15)') Isurf(1:Mvol - 1), 0.0
        write (iunit, '(" Lconstraint = ",i9        )') Lconstraint
        write (iunit, '(" pl          = ",257i23    )') pl(0:Mvol)
        write (iunit, '(" ql          = ",257i23    )') ql(0:Mvol)
        write (iunit, '(" pr          = ",257i23    )') pr(0:Mvol)
        write (iunit, '(" qr          = ",257i23    )') qr(0:Mvol)
        write (iunit, '(" iota        = ",257es23.15)') iota(0:Mvol)
        write (iunit, '(" lp          = ",257i23    )') lp(0:Mvol)
        write (iunit, '(" lq          = ",257i23    )') lq(0:Mvol)
        write (iunit, '(" rp          = ",257i23    )') rp(0:Mvol)
        write (iunit, '(" rq          = ",257i23    )') rq(0:Mvol)
        write (iunit, '(" oita        = ",257es23.15)') oita(0:Mvol)
        write (iunit, '(" mupftol     = ",es23.15   )') mupftol
        write (iunit, '(" mupfits     = ",i9        )') mupfits
        write (iunit, '(" Lreflect    = ",i9        )') Lreflect
        write (iunit, '(" rpol        = ",es23.15   )') rpol
        write (iunit, '(" rtor        = ",es23.15   )') rtor

        write (iunit, '(" Rac         = ",99es23.15)') iRbc(1:Ntor + 1, 0)
        write (iunit, '(" Zas         = ",99es23.15)') iZbs(1:Ntor + 1, 0)
        write (iunit, '(" Ras         = ",99es23.15)') iRbs(1:Ntor + 1, 0)
        write (iunit, '(" Zac         = ",99es23.15)') iZbc(1:Ntor + 1, 0)

        do mm = 0, Mpol
            do nn = -Ntor, Ntor

                if (mm == 0 .and. nn < 0) cycle

                select case (mm)
                case (0:9)
                    if (nn < -9 .and. nn > -99) write (iunit, 1000) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1001) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1002) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1001) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                case (10:99)
                    if (nn < -9 .and. nn > -99) write (iunit, 1003) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1004) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1005) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1004) nn, mm, Rbc(nn, mm), nn, mm, Zbs(nn, mm), nn, mm, Rbs(nn, mm), nn, mm, Zbc(nn, mm)
                end select

            end do
        end do

        do mm = 0, Mpol
            do nn = -Ntor, Ntor

                if (mm == 0 .and. nn < 0) cycle

                select case (mm)
                case (0:9)
                    if (nn < -9 .and. nn > -99) write (iunit, 1010) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1011) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1012) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1011) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                case (10:99)
                    if (nn < -9 .and. nn > -99) write (iunit, 1013) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1014) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1015) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1014) nn, mm, Rwc(nn, mm), nn, mm, Zws(nn, mm), nn, mm, Rws(nn, mm), nn, mm, Zwc(nn, mm)
                end select

            end do
        end do

        do mm = 0, Mpol
            do nn = -Ntor, Ntor

                if (mm == 0 .and. nn < 0) cycle

                select case (mm)
                case (0:9)
                    if (nn < -9 .and. nn > -99) write (iunit, 1020) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1021) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1022) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1021) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                case (10:99)
                    if (nn < -9 .and. nn > -99) write (iunit, 1023) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn < 0 .and. nn >= -9) write (iunit, 1024) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn >= 0 .and. nn <= 9) write (iunit, 1025) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                    if (nn > 9 .and. nn <= 99) write (iunit, 1024) nn, mm, Vns(nn, mm), nn, mm, Bns(nn, mm), nn, mm, Vnc(nn, mm), nn, mm, Bnc(nn, mm)
                end select

            end do
        end do

1000    format("Rbc(", i3, ",", i1, ")", 2x, "=", es23.15, " Zbs(", i3, ",", i1, ")", 2x, "=", es23.15, " Rbs(", i3, ",", i1, ")", 2x, "=", es23.15, " Zbc(", i3, ",", i1, ")", 2x, "=", es23.15)
1001    format("Rbc(", i2, ",", i1, ")", 3x, "=", es23.15, " Zbs(", i2, ",", i1, ")", 3x, "=", es23.15, " Rbs(", i2, ",", i1, ")", 3x, "=", es23.15, " Zbc(", i2, ",", i1, ")", 3x, "=", es23.15)
1002    format("Rbc(", i1, ",", i1, ")", 4x, "=", es23.15, " Zbs(", i1, ",", i1, ")", 4x, "=", es23.15, " Rbs(", i1, ",", i1, ")", 4x, "=", es23.15, " Zbc(", i1, ",", i1, ")", 4x, "=", es23.15)
1003    format("Rbc(", i3, ",", i2, ")", 1x, "=", es23.15, " Zbs(", i3, ",", i2, ")", 1x, "=", es23.15, " Rbs(", i3, ",", i2, ")", 1x, "=", es23.15, " Zbc(", i3, ",", i2, ")", 1x, "=", es23.15)
1004    format("Rbc(", i2, ",", i2, ")", 2x, "=", es23.15, " Zbs(", i2, ",", i2, ")", 2x, "=", es23.15, " Rbs(", i2, ",", i2, ")", 2x, "=", es23.15, " Zbc(", i2, ",", i2, ")", 2x, "=", es23.15)
1005    format("Rbc(", i1, ",", i2, ")", 3x, "=", es23.15, " Zbs(", i1, ",", i2, ")", 3x, "=", es23.15, " Rbs(", i1, ",", i2, ")", 3x, "=", es23.15, " Zbc(", i1, ",", i2, ")", 3x, "=", es23.15)

1010    format("Rwc(", i3, ",", i1, ")", 2x, "=", es23.15, " Zws(", i3, ",", i1, ")", 2x, "=", es23.15, " Rws(", i3, ",", i1, ")", 2x, "=", es23.15, " Zwc(", i3, ",", i1, ")", 2x, "=", es23.15)
1011    format("Rwc(", i2, ",", i1, ")", 3x, "=", es23.15, " Zws(", i2, ",", i1, ")", 3x, "=", es23.15, " Rws(", i2, ",", i1, ")", 3x, "=", es23.15, " Zwc(", i2, ",", i1, ")", 3x, "=", es23.15)
1012    format("Rwc(", i1, ",", i1, ")", 4x, "=", es23.15, " Zws(", i1, ",", i1, ")", 4x, "=", es23.15, " Rws(", i1, ",", i1, ")", 4x, "=", es23.15, " Zwc(", i1, ",", i1, ")", 4x, "=", es23.15)
1013    format("Rwc(", i3, ",", i2, ")", 1x, "=", es23.15, " Zws(", i3, ",", i2, ")", 1x, "=", es23.15, " Rws(", i3, ",", i2, ")", 1x, "=", es23.15, " Zwc(", i3, ",", i2, ")", 1x, "=", es23.15)
1014    format("Rwc(", i2, ",", i2, ")", 2x, "=", es23.15, " Zws(", i2, ",", i2, ")", 2x, "=", es23.15, " Rws(", i2, ",", i2, ")", 2x, "=", es23.15, " Zwc(", i2, ",", i2, ")", 2x, "=", es23.15)
1015    format("Rwc(", i1, ",", i2, ")", 3x, "=", es23.15, " Zws(", i1, ",", i2, ")", 3x, "=", es23.15, " Rws(", i1, ",", i2, ")", 3x, "=", es23.15, " Zwc(", i1, ",", i2, ")", 3x, "=", es23.15)

1020    format("Vns(", i3, ",", i1, ")", 2x, "=", es23.15, " Bns(", i3, ",", i1, ")", 2x, "=", es23.15, " Vnc(", i3, ",", i1, ")", 2x, "=", es23.15, " Bnc(", i3, ",", i1, ")", 2x, "=", es23.15)
1021    format("Vns(", i2, ",", i1, ")", 3x, "=", es23.15, " Bns(", i2, ",", i1, ")", 3x, "=", es23.15, " Vnc(", i2, ",", i1, ")", 3x, "=", es23.15, " Bnc(", i2, ",", i1, ")", 3x, "=", es23.15)
1022    format("Vns(", i1, ",", i1, ")", 4x, "=", es23.15, " Bns(", i1, ",", i1, ")", 4x, "=", es23.15, " Vnc(", i1, ",", i1, ")", 4x, "=", es23.15, " Bnc(", i1, ",", i1, ")", 4x, "=", es23.15)
1023    format("Vns(", i3, ",", i2, ")", 1x, "=", es23.15, " Bns(", i3, ",", i2, ")", 1x, "=", es23.15, " Vnc(", i3, ",", i2, ")", 1x, "=", es23.15, " Bnc(", i3, ",", i2, ")", 1x, "=", es23.15)
1024    format("Vns(", i2, ",", i2, ")", 2x, "=", es23.15, " Bns(", i2, ",", i2, ")", 2x, "=", es23.15, " Vnc(", i2, ",", i2, ")", 2x, "=", es23.15, " Bnc(", i2, ",", i2, ")", 2x, "=", es23.15)
1025    format("Vns(", i1, ",", i2, ")", 3x, "=", es23.15, " Bns(", i1, ",", i2, ")", 3x, "=", es23.15, " Vnc(", i1, ",", i2, ")", 3x, "=", es23.15, " Bnc(", i1, ",", i2, ")", 3x, "=", es23.15)

        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput - cpus, myid
        end if

        write (iunit, '("&numericlist")')
        write (iunit, '(" Linitialize = ",i9            )') Linitialize
        write (iunit, '(" LautoinitBn = ",i9            )') LautoinitBn
        write (iunit, '(" Lzerovac    = ",i9            )') Lzerovac
        write (iunit, '(" Ndiscrete   = ",i9            )') Ndiscrete
        write (iunit, '(" Nquad       = ",i9            )') Nquad
        write (iunit, '(" iMpol       = ",i9            )') iMpol
        write (iunit, '(" iNtor       = ",i9            )') iNtor
        write (iunit, '(" Lsparse     = ",i9            )') Lsparse
        write (iunit, '(" Lsvdiota    = ",i9            )') Lsvdiota
        write (iunit, '(" imethod     = ",i9            )') imethod
        write (iunit, '(" iorder      = ",i9            )') iorder
        write (iunit, '(" iprecon     = ",i9            )') iprecon
        write (iunit, '(" iotatol     = ",es23.15       )') iotatol
        write (iunit, '(" Lextrap     = ",i9            )') Lextrap
        write (iunit, '(" Mregular    = ",i9            )') Mregular
        write (iunit, '(" Lrzaxis     = ",i9            )') Lrzaxis
        write (iunit, '(" Ntoraxis    = ",i9            )') Ntoraxis
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput - cpus, myid
        end if

        write (iunit, '("&locallist")')
        write (iunit, '(" LBeltrami   = ",i9            )') LBeltrami
        write (iunit, '(" Linitgues   = ",i9            )') Linitgues
        write (iunit, '(" Lmatsolver  = ",i9            )') Lmatsolver
        write (iunit, '(" NiterGMRES  = ",i9            )') NiterGMRES
        write (iunit, '(" LGMRESprec  = ",i9            )') LGMRESprec
        write (iunit, '(" epsGMRES    = ",es23.15       )') epsGMRES
        write (iunit, '(" epsILU      = ",es23.15       )') epsILU

        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing globallist ;")') cput - cpus, myid
        end if

        write (iunit, '("&globallist")')
        write (iunit, '(" Lfindzero   = ",i9            )') Lfindzero
        write (iunit, '(" escale      = ",es23.15       )') escale
        write (iunit, '(" opsilon     = ",es23.15       )') opsilon
        write (iunit, '(" pcondense   = ",es23.15       )') pcondense
        write (iunit, '(" epsilon     = ",es23.15       )') epsilon
        write (iunit, '(" wpoloidal   = ",es23.15       )') wpoloidal
        write (iunit, '(" upsilon     = ",es23.15       )') upsilon
        write (iunit, '(" forcetol    = ",es23.15       )') forcetol
        write (iunit, '(" c05xmax     = ",es23.15       )') c05xmax
        write (iunit, '(" c05xtol     = ",es23.15       )') c05xtol
        write (iunit, '(" c05factor   = ",es23.15       )') c05factor
        write (iunit, '(" LreadGF     = ",L9            )') LreadGF
        write (iunit, '(" mfreeits    = ",i9            )') mfreeits
        write (iunit, '(" gBntol      = ",es23.15       )') gBntol
        write (iunit, '(" gBnbld      = ",es23.15       )') gBnbld
        write (iunit, '(" vcasingeps  = ",es23.15       )') vcasingeps
        write (iunit, '(" vcasingtol  = ",es23.15       )') vcasingtol
        write (iunit, '(" vcasingits  = ",i9            )') vcasingits
        write (iunit, '(" vcasingper  = ",i9            )') vcasingper
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput - cpus, myid
        end if

        write (iunit, '("&diagnosticslist")')
        write (iunit, '(" odetol      = ",es23.15       )') odetol
        write (iunit, '(" nPpts       = ",i9            )') nPpts
        write (iunit, '(" Ppts        = ",es23.15       )') Ppts
        write (iunit, '(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
        write (iunit, '(" LHevalues   = ",L9            )') LHevalues
        write (iunit, '(" LHevectors  = ",L9            )') LHevectors
        write (iunit, '(" LHmatrix    = ",L9            )') LHmatrix
        write (iunit, '(" Lperturbed  = ",i9            )') Lperturbed
        write (iunit, '(" dpp         = ",i9            )') dpp
        write (iunit, '(" dqq         = ",i9            )') dqq
        write (iunit, '(" dRZ         = ",es23.15       )') dRZ
        write (iunit, '(" Lcheck      = ",i9            )') Lcheck
        write (iunit, '(" Ltiming     = ",L9            )') Ltiming
        write (iunit, '("/")')

        if (Wwrtend) then; cput = MPI_WTIME(); write (ounit, '("wrtend : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput - cpus, myid
        end if

        write (iunit, '("&screenlist")')
        if (Wreadin) write (iunit, '(" Wreadin = ",L1                )') Wreadin
        if (Wwrtend) write (iunit, '(" Wwrtend = ",L1                )') Wwrtend
        if (Wmacros) write (iunit, '(" Wmacros = ",L1                )') Wmacros
        write (iunit, '("/")')

        do imn = 1, mn; write (iunit, '(2i6,1024es23.15)') im(imn), in(imn)/Nfp, (iRbc(imn, vvol), iZbs(imn, vvol), iRbs(imn, vvol), iZbc(imn, vvol), vvol=1, Nvol)
        end do

        close (iunit)

    end subroutine wrtend

    subroutine IsMyVolume(vvol)

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        integer, intent(in) :: vvol

        IsMyVolumeValue = -1
        if (myid /= modulo(vvol - 1, ncpu)) then
            IsMyVolumeValue = 0
        else
            IsMyVolumeValue = 1
        end if

    end subroutine IsMyVolume

    subroutine WhichCpuID(vvol, cpu_id)

        use mpi
        implicit none
        integer :: ierr, astat, ios, nthreads, ithread
        real(8) :: cput, cpui, cpuo = 0

        integer :: vvol, cpu_id

        cpu_id = modulo(vvol - 1, ncpu)

    end subroutine WhichCpuID

end module allglobal

module fftw_interface

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    type(c_ptr) :: planf
    type(c_ptr) :: planb
    complex(c_double_complex), allocatable :: cplxin(:, :, :)
    complex(c_double_complex), allocatable :: cplxout(:, :, :)

end module fftw_interface

