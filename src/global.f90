
module constants

  implicit none

  REAL, parameter :: zero       =    0.0 !< 0
  REAL, parameter :: one        =    1.0 !< 1
  REAL, parameter :: two        =    2.0 !< 2
  REAL, parameter :: three      =    3.0 !< 3
  REAL, parameter :: four       =    4.0 !< 4
  REAL, parameter :: five       =    5.0 !< 5
  REAL, parameter :: six        =    6.0 !< 6
  REAL, parameter :: seven      =    7.0 !< 7
  REAL, parameter :: eight      =    8.0 !< 8
  REAL, parameter :: nine       =    9.0 !< 9
  REAL, parameter :: ten        =   10.0 !< 10

  REAL, parameter :: eleven     =   11.0 !< 11
  REAL, parameter :: twelve     =   12.0 !< 12

  REAL, parameter :: hundred    =  100.0 !< 100
  REAL, parameter :: thousand   = 1000.0 !< 1000

  REAL, parameter :: half       =   one / two   !< 1/2
  REAL, parameter :: third      =   one / three !< 1/3
  REAL, parameter :: quart      =   one / four  !< 1/4
  REAL, parameter :: fifth      =   one / five  !< 1/5
  REAL, parameter :: sixth      =   one / six   !< 1/6

  REAL, parameter :: pi2        =   6.28318530717958623 !< \f$2\pi\f$
  REAL, parameter :: pi         =   pi2 / two           !< \f$\pi\f$
  REAL, parameter :: mu0        =   2.0E-07 * pi2       !< \f$4\pi\cdot10^{-7}\f$
  REAL, parameter :: goldenmean =   1.618033988749895   !< golden mean = \f$( 1 + \sqrt 5 ) / 2\f$ ;

  REAL, parameter :: version    =   3.23  !< version of SPEC

end module constants



module numerical

  implicit none

  REAL, parameter :: machprec = 1.11e-16           !< machine precision: 0.5*epsilon(one) for 64 bit double precision
  REAL, parameter :: vsmall = 100*machprec         !< very small number
  REAL, parameter :: small = 10000*machprec        !< small number
  REAL, parameter :: sqrtmachprec = sqrt(machprec) !< square root of machine precision
  REAL, parameter :: logtolerance = 1.0e-32        !< this is used to avoid taking alog10(zero); see e.g. dforce;

end module numerical



module fileunits

  implicit none

  INTEGER :: iunit = 10 !< input; used in global/readin:ext.sp, global/wrtend:ext.sp.end
  INTEGER :: ounit =  6 !< screen output;
  INTEGER :: gunit = 13 !< wall geometry; used in wa00aa

  INTEGER :: aunit = 11 !< vector potential; used in ra00aa:.ext.AtAzmn;
  INTEGER :: dunit = 12 !< derivative matrix; used in newton:.ext.GF;
  INTEGER :: hunit = 14 !< eigenvalues of Hessian; under re-construction;
  INTEGER :: munit = 14 !< matrix elements of Hessian;
  INTEGER :: lunit = 20 !< local unit; used in lunit+myid: pp00aa:.ext.poincare,.ext.transform;
  INTEGER :: vunit = 15 !< for examination of adaptive quadrature; used in casing:.ext.vcint;

  contains
    subroutine mute(action)
      implicit none

      INTEGER,intent(in) :: action
      INTEGER, parameter :: iopen = 1, iclose = 0, null = 37
      INTEGER            :: ios
      character(len=*), parameter :: nullfile="/dev/null"

      if (action == iopen) then
        ounit = null
        open(ounit, file=nullfile, status="unknown", action="write", iostat=ios) ! create a scratch file?
        if (ios.ne.0) print *, "something wrong with open a tmp file in focuspy.mute. IOSTAT=", ios
      else
        close(ounit)
        ounit = 6 ! recover to screen output
      endif
      return
    end subroutine mute

end module fileunits

module cputiming

  REAL    :: Tmanual = 0.0, manualT = 0.0
  REAL    :: Trzaxis = 0.0, rzaxisT = 0.0
  REAL    :: Tpackxi = 0.0, packxiT = 0.0
  REAL    :: Tvolume = 0.0, volumeT = 0.0
  REAL    :: Tcoords = 0.0, coordsT = 0.0
  REAL    :: Tbasefn = 0.0, basefnT = 0.0
  REAL    :: Tmemory = 0.0, memoryT = 0.0
  REAL    :: Tmetrix = 0.0, metrixT = 0.0
  REAL    :: Tma00aa = 0.0, ma00aaT = 0.0
  REAL    :: Tmatrix = 0.0, matrixT = 0.0
  REAL    :: Tmp00ac = 0.0, mp00acT = 0.0
  REAL    :: Tma02aa = 0.0, ma02aaT = 0.0
  REAL    :: Tpackab = 0.0, packabT = 0.0
  REAL    :: Ttr00ab = 0.0, tr00abT = 0.0
  REAL    :: Tcurent = 0.0, curentT = 0.0
  REAL    :: Tdf00ab = 0.0, df00abT = 0.0
  REAL    :: Tlforce = 0.0, lforceT = 0.0
  REAL    :: Tforce_real = 0.0, force_realT = 0.0
  REAL    :: Tforce_real_helper = 0.0, force_real_helperT = 0.0
  REAL    :: Tforce_real_dforce = 0.0, force_real_dforceT = 0.0
  REAL    :: Tforce_real_dfp200 = 0.0, force_real_dfp200T = 0.0
  REAL    :: Tforce_real_lforce = 0.0, force_real_lforceT = 0.0
  REAL    :: Tbcast_freal_jac = 0.0, bcast_freal_jacT = 0.0
  REAL    :: Tget_jac = 0.0, get_jacT = 0.0
  REAL    :: Tget_angle_constraint = 0.0, get_angle_constraintT = 0.0
  REAL    :: Tangconst_dfp200 = 0.0, angconst_dfp200T = 0.0
  REAL    :: Tintghs = 0.0, intghsT = 0.0
  REAL    :: Tmtrxhs = 0.0, mtrxhsT = 0.0
  REAL    :: Tlbpol = 0.0, lbpolT = 0.0
  REAL    :: Tbrcast = 0.0, brcastT = 0.0
  REAL    :: Tdfp100 = 0.0, dfp100T = 0.0
  REAL    :: Tdfp200 = 0.0, dfp200T = 0.0
  REAL    :: Tdforce = 0.0, dforceT = 0.0
  REAL    :: Tjo00aa = 0.0, jo00aaT = 0.0
  REAL    :: Tpp00aa = 0.0, pp00aaT = 0.0
  REAL    :: Tpp00ab = 0.0, pp00abT = 0.0
  REAL    :: Tbfield = 0.0, bfieldT = 0.0
  REAL    :: Tstzxyz = 0.0, stzxyzT = 0.0
  REAL    :: Thesian = 0.0, hesianT = 0.0
  REAL    :: Tra00aa = 0.0, ra00aaT = 0.0
  REAL    :: Tnumrec = 0.0, numrecT = 0.0
  REAL    :: Tdcuhre = 0.0, dcuhreT = 0.0
  REAL    :: Tminpack = 0.0, minpackT = 0.0
  REAL    :: Tiqpack = 0.0, iqpackT = 0.0
  REAL    :: Trksuite = 0.0, rksuiteT = 0.0
  REAL    :: Ti1mach = 0.0, i1machT = 0.0
  REAL    :: Td1mach = 0.0, d1machT = 0.0
  REAL    :: Tilut = 0.0, ilutT = 0.0
  REAL    :: Titers = 0.0, itersT = 0.0
  REAL    :: Tsphdf5 = 0.0, sphdf5T = 0.0
  REAL    :: Tpreset = 0.0, presetT = 0.0
  REAL    :: Tglobal = 0.0, globalT = 0.0
  REAL    :: Txspech = 0.0, xspechT = 0.0
  REAL    :: Tinputlist = 0.0, inputlistT = 0.0

  REAL :: Treadin = 0.0
  REAL :: Twrtend = 0.0

end module cputiming



module typedefns

  type subgrid
    REAL,    allocatable :: s(:) !< coefficients
    INTEGER, allocatable :: i(:) !< indices
  end type subgrid

  type MatrixLU
    REAL, allocatable :: mat(:,:)
    INTEGER, allocatable :: ipivot(:)
  end type MatrixLU

  type derivative
     LOGICAL :: L      !< what is this?
     INTEGER :: vol    !< Used in coords(); required for global constraint force gradient evaluation
     INTEGER :: innout !< what is this?
     INTEGER :: ii     !< what is this?
     INTEGER :: irz    !< what is this?
     INTEGER :: issym  !< what is this?
  end type derivative

end module typedefns



module allglobal

  use constants
  use typedefns

  implicit none


  INTEGER              :: myid !< MPI rank of current CPU
  INTEGER              :: ncpu !< number of MPI tasks
  INTEGER              :: IsMyVolumeValue !< flag to indicate if a CPU is operating on its assigned volume
  REAL                 :: cpus !< initial time
  INTEGER              :: MPI_COMM_SPEC !< SPEC MPI communicator


  LOGICAL              :: skip_write = .false. ! flag to disable any HDF5-related calls

  REAL                 :: pi2nfp           !       pi2/nfp     ; assigned in readin;
  REAL                 :: pi2pi2nfp
  REAL                 :: pi2pi2nfphalf
  REAL                 :: pi2pi2nfpquart



  CHARACTER(LEN=1000)  :: ext ! extension of input filename, i.e., "G3V01L1Fi.001" for an input file G3V01L1Fi.001.sp

  REAL                 :: ForceErr !< total force-imbalance
  REAL                 :: Energy   !< MHD energy
  REAL                 :: BnsErr   !< (in freeboundary) error in self-consistency of field on plasma boundary (Picard iteration)

  REAL   , allocatable :: IPDt(:), IPDtDpf(:,:)  !< Toroidal pressure-driven current

  INTEGER              :: Mvol

  LOGICAL              :: YESstellsym !< internal shorthand copies of Istellsym, which is an integer input;
  LOGICAL              :: NOTstellsym !< internal shorthand copies of Istellsym, which is an integer input;

  REAL   , allocatable :: cheby(:,:) !< local workspace for evaluation of Chebychev polynomials
  REAL   , allocatable :: zernike(:,:,:) !< local workspace for evaluation of Zernike polynomials

  REAL   , allocatable :: TT(:,:,:)    !< derivatives of Chebyshev polynomials at the inner and outer interfaces;
  REAL   , allocatable :: RTT(:,:,:,:) !< derivatives of Zernike   polynomials at the inner and outer interfaces;

  REAL   , allocatable :: RTM(:,:) !< \f$r^m\f$ term of Zernike polynomials at the origin
  REAL   , allocatable :: ZernikeDof(:) !< Zernike degree of freedom for each \f$m\f$



  INTEGER              :: mne    !< enhanced resolution for metric elements
  INTEGER, allocatable :: ime(:) !< enhanced poloidal mode numbers for metric elements
  INTEGER, allocatable :: ine(:) !< enhanced toroidal mode numbers for metric elements

  INTEGER              :: mns    !< enhanced resolution for straight field line transformation
  INTEGER, allocatable :: ims(:) !< enhanced poloidal mode numbers for straight field line transformation
  INTEGER, allocatable :: ins(:) !< enhanced toroidal mode numbers for straight field line transformation

  INTEGER              :: lMpol !< what is this?
  INTEGER              :: lNtor !< what is this?
  INTEGER              :: sMpol !< what is this?
  INTEGER              :: sNtor !< what is this?



  REAL                 :: xoffset = 1.0 !< used to normalize NAG routines (which ones exacly where?)



  LOGICAL, allocatable :: ImagneticOK(:) !< used to indicate if Beltrami fields have been correctly constructed;

  LOGICAL              :: IconstraintOK !< Used to break iteration loops of slaves in the global constraint minimization.

  REAL   , allocatable :: beltramierror(:,:)  !< to store the integral of |curlB-mu*B| computed by jo00aa;



  INTEGER              :: mn    !< total number of Fourier harmonics for coordinates/fields; calculated from Mpol, Ntor in readin()
  INTEGER, allocatable :: im(:) !< poloidal mode numbers for Fourier representation
  INTEGER, allocatable :: in(:) !< toroidal mode numbers for Fourier representation

  REAL,    allocatable :: halfmm(:) !< I saw this already somewhere...
  REAL,    allocatable :: regumm(:) !< I saw this already somewhere...

  REAL                 :: Rscale    !< no idea
  REAL,    allocatable :: psifactor(:,:) !< no idea
  REAL,    allocatable :: inifactor(:,:) !< no idea

  REAL,    allocatable :: BBweight(:) !< weight on force-imbalance harmonics; used in dforce()

  REAL,    allocatable :: mmpp(:) !< spectral condensation factors



  REAL,    allocatable :: iRbc(:,:) !< cosine R harmonics of interface surface geometry;     stellarator symmetric
  REAL,    allocatable :: iZbs(:,:) !<   sine Z harmonics of interface surface geometry;     stellarator symmetric
  REAL,    allocatable :: iRbs(:,:) !<   sine R harmonics of interface surface geometry; non-stellarator symmetric
  REAL,    allocatable :: iZbc(:,:) !< cosine Z harmonics of interface surface geometry; non-stellarator symmetric

  REAL,    allocatable :: dRbc(:,:) !< cosine R harmonics of interface surface geometry;     stellarator symmetric; linear deformation
  REAL,    allocatable :: dZbs(:,:) !<   sine Z harmonics of interface surface geometry;     stellarator symmetric; linear deformation
  REAL,    allocatable :: dRbs(:,:) !<   sine R harmonics of interface surface geometry; non-stellarator symmetric; linear deformation
  REAL,    allocatable :: dZbc(:,:) !< cosine Z harmonics of interface surface geometry; non-stellarator symmetric; linear deformation

  REAL,    allocatable :: iRij(:,:) !< interface surface geometry; real space
  REAL,    allocatable :: iZij(:,:) !< interface surface geometry; real space
  REAL,    allocatable :: dRij(:,:) !< interface surface geometry; real space
  REAL,    allocatable :: dZij(:,:) !< interface surface geometry; real space
  REAL,    allocatable :: tRij(:,:) !< interface surface geometry; real space
  REAL,    allocatable :: tZij(:,:) !< interface surface geometry; real space

  REAL,    allocatable :: iVns(:)   !<   sine harmonics of vacuum normal magnetic field on interfaces;     stellarator symmetric
  REAL,    allocatable :: iBns(:)   !<   sine harmonics of plasma normal magnetic field on interfaces;     stellarator symmetric
  REAL,    allocatable :: iVnc(:)   !< cosine harmonics of vacuum normal magnetic field on interfaces; non-stellarator symmetric
  REAL,    allocatable :: iBnc(:)   !< cosine harmonics of plasma normal magnetic field on interfaces; non-stellarator symmetric

  REAL,    allocatable :: lRbc(:)   !< local workspace
  REAL,    allocatable :: lZbs(:)   !< local workspace
  REAL,    allocatable :: lRbs(:)   !< local workspace
  REAL,    allocatable :: lZbc(:)   !< local workspace

  INTEGER              :: num_modes
  INTEGER, allocatable :: mmRZRZ(:), nnRZRZ(:)
  REAL,    allocatable :: allRZRZ(:,:,:)



  INTEGER              :: Nt  !< discrete resolution along \f$\theta\f$ of grid in real space
  INTEGER              :: Nz  !< discrete resolution along \f$\zeta\f$  of grid in real space
  INTEGER              :: Ntz !< discrete resolution; Ntz=Nt*Nz shorthand
  INTEGER              :: hNt !< discrete resolution; Ntz=Nt*Nz shorthand
  INTEGER              :: hNz !< discrete resolution; Ntz=Nt*Nz shorthand
  REAL                 :: soNtz !< one / sqrt (one*Ntz); shorthand

  REAL   , allocatable :: Rij(:,:,:) !< real-space grid; R
  REAL   , allocatable :: Zij(:,:,:) !< real-space grid; Z
  REAL   , allocatable :: Xij(:,:,:) !< what is this?
  REAL   , allocatable :: Yij(:,:,:) !< what is this?
  REAL   , allocatable :: sg(:,:)    !< real-space grid; jacobian and its derivatives
  REAL   , allocatable :: guvij(:,:,:,:) !< real-space grid; metric elements
  REAL   , allocatable :: gvuij(:,:,:)   !< real-space grid; metric elements (?); 10 Dec 15;
  REAL   , allocatable :: guvijsave(:,:,:,:) !< what is this?

  INTEGER, allocatable :: ki(:,:)     !< identification of Fourier modes
  INTEGER, allocatable :: kijs(:,:,:) !< identification of Fourier modes
  INTEGER, allocatable :: kija(:,:,:) !< identification of Fourier modes

  INTEGER, allocatable :: iotakkii(:)   !< identification of Fourier modes
  INTEGER, allocatable :: iotaksub(:,:) !< identification of Fourier modes
  INTEGER, allocatable :: iotakadd(:,:) !< identification of Fourier modes
  INTEGER, allocatable :: iotaksgn(:,:) !< identification of Fourier modes

  REAL   , allocatable :: efmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: ofmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: cfmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: sfmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: evmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: odmn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: comn(:) !< Fourier harmonics; dummy workspace
  REAL   , allocatable :: simn(:) !< Fourier harmonics; dummy workspace

  REAL   , allocatable :: ijreal(:) !< what is this ?
  REAL   , allocatable :: ijimag(:) !< what is this ?
  REAL   , allocatable :: jireal(:) !< what is this ?
  REAL   , allocatable :: jiimag(:) !< what is this ?

  REAL   , allocatable :: jkreal(:) !< what is this ?
  REAL   , allocatable :: jkimag(:) !< what is this ?
  REAL   , allocatable :: kjreal(:) !< what is this ?
  REAL   , allocatable :: kjimag(:) !< what is this ?

  REAL   , allocatable :: Bsupumn(:,:,:) !< tangential field on interfaces; \f$\theta\f$-component; required for virtual casing construction of field; 11 Oct 12
  REAL   , allocatable :: Bsupvmn(:,:,:) !< tangential field on interfaces; \f$\zeta\f$ -component; required for virtual casing construction of field; 11 Oct 12



  REAL   , allocatable :: goomne(:,:) !< described in preset()
  REAL   , allocatable :: goomno(:,:) !< described in preset()
  REAL   , allocatable :: gssmne(:,:) !< described in preset()
  REAL   , allocatable :: gssmno(:,:) !< described in preset()
  REAL   , allocatable :: gstmne(:,:) !< described in preset()
  REAL   , allocatable :: gstmno(:,:) !< described in preset()
  REAL   , allocatable :: gszmne(:,:) !< described in preset()
  REAL   , allocatable :: gszmno(:,:) !< described in preset()
  REAL   , allocatable :: gttmne(:,:) !< described in preset()
  REAL   , allocatable :: gttmno(:,:) !< described in preset()
  REAL   , allocatable :: gtzmne(:,:) !< described in preset()
  REAL   , allocatable :: gtzmno(:,:) !< described in preset()
  REAL   , allocatable :: gzzmne(:,:) !< described in preset()
  REAL   , allocatable :: gzzmno(:,:) !< described in preset()



  REAL,    allocatable :: DToocc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DToocs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DToosc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DTooss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TTsscc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TTsscs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TTsssc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TTssss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDstcc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDstcs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDstsc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDstss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDszcc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDszcs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDszsc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: TDszss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDttcc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDttcs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDttsc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDttss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDtzcc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDtzcs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDtzsc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDtzss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDzzcc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDzzcs(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDzzsc(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()
  REAL,    allocatable :: DDzzss(:,:,:,:) !< volume-integrated Chebychev-metrics; see matrix()



  REAL,    allocatable :: Tsc(:,:) !< what is this?
  REAL,    allocatable :: Tss(:,:) !< what is this?
  REAL,    allocatable :: Dtc(:,:) !< what is this?
  REAL,    allocatable :: Dts(:,:) !< what is this?
  REAL,    allocatable :: Dzc(:,:) !< what is this?
  REAL,    allocatable :: Dzs(:,:) !< what is this?
  REAL,    allocatable :: Ttc(:,:) !< what is this?
  REAL,    allocatable :: Tzc(:,:) !< what is this?
  REAL,    allocatable :: Tts(:,:) !< what is this?
  REAL,    allocatable :: Tzs(:,:) !< what is this?



  REAL,    allocatable :: dtflux(:) !< \f$\delta \psi_{toroidal}\f$ in each annulus
  REAL,    allocatable :: dpflux(:) !< \f$\delta \psi_{poloidal}\f$ in each annulus



  REAL,    allocatable :: sweight(:) !< minimum poloidal length constraint weight



  INTEGER, allocatable :: NAdof(:) !< degrees of freedom in Beltrami fields in each annulus
  INTEGER, allocatable :: Nfielddof(:) !< degrees of freedom in Beltrami fields in each annulus, field only, no Lagrange multipliers

  type(subgrid), allocatable :: Ate(:,:,:) !< magnetic vector potential cosine Fourier harmonics;     stellarator-symmetric
  type(subgrid), allocatable :: Aze(:,:,:) !< magnetic vector potential cosine Fourier harmonics;     stellarator-symmetric
  type(subgrid), allocatable :: Ato(:,:,:) !< magnetic vector potential   sine Fourier harmonics; non-stellarator-symmetric
  type(subgrid), allocatable :: Azo(:,:,:) !< magnetic vector potential   sine Fourier harmonics; non-stellarator-symmetric

  INTEGER      , allocatable :: Lma(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmb(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmc(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmd(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lme(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmf(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmg(:,:) !< Lagrange multipliers (?)
  INTEGER      , allocatable :: Lmh(:,:) !< Lagrange multipliers (?)

  REAL         , allocatable :: Lmavalue(:,:) !< what is this?
  REAL         , allocatable :: Lmbvalue(:,:) !< what is this?
  REAL         , allocatable :: Lmcvalue(:,:) !< what is this?
  REAL         , allocatable :: Lmdvalue(:,:) !< what is this?
  REAL         , allocatable :: Lmevalue(:,:) !< what is this?
  REAL         , allocatable :: Lmfvalue(:,:) !< what is this?
  REAL         , allocatable :: Lmgvalue(:,:) !< what is this?
  REAL         , allocatable :: Lmhvalue(:,:) !< what is this?

  INTEGER      , allocatable :: Fso(:,:) !< what is this?
  INTEGER      , allocatable :: Fse(:,:) !< what is this?

  LOGICAL                    :: Lcoordinatesingularity !< set by \c LREGION macro; true if inside the innermost volume
  LOGICAL                    :: Lplasmaregion          !< set by \c LREGION macro; true if inside the plasma region
  LOGICAL                    :: Lvacuumregion          !< set by \c LREGION macro; true if inside the vacuum region
  LOGICAL                    :: Lsavedguvij            !< flag used in matrix free
  LOGICAL                    :: Localconstraint        !< what is this?


  REAL,   allocatable :: Remn_ext(:,:,:)
  REAL,   allocatable :: Romn_ext(:,:,:)
  REAL,   allocatable :: Zemn_ext(:,:,:)
  REAL,   allocatable :: Zomn_ext(:,:,:)

  LOGICAL             :: use_ext_mesh = .false.
  


   REAL,   allocatable :: dMA(:,:) !< energy and helicity matrices; quadratic forms
   REAL,   allocatable :: dMB(:,:) !< energy and helicity matrices; quadratic forms
   REAL,   allocatable :: dMD(:,:) !< energy and helicity matrices; quadratic forms

   REAL,   allocatable :: dMAS(:)     !< sparse version of dMA, data
   REAL,   allocatable :: dMDS(:)     !< sparse version of dMD, data
   INTEGER,allocatable :: idMAS(:)    !< sparse version of dMA and dMD, indices
   INTEGER,allocatable :: jdMAS(:)    !< sparse version of dMA and dMD, indices
   INTEGER,allocatable :: NdMASmax(:) !< number of elements for sparse matrices
   INTEGER,allocatable :: NdMAS(:)    !< number of elements for sparse matrices

   REAL,   allocatable :: dMG(:  ) !< what is this?

   REAL,   allocatable :: AdotX(:) !< the matrix-vector product
   REAL,   allocatable :: DdotX(:) !< the matrix-vector product

   REAL,   allocatable :: solution(:,:) !< this is allocated in dforce; used in mp00ac and ma02aa; and is passed to packab

   REAL,   allocatable :: GMRESlastsolution(:,:,:) !< used to store the last solution for restarting GMRES

   REAL,   allocatable :: MBpsi(:)      !< matrix vector products

   REAL,   allocatable :: BeltramiInverse(:,:) !< Beltrami inverse matrix



  REAL   , allocatable :: diotadxup(:,:,:) !< measured rotational transform on inner/outer interfaces for each volume;          d(transform)/dx; (see dforce)
  REAL   , allocatable :: dItGpdxtp(:,:,:) !< measured toroidal and poloidal current on inner/outer interfaces for each volume; d(Itor,Gpol)/dx; (see dforce)

  REAL   , allocatable :: glambda(:,:,:,:) !< save initial guesses for iterative calculation of rotational-transform

  INTEGER              :: lmns !< number of independent degrees of freedom in angle transformation;

  REAL,    allocatable :: dlambdaout(:,:,:)



  REAL,    allocatable ::  Bemn(:,:,:) !< force vector;     stellarator-symmetric (?)
  REAL,    allocatable ::  Iomn(:,:)   !< force vector;     stellarator-symmetric (?)
  REAL,    allocatable ::  Somn(:,:,:) !< force vector; non-stellarator-symmetric (?)
  REAL,    allocatable ::  Pomn(:,:,:) !< force vector; non-stellarator-symmetric (?)

  REAL,    allocatable ::  Bomn(:,:,:) !< force vector;     stellarator-symmetric (?)
  REAL,    allocatable ::  Iemn(:,:)   !< force vector;     stellarator-symmetric (?)
  REAL,    allocatable ::  Semn(:,:,:) !< force vector; non-stellarator-symmetric (?)
  REAL,    allocatable ::  Pemn(:,:,:) !< force vector; non-stellarator-symmetric (?)

  REAL,    allocatable ::  BBe(:) !< force vector (?);     stellarator-symmetric (?)
  REAL,    allocatable ::  IIo(:) !< force vector (?);     stellarator-symmetric (?)
  REAL,    allocatable ::  BBo(:) !< force vector (?); non-stellarator-symmetric (?)
  REAL,    allocatable ::  IIe(:) !< force vector (?); non-stellarator-symmetric (?)


  REAL,    allocatable ::  Btemn(:,:,:) !< covariant \f$\theta\f$ cosine component of the tangential field on interfaces;     stellarator-symmetric
  REAL,    allocatable ::  Bzemn(:,:,:) !< covariant \f$\zeta\f$  cosine component of the tangential field on interfaces;     stellarator-symmetric
  REAL,    allocatable ::  Btomn(:,:,:) !< covariant \f$\theta\f$   sine component of the tangential field on interfaces; non-stellarator-symmetric
  REAL,    allocatable ::  Bzomn(:,:,:) !< covariant \f$\zeta\f$    sine component of the tangential field on interfaces; non-stellarator-symmetric


  REAL,    allocatable ::  Bloweremn(:,:) !< covariant field for Hessian computation
  REAL,    allocatable ::  Bloweromn(:,:) !< covariant field for Hessian computation


  INTEGER              :: LGdof !<       geometrical degrees of freedom associated with each interface
  INTEGER              :: NGdof !< total geometrical degrees of freedom


  REAL,    allocatable :: dBBdRZ(:,:,:) !< derivative of magnetic field w.r.t. geometry (?)
  REAL,    allocatable :: dIIdRZ(:  ,:) !< derivative of spectral constraints w.r.t. geometry (?)

  REAL,    allocatable :: dFFdRZ(:,:,:,:,:) !< derivatives of B^2 at the interfaces wrt geometry
  REAL,    allocatable :: dBBdmp(:,:,:,:  ) !< derivatives of B^2 at the interfaces wrt mu and dpflux
  
  REAL,    allocatable :: freal_jac(:,:,:,:,:)
  REAL,    allocatable :: freal_dBBdmp(:,:,:,:  )
  INTEGER   :: dbbdmp_in_frealjac = 1

  REAL,    allocatable :: HdFFdRZ(:,:,:,:,:) !< derivatives of B^2 at the interfaces wrt geometry 2D Hessian; 

  REAL,    allocatable :: denergydrr(:,:,:,:,:) !< derivatives of energy at the interfaces wrt geometry 3D Hessian; 
  REAL,    allocatable :: denergydrz(:,:,:,:,:) !< derivatives of energy at the interfaces wrt geometry 3D Hessian; 
  REAL,    allocatable :: denergydzr(:,:,:,:,:) !< derivatives of energy at the interfaces wrt geometry 3D Hessian; 
  REAL,    allocatable :: denergydzz(:,:,:,:,:) !< derivatives of energy at the interfaces wrt geometry 3D Hessian; 



  REAL,    allocatable :: dmupfdx(:,:,:,:,:)  !< derivatives of mu and dpflux wrt geometry at constant interface transform



  REAL,    allocatable :: freal_jac_full(:,:)
  REAL,    allocatable :: dessian(:,:)      !< derivative of force gradient matrix (?)

  REAL,    allocatable :: force_final(:) !< Final force on the interfaces [inface*mode]


  REAL   , allocatable :: cosi(:,:) !< some precomputed cosines
  REAL   , allocatable :: sini(:,:) !< some precomputed sines
  REAL   , allocatable :: gteta(:)  !< something related to \f$\sqrt g\f$ and \f$\theta\f$ ?
  REAL   , allocatable :: gzeta(:)  !< something related to \f$\sqrt g\f$ and \f$\zeta\f$ ?

  REAL   , allocatable :: ajk(:)    !< definition of coordinate axis

  REAL   , allocatable :: dRadR(:,:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dRadZ(:,:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dZadR(:,:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dZadZ(:,:,:,:) !< derivatives of coordinate axis

  REAL   , allocatable :: dRodR(:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dRodZ(:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dZodR(:,:,:) !< derivatives of coordinate axis
  REAL   , allocatable :: dZodZ(:,:,:) !< derivatives of coordinate axis

  INTEGER, allocatable :: djkp(:,:) !< for calculating cylindrical volume
  INTEGER, allocatable :: djkm(:,:) !< for calculating cylindrical volume



  REAL   , allocatable :: lBBintegral(:) !< B.B integral
  REAL   , allocatable :: lABintegral(:) !< A.B integral

  REAL   , allocatable :: vvolume(:) !< volume integral of \f$\sqrt g\f$; computed in volume
  REAL                 :: dvolume    !< derivative of volume w.r.t. interface geometry



  INTEGER              :: ivol !< labels volume; some subroutines (called by NAG) are fixed argument list but require the volume label

  REAL                 :: gBzeta !< toroidal (contravariant) field; calculated in bfield; required to convert \f$\dot \theta\f$ to \f$B^\theta\f$, \f$\dot s\f$ to \f$B^s\f$

  INTEGER, allocatable :: Iquad(:) !< internal copy of Nquad

  REAL   , allocatable :: gaussianweight(:,:)    !<   weights for Gaussian quadrature
  REAL   , allocatable :: gaussianabscissae(:,:) !< abscissae for Gaussian quadrature

  LOGICAL              :: LBlinear !< controls selection of Beltrami field solver; depends on LBeltrami
  LOGICAL              :: LBnewton !< controls selection of Beltrami field solver; depends on LBeltrami
  LOGICAL              :: LBsequad !< controls selection of Beltrami field solver; depends on LBeltrami

  REAL                 :: oRZp(1:3) !< used in mg00aa() to determine \f$(s,\theta,\zeta)\f$ given \f$(R,Z,\varphi)\f$


  type(derivative)     :: dBdX !< \f${\rm d}\mathbf{B}/{\rm d}\mathbf{X}\f$ (?)

  INTEGER              :: globaljk  !< labels position
  REAL, allocatable    :: Dxyz(:,:) !< computational boundary; position
  REAL, allocatable    :: Nxyz(:,:) !< computational boundary; normal
  REAL, allocatable    :: Jxyz(:,:) !< plasma        boundary; surface current

  REAL                 :: tetazeta(1:2) !< what is this?

  REAL                 :: virtualcasingfactor = -one / (four*pi) !< this agrees with diagno

  INTEGER              :: IBerror !< for computing error in magnetic field

  INTEGER              :: nfreeboundaryiterations !< number of free-boundary iterations already performed



  INTEGER, parameter   :: Node = 2 !< best to make this global for consistency between calling and called routines



  LOGICAL              :: first_free_bound = .false. !< flag to indicate that this is the first free-boundary iteration

contains



subroutine build_vector_potential(lvol, iocons, aderiv, tderiv)


  use constants, only: zero, half

  use fileunits, only: ounit

  use inputlist, only: Lrad, Wbuild_vector_potential, Wmacros

  use cputiming



  LOCALS

  INTEGER              :: aderiv    ! Derivative of A. -1: w.r.t geometrical degree of freedom
  INTEGER              :: tderiv    ! Derivative of Chebyshev polynomialc. 0: no derivatives
  INTEGER              :: ii,  &    ! Loop index on Fourier harmonics
                          ll,  &    ! Loop index on radial resolution
                          mi,  &    ! Poloidal mode number
                          lvol,&    ! Volume number
                          iocons    ! inner (0) or outer (1) side of the volume
  REAL                 :: mfactor   ! Regularization factor when LcoordinateSingularity



  BEGIN(build_vector_potential)

  efmn(1:mn) = zero ; sfmn(1:mn) = zero ; cfmn(1:mn) = zero ; ofmn(1:mn) = zero

  do ii = 1, mn ! loop over Fourier harmonics; 13 Sep 13;

   if( Lcoordinatesingularity ) then
    mi = im(ii)
    do ll = mi, Lrad(lvol),2 ! loop over Zernike polynomials; Lrad is the radial resolution; 01 Jul 19;
      ;                      ; efmn(ii) = efmn(ii) +          Ate(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      ;                      ; cfmn(ii) = cfmn(ii) +          Aze(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) +          Ato(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      ;                      ; sfmn(ii) = sfmn(ii) +          Azo(lvol,aderiv,ii)%s(ll) * RTT(ll,mi,iocons,1) * half
      endif
    enddo ! end of do ll; 20 Feb 13;
   else
    do ll = 0, Lrad(lvol) ! loop over Chebyshev polynomials; Lrad is the radial resolution;
      ;                      ; efmn(ii) = efmn(ii) +          Ate(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1) ! aderiv labels deriv. wrt mu, pflux;
      ;                      ; cfmn(ii) = cfmn(ii) +          Aze(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      if( NOTstellsym ) then ; ofmn(ii) = ofmn(ii) +          Ato(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      ;                     ; sfmn(ii) = sfmn(ii) +          Azo(lvol,aderiv,ii)%s(ll) * TT(ll,iocons,1)
      endif
    enddo ! end of do ll; 20 Feb 13;
   end if ! Lcoordinatesingularity; 01 Jul 19;
  enddo ! end of do ii; 20 Feb 13;

end subroutine build_vector_potential



subroutine set_mpi_comm(comm)

  implicit none

  integer, intent(in) :: comm
  integer             :: ierr

  MPI_COMM_SPEC = comm

  myid = 0 ; ncpu = 1

  ierr=0
  call MPI_COMM_RANK( MPI_COMM_SPEC, myid, ierr )
  if (ierr.ne.0) write(*,*) "error in call to MPI_COMM_RANK"

  ierr=0
  call MPI_COMM_SIZE( MPI_COMM_SPEC, ncpu, ierr )
  if (ierr.ne.0) write(*,*) "error in call to MPI_COMM_SIZE"

end subroutine



subroutine read_inputlists_from_file()

   use constants
   use fileunits
   use inputlist

#ifdef IFORT
   use ifport ! for fseek, ftell with Intel compiler
#endif

   LOCALS

   LOGICAL              :: Lspexist
   integer :: filepos, seek_status, cpfile, instat, idx_mode

   character(len=1000) :: line

   INTEGER              :: mm, nn, MNMAX
   REAL,    allocatable :: RZRZ(:,:) ! local array used for reading interface Fourier harmonics from file;

   inquire( file=trim(ext)//".sp", exist=Lspexist ) ! check if file exists;
   FATAL( readin, .not.Lspexist, the input file does not exist ) ! if not, abort;



   open( iunit, file=trim(ext)//".sp", status="old")



   instat = 0 ! initially, no error

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading physicslist     from ext.sp ;")') cput-cpus
   endif

   read(iunit, physicslist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in physicslist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    physicslist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading numericlist     from ext.sp ;")') cput-cpus
   endif

   read(iunit, numericlist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in numericlist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    numericlist     from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading locallist      from ext.sp ;")') cput-cpus
   endif

   read(iunit, locallist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in locallist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    locallist      from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading globallist   from ext.sp ;")') cput-cpus
   endif

   read(iunit, globallist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in globallist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    globallist   from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading diagnosticslist from ext.sp ;")') cput-cpus
   endif

   read(iunit, diagnosticslist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in diagnosticslist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    diagnosticslist from ext.sp ;")') cput-cpus
   endif

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : reading screenlist      from ext.sp ;")') cput-cpus
   endif

   read(iunit, screenlist, iostat=instat)
   if (instat .ne. 0) then
     backspace(iunit)
     read(iunit,fmt='(A)') line
     write(*,'(A)') 'Invalid line in screenlist: '//trim(line)
   end if

   if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : read    screenlist      from ext.sp ;")') cput-cpus
   endif


   instat = 0

   num_modes = 0

   MNMAX = MNtor + 1 + MMpol * ( 2 * MNtor + 1 )
   if(allocated(mmRZRZ)) deallocate(mmRZRZ, nnRZRZ, allRZRZ)
   allocate(mmRZRZ(1:MNMAX), nnRZRZ(1:MNMAX), allRZRZ(1:4,1:Nvol,1:MNMAX))

   if (Linitialize .le. 0) then

     FATAL( readin, Nvol.lt.1 .or. Nvol.gt.MNvol,            invalid Nvol: may need to recompile with higher MNvol )

     SALLOCATE( RZRZ, (1:4,1:Nvol), zero ) ! temp array for reading input;

#ifdef IFORT
     filepos = ftell(iunit)+1
#else
     call ftell(iunit, filepos)
#endif
     do ! will read in Fourier harmonics until the end of file is reached;
       read(iunit,*,iostat=instat) mm, nn, RZRZ(1:4,1:Nvol)   !if change of angle applies, transformation assumes m>=0 and for m=0 only n>=0;
       if( instat.ne.0 ) exit

       num_modes = num_modes + 1
     enddo

     rewind(iunit)

#ifdef IFORT
     seek_status = fseek(iunit, filepos, 0)
#else
     call fseek(iunit, filepos, 0, seek_status)
#endif
     FATAL(inplst, seek_status.ne.0, failed to seek to end of input namelists )

     do idx_mode = 1, num_modes
       read(iunit,*,iostat=instat) mmRZRZ(idx_mode), nnRZRZ(idx_mode), allRZRZ(1:4,1:Nvol, idx_mode)
     enddo

     DALLOCATE(RZRZ)

    end if ! Linitialize .le. 0

    close(iunit)

end subroutine ! read_inputlists_from_file





subroutine write_spec_namelist()
  use constants
  use fileunits
  use inputlist

  LOCALS

  LOGICAL :: exist
  CHARACTER(LEN=100), PARAMETER :: example = 'example.sp'

  if( myid == 0 ) then
     inquire(file=trim(example), EXIST=exist) ! inquire if inputfile existed;
     FATAL( global, exist, example input file example.sp already existed )
     open(iunit, file=trim(example), status='unknown', action='write')
     write(iunit, physicslist)
     write(iunit, numericlist)
     write(iunit, locallist)
     write(iunit, globallist)
     write(iunit, diagnosticslist)
     write(iunit, screenlist)
     close(iunit)
  endif

  return
end subroutine



subroutine check_inputs()

   use numerical
   use constants
   use fileunits
   use inputlist
   use cputiming, only: Treadin

   LOCALS

   INTEGER              :: vvol
   REAL                 :: xx, toroidalflux, toroidalcurrent

   BEGIN(readin)

   Mvol = Nvol + Lfreebound ! this is just for screen output and initial check; true assignment of Mvol appears outside if( myid.eq.0 ) then ;

   write(ounit,'("readin : ", 10x ," : ")')

   cput = GETTIME

   write(ounit,1010) cput-cpus, Igeometry, Istellsym, Lreflect
   write(ounit,1011)            Lfreebound, phiedge, curtor, curpol
   write(ounit,1012)            gamma
   write(ounit,1013)            Nfp, Nvol, Mvol, Mpol, Ntor
   write(ounit,1014)            pscale, Ladiabatic, Lconstraint, mupftol, mupfits
   write(ounit,1015)            Lrad(1:min(Mvol,32))

1010 format("readin : ",f10.2," : Igeometry=",i3," ; Istellsym=",i3," ; Lreflect="i3" ;")
1011 format("readin : ", 10x ," : Lfreebound=",i3," ; phiedge="es23.15" ; curtor="es23.15" ; curpol="es23.15" ;")
1012 format("readin : ", 10x ," : gamma="es23.15" ;")
1013 format("readin : ", 10x ," : Nfp=",i3," ; Nvol=",i3," ; Mvol=",i3," ; Mpol=",i3," ; Ntor=",i3," ;")
1014 format("readin : ", 10x ," : pscale="es13.5" ; Ladiabatic="i2" ; Lconstraint="i3" ; mupf: tol,its="es10.2" ,"i4" ;")
1015 format("readin : ", 10x ," : Lrad = "257(i2,",",:))

   FATAL( readin, Igeometry.lt.1 .or. Igeometry.gt.3, invalid geometry )
   FATAL( readin, Nfp.le.0, invalid Nfp )
   FATAL( readin, Mpol.lt.0 .or. Mpol.gt.MMpol, invalid poloidal resolution: may need to recompile with higher MMpol )
   FATAL( readin, Ntor.lt.0 .or. Ntor.gt.MNtor, invalid toroidal resolution: may need to recompile with higher MNtor )
   FATAL( readin, Nvol.lt.1 .or. Nvol.gt.MNvol, invalid Nvol: may need to recompile with higher MNvol )
   FATAL( readin, mupftol.le.zero, mupftol is too small )
   FATAL( readin, abs(one+gamma).lt.vsmall, 1+gamma appears in denominator in dforce ) !< \todo Please check this; SRH: 27 Feb 18;
   FATAL( readin, abs(one-gamma).lt.vsmall, 1-gamma appears in denominator in fu00aa ) !< \todo Please check this; SRH: 27 Feb 18;
   FATAL( readin, Lconstraint.lt.-1 .or. Lconstraint.gt.3, illegal Lconstraint )
   FATAL( readin, Igeometry.eq.1 .and. rpol.lt.vsmall, poloidal extent of slab too small or negative )
   FATAL( readin, Igeometry.eq.1 .and. rtor.lt.vsmall, toroidal extent of slab too small or negative )

   if( Istellsym.eq.1 ) then
    Rbs(-MNtor:MNtor,-MMpol:MMpol) = zero
    Zbc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Rws(-MNtor:MNtor,-MMpol:MMpol) = zero
    Zwc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Vnc(-MNtor:MNtor,-MMpol:MMpol) = zero
    Bnc(-MNtor:MNtor,-MMpol:MMpol) = zero
   endif




   FATAL( readin, abs(tflux(Nvol)).lt. vsmall, enclosed toroidal flux cannot be zero )

   toroidalflux = tflux(Nvol) ! toroidal flux is a local variable; SRH: 27 Feb 18

   tflux(1:Mvol) = tflux(1:Mvol) / toroidalflux ! normalize toroidal flux
   pflux(1:Mvol) = pflux(1:Mvol) / toroidalflux ! normalize poloidal flux

   FATAL( readin, tflux(1).lt.zero, enclosed toroidal flux cannot be zero )
   do vvol = 2, Mvol
   enddo

   do vvol = 1, Mvol
    FATAL( readin, Lrad(vvol ).lt.2, require Chebyshev resolution Lrad > 2 so that Lagrange constraints can be satisfied )
   enddo

   if (Igeometry.ge.2 .and. Lrad(1).lt.Mpol) then
     write(ounit,'("readin : ",f10.2," : Minimum Lrad(1) is Mpol, automatically adjusted it to Mpol+4")') cput-cpus
     Lrad(1) = Mpol + 4
   endif
   FATAL( readin, mupfits.le.0, must give ma01aa:hybrj a postive integer value for the maximum iterations = mupfits given on input )




   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1020) cput-cpus, Linitialize, LautoinitBn, Lzerovac, Ndiscrete
   write(ounit,1021)            Nquad, iMpol, iNtor
   write(ounit,1022)            Lsparse, Lsvdiota, imethod, iorder, iprecon, iotatol
   write(ounit,1023)            Lextrap, Mregular, Lrzaxis, Ntoraxis

1020 format("readin : ",f10.2," : Linitialize=",i3," ;LautoinitBn=",i3," ; Lzerovac=",i2," ; Ndiscrete="i2" ;")
1021 format("readin : ", 10x ," : Nquad="i4" ; iMpol="i4" ; iNtor="i4" ;")
1022 format("readin : ", 10x ," : Lsparse="i2" ; Lsvdiota="i2" ; imethod="i2" ; iorder="i2" ; iprecon="i2" ; iotatol="es13.5" ;")
1023 format("readin : ", 10x ," : Lextrap="i2" ; Mregular="i3" ; Lrzaxis="i2" ; Ntoraxis="i2" ;")

   FATAL( readin, Ndiscrete.le.0, error )


   FATAL( readin, iotatol.gt.one, illegal value for sparse tolerance ) ! I think that the sparse iota solver is no longer implemented; SRH: 27 Feb 18;





   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1030) cput-cpus, LBeltrami, Linitgues, Lmatsolver, LGMRESprec, NiterGMRES, epsGMRES, epsILU

1030 format("readin : ",f10.2," : LBeltrami="i2" ; Linitgues="i2" ; Lmatsolver="i2" ; LGMRESprec="i2" ; NiterGMRES="i4" ; epsGMRES="es13.5" ; epsILU="es13.5" ;" )

   FATAL( readin, LBeltrami.lt.0 .or. LBeltrami.gt.7, error )
   FATAL( readin, LGMRESprec.lt.0 .or. LGMRESprec.gt.1, error )
   FATAL( readin, NiterGMRES.lt.0, error )
   FATAL( readin, abs(epsGMRES).le.machprec , error )
   FATAL( readin, abs(epsILU).le.machprec , error )




   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1040) cput-cpus, Lfindzero
   write(ounit,1041)            escale, opsilon, pcondense, epsilon, wpoloidal, upsilon
   write(ounit,1042)            forcetol, c05xmax, c05xtol, c05factor, LreadGF
   write(ounit,1043)            mfreeits, gBntol, gBnbld
   write(ounit,1044)            vcasingeps, vcasingtol, vcasingits, vcasingper

1040 format("readin : ",f10.2," : Lfindzero="i2" ;")
1041 format("readin : ", 10x ," : escale="es13.5" ; opsilon="es13.5" ; pcondense="f7.3" ; epsilon="es13.5" ; wpoloidal="f7.4" ; upsilon="es13.5" ;")
1042 format("readin : ", 10x ," : forcetol="es13.5" ; c05xmax="es13.5" ; c05xtol="es13.5" ; c05factor="es13.5" ; LreadGF="L2" ; ")
1043 format("readin : ", 10x ," : mfreeits="i4" ; gBntol="es13.5" ; gBnbld="es13.5" ;")
1044 format("readin : ", 10x ," : vcasingeps="es13.5" ; vcasingtol="es13.5" ; vcasingits="i6" ; vcasingper="i6" ;")

   FATAL( readin, escale      .lt.zero     , error )
   FATAL( readin, pcondense   .lt.one      , error )
   FATAL( readin, abs(c05xtol).le.machprec , error )
   FATAL( readin, c05factor   .le.zero     , error )

   FATAL( readin, Igeometry.eq.3 .and. pcondense.le.zero, pcondense must be positive )




   write(ounit,'("readin : ", 10x ," : ")')

   write(ounit,1050) cput-cpus, odetol, nPpts
   write(ounit,1051)            LHevalues, LHevectors, LHmatrix, Lperturbed, dpp, dqq, dRZ, Lcheck, Ltiming

1050 format("readin : ",f10.2," : odetol="es10.2" ; nPpts="i6" ;")
1051 format("readin : ", 10x ," : LHevalues="L2" ; LHevectors="L2" ; LHmatrix="L2" ; Lperturbed="i2" ; dpp="i3" ; dqq="i3" ; dRZ="es16.8" ; Lcheck="i3" ; Ltiming="L2" ;")

   FATAL( readin, odetol.le.zero, input error )






   write(ounit,'("readin : ", 10x ," : ")')

   RETURN(readin)

end subroutine ! check_inputs



subroutine broadcast_inputs

  use fileunits
  use inputlist

  LOCALS

  ClBCAST( ext        ,     100, 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting physicslist     from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Igeometry  ,       1, 0 )
  IlBCAST( Istellsym  ,       1, 0 )
  IlBCAST( Lfreebound ,       1, 0 )
  RlBCAST( phiedge    ,       1, 0 )
  RlBCAST( curtor     ,       1, 0 )
  RlBCAST( curpol     ,       1, 0 )
  RlBCAST( gamma      ,       1, 0 )
  IlBCAST( Nfp        ,       1, 0 )
  IlBCAST( Nvol       ,       1, 0 )
  IlBCAST( Mpol       ,       1, 0 )
  IlBCAST( Ntor       ,       1, 0 )
  IlBCAST( Lrad       , MNvol+1, 0 )
  RlBCAST( tflux      , MNvol+1, 0 )
  RlBCAST( pflux      , MNvol+1, 0 )
  RlBCAST( helicity   , MNvol  , 0 )
  RlBCAST( pscale     ,       1, 0 )
  RlBCAST( pressure   , MNvol+1, 0 )
  IlBCAST( Ladiabatic ,       1, 0 )
  RlBCAST( adiabatic  , MNvol+1, 0 )
  RlBCAST( mu         , MNvol+1, 0 )
  RlBCAST( Ivolume    , MNvol+1, 0 )
  RlBCAST( Isurf      , MNvol+1, 0 )
  IlBCAST( Lconstraint,       1, 0 )
  IlBCAST( pl         , MNvol  , 0 )
  IlBCAST( ql         , MNvol  , 0 )
  IlBCAST( pr         , MNvol  , 0 )
  IlBCAST( qr         , MNvol  , 0 )
  RlBCAST( iota       , MNvol  , 0 )
  IlBCAST( lp         , MNvol  , 0 )
  IlBCAST( lq         , MNvol  , 0 )
  IlBCAST( rp         , MNvol  , 0 )
  IlBCAST( rq         , MNvol  , 0 )
  RlBCAST( oita       , MNvol  , 0 )
  RlBCAST( mupftol    ,       1, 0 )
  IlBCAST( mupfits    ,       1, 0 )
  IlBCAST( Lreflect   ,       1, 0 )
  RlBCAST( rpol       ,       1, 0 )
  RlBCAST( rtor       ,       1, 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting numericlist     from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Linitialize, 1, 0 )
  IlBCAST( LautoinitBn, 1, 0 )
  IlBCAST( Lzerovac   , 1, 0 )
  IlBCAST( Ndiscrete  , 1, 0 )
  IlBCAST( Nquad      , 1, 0 )
  IlBCAST( iMpol      , 1, 0 )
  IlBCAST( iNtor      , 1, 0 )
  IlBCAST( Lsparse    , 1, 0 )
  IlBCAST( Lsvdiota   , 1, 0 )
  IlBCAST( imethod    , 1, 0 )
  IlBCAST( iorder     , 1, 0 )
  IlBCAST( iprecon    , 1, 0 )
  RlBCAST( iotatol    , 1, 0 )
  IlBCAST( Lextrap    , 1, 0 )
  IlBCAST( Mregular   , 1, 0 )
  IlBCAST( Lrzaxis    , 1, 0 )
  IlBCAST( Ntoraxis   , 1, 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting globallist      from ext.sp ;")') cput-cpus
  endif

  IlBCAST( Lfindzero , 1 , 0 )
  RlBCAST( escale    , 1 , 0 )
  RlBCAST( opsilon   , 1 , 0 )
  RlBCAST( pcondense , 1 , 0 )
  RlBCAST( epsilon   , 1 , 0 )
  RlBCAST( wpoloidal , 1 , 0 )
  RlBCAST( upsilon   , 1 , 0 )
  RlBCAST( forcetol  , 1 , 0 )
  RlBCAST( c05xmax   , 1 , 0 )
  RlBCAST( c05xtol   , 1 , 0 )
  RlBCAST( c05factor , 1 , 0 )
  LlBCAST( LreadGF   , 1 , 0 )
  IlBCAST( mfreeits  , 1 , 0 )
  RlBCAST( gBntol    , 1 , 0 )
  RlBCAST( gBnbld    , 1 , 0 )
  RlBCAST( vcasingeps, 1 , 0 )
  RlBCAST( vcasingtol, 1 , 0 )
  IlBCAST( vcasingits, 1 , 0 )
  IlBCAST( vcasingper, 1 , 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting locallist       from ext.sp ;")') cput-cpus
  endif

  IlBCAST( LBeltrami    , 1, 0 )
  IlBCAST( Linitgues    , 1, 0 )
  RlBCAST( maxrndgues   , 1, 0)
  IlBCAST( Lmatsolver   , 1, 0 )
  IlBCAST( NiterGMRES   , 1, 0 )
  RlBCAST( epsGMRES     , 1, 0 )
  IlBCAST( LGMRESprec   , 1, 0 )
  RlBCAST( epsILU       , 1, 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting diagnosticslist from ext.sp ;")') cput-cpus
  endif

  RlBCAST( odetol    , 1      , 0 )
  IlBCAST( nPpts     , 1      , 0 )
  RlBCAST( Ppts      , 1      , 0 )
  IlBCAST( nPtrj     , MNvol+1, 0 )
  LlBCAST( LHevalues , 1      , 0 )
  LlBCAST( LHevectors, 1      , 0 )
  LlBCAST( Ltransform, 1      , 0 )
  LlBCAST( LHmatrix  , 1      , 0 )
  IlBCAST( Lperturbed, 1      , 0 )
  IlBCAST( dpp       , 1      , 0 )
  IlBCAST( dqq       , 1      , 0 )
  IlBCAST( Lerrortype, 1      , 0 )
  IlBCAST( Ngrid     , 1      , 0 )
  RlBCAST( dRZ       , 1      , 0 )
  IlBCAST( Lcheck    , 1      , 0 )
  LlBCAST( Ltiming   , 1      , 0 )




  if( Wreadin ) then ; cput = GETTIME ; write(ounit,'("readin : ",f10.2," : broadcasting screenlist      from ext.sp ;")') cput-cpus
  endif

  LlBCAST( Wreadin, 1, 0 )
  LlBCAST( Wwrtend, 1, 0 )
  LlBCAST( Wmacros, 1, 0 )

end subroutine ! broadcast_inputs



subroutine wrtend

  use constants, only :

  use numerical, only : machprec

  use fileunits, only : ounit, iunit

  use cputiming, only : Twrtend

  use inputlist



  LOCALS

  INTEGER              :: vvol !< iteration variable over all nested volumes
  INTEGER              :: imn  !< iteration variable for all Fourier harmonics
  INTEGER              :: ii   !< iteration variable for all Fourier harmonics
  INTEGER              :: jj   !< iteration variable
  INTEGER              :: kk   !< iteration variable
  INTEGER              :: jk   !< iteration variable
  INTEGER              :: Lcurvature !< curvature flag (?)
  INTEGER              :: mm   !< current poloidal mode number
  INTEGER              :: nn   !< current toroidal mode number

  REAL                 :: lss !< (?)
  REAL                 :: teta !< (?)
  REAL                 :: zeta !< (?)
  REAL                 :: st(1:Node) !< (?)
  REAL                 :: Bst(1:Node) !< (?)
  REAL                 :: BR !< (?)
  REAL                 :: BZ !< (?)
  REAL                 :: BP !< (?)

  BEGIN(wrtend)



  if( myid.ne.0 ) goto 9999


  open(iunit,file=trim(ext)//".sp.end",status="unknown") ! restart input file;

  write(iunit,'("&physicslist")')
  write(iunit,'(" Igeometry   = ",i9        )') Igeometry
  write(iunit,'(" Istellsym   = ",i9        )') Istellsym
  write(iunit,'(" Lfreebound  = ",i9        )') Lfreebound
  write(iunit,'(" phiedge     = ",es23.15   )') phiedge
  write(iunit,'(" curtor      = ",es23.15   )') curtor
  write(iunit,'(" curpol      = ",es23.15   )') curpol
  write(iunit,'(" gamma       = ",es23.15   )') gamma
  write(iunit,'(" Nfp         = ",i9        )') Nfp
  write(iunit,'(" Nvol        = ",i9        )') Nvol
  write(iunit,'(" Mpol        = ",i9        )') Mpol
  write(iunit,'(" Ntor        = ",i9        )') Ntor
  write(iunit,'(" Lrad        = ",257i23    )') Lrad(1:Mvol)
  write(iunit,'(" tflux       = ",257es23.15)') tflux(1:Mvol)
  write(iunit,'(" pflux       = ",257es23.15)') pflux(1:Mvol)
  write(iunit,'(" helicity    = ",256es23.15)') helicity(1:Mvol)
  write(iunit,'(" pscale      = ",es23.15   )') pscale
  write(iunit,'(" Ladiabatic  = ",i9        )') Ladiabatic
  write(iunit,'(" pressure    = ",257es23.15)') pressure(1:Mvol)
  write(iunit,'(" adiabatic   = ",257es23.15)') adiabatic(1:Mvol)
  write(iunit,'(" mu          = ",257es23.15)') mu(1:Mvol)
  write(iunit,'(" Ivolume     = ",257es23.15)') Ivolume(1:Mvol)
  write(iunit,'(" Isurf       = ",257es23.15)') Isurf(1:Mvol-1), 0.0
  write(iunit,'(" Lconstraint = ",i9        )') Lconstraint
  write(iunit,'(" pl          = ",257i23    )') pl(0:Mvol)
  write(iunit,'(" ql          = ",257i23    )') ql(0:Mvol)
  write(iunit,'(" pr          = ",257i23    )') pr(0:Mvol)
  write(iunit,'(" qr          = ",257i23    )') qr(0:Mvol)
  write(iunit,'(" iota        = ",257es23.15)') iota(0:Mvol)
  write(iunit,'(" lp          = ",257i23    )') lp(0:Mvol)
  write(iunit,'(" lq          = ",257i23    )') lq(0:Mvol)
  write(iunit,'(" rp          = ",257i23    )') rp(0:Mvol)
  write(iunit,'(" rq          = ",257i23    )') rq(0:Mvol)
  write(iunit,'(" oita        = ",257es23.15)') oita(0:Mvol)
  write(iunit,'(" mupftol     = ",es23.15   )') mupftol
  write(iunit,'(" mupfits     = ",i9        )') mupfits
  write(iunit,'(" Lreflect    = ",i9        )') Lreflect
  write(iunit,'(" rpol        = ",es23.15   )') rpol
  write(iunit,'(" rtor        = ",es23.15   )') rtor

 write(iunit,'(" Rac         = ",99es23.15)') iRbc(1:Ntor+1,0)
 write(iunit,'(" Zas         = ",99es23.15)') iZbs(1:Ntor+1,0)
 write(iunit,'(" Ras         = ",99es23.15)') iRbs(1:Ntor+1,0)
 write(iunit,'(" Zac         = ",99es23.15)') iZbc(1:Ntor+1,0)

  do mm = 0, Mpol ! will write out the plasma boundary harmonics;
   do nn = -Ntor, Ntor

    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;

    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1000) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1002) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1001) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 )     write(iunit,1003) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 )     write(iunit,1005) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 )     write(iunit,1004) nn, mm, Rbc(nn,mm), nn, mm, Zbs(nn,mm), nn, mm, Rbs(nn,mm), nn, mm, Zbc(nn,mm)
    end select ! end of select case( mm );

   enddo ! end of do nn;
  enddo ! end of do mm;

  do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
   do nn = -Ntor, Ntor

    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;

    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1010) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1012) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1011) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1013) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1015) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1014) nn, mm, Rwc(nn,mm), nn, mm, Zws(nn,mm), nn, mm, Rws(nn,mm), nn, mm, Zwc(nn,mm)
    end select ! end of select case( mm );

   enddo ! end of do nn;
  enddo ! end of do mm;

  do mm = 0, Mpol ! will write out the computation domain harmonics; (only relevant in free-boundary case);
   do nn = -Ntor, Ntor

    if( mm.eq.0 .and. nn.lt.0 ) cycle ! these modes are always excluded; 13 Oct 12;

    select case( mm )
    case(   0:  9 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1020) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1022) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1021) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
    case(  10: 99 )
     if( nn.lt.- 9 .and. nn.gt.-99 ) write(iunit,1023) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.lt.  0 .and. nn.ge.- 9 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.ge.  0 .and. nn.le.  9 ) write(iunit,1025) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
     if( nn.gt.  9 .and. nn.le. 99 ) write(iunit,1024) nn, mm, Vns(nn,mm), nn, mm, Bns(nn,mm), nn, mm, Vnc(nn,mm), nn, mm, Bnc(nn,mm)
    end select ! end of select case( mm );

   enddo ! end of do nn;
  enddo ! end of do mm;

1000 format("Rbc(",i3,",",i1,")",2x,"=",es23.15," Zbs(",i3,",",i1,")",2x,"=",es23.15," Rbs(",i3,",",i1,")",2x,"=",es23.15," Zbc(",i3,",",i1,")",2x,"=",es23.15)
1001 format("Rbc(",i2,",",i1,")",3x,"=",es23.15," Zbs(",i2,",",i1,")",3x,"=",es23.15," Rbs(",i2,",",i1,")",3x,"=",es23.15," Zbc(",i2,",",i1,")",3x,"=",es23.15)
1002 format("Rbc(",i1,",",i1,")",4x,"=",es23.15," Zbs(",i1,",",i1,")",4x,"=",es23.15," Rbs(",i1,",",i1,")",4x,"=",es23.15," Zbc(",i1,",",i1,")",4x,"=",es23.15)
1003 format("Rbc(",i3,",",i2,")",1x,"=",es23.15," Zbs(",i3,",",i2,")",1x,"=",es23.15," Rbs(",i3,",",i2,")",1x,"=",es23.15," Zbc(",i3,",",i2,")",1x,"=",es23.15)
1004 format("Rbc(",i2,",",i2,")",2x,"=",es23.15," Zbs(",i2,",",i2,")",2x,"=",es23.15," Rbs(",i2,",",i2,")",2x,"=",es23.15," Zbc(",i2,",",i2,")",2x,"=",es23.15)
1005 format("Rbc(",i1,",",i2,")",3x,"=",es23.15," Zbs(",i1,",",i2,")",3x,"=",es23.15," Rbs(",i1,",",i2,")",3x,"=",es23.15," Zbc(",i1,",",i2,")",3x,"=",es23.15)

1010 format("Rwc(",i3,",",i1,")",2x,"=",es23.15," Zws(",i3,",",i1,")",2x,"=",es23.15," Rws(",i3,",",i1,")",2x,"=",es23.15," Zwc(",i3,",",i1,")",2x,"=",es23.15)
1011 format("Rwc(",i2,",",i1,")",3x,"=",es23.15," Zws(",i2,",",i1,")",3x,"=",es23.15," Rws(",i2,",",i1,")",3x,"=",es23.15," Zwc(",i2,",",i1,")",3x,"=",es23.15)
1012 format("Rwc(",i1,",",i1,")",4x,"=",es23.15," Zws(",i1,",",i1,")",4x,"=",es23.15," Rws(",i1,",",i1,")",4x,"=",es23.15," Zwc(",i1,",",i1,")",4x,"=",es23.15)
1013 format("Rwc(",i3,",",i2,")",1x,"=",es23.15," Zws(",i3,",",i2,")",1x,"=",es23.15," Rws(",i3,",",i2,")",1x,"=",es23.15," Zwc(",i3,",",i2,")",1x,"=",es23.15)
1014 format("Rwc(",i2,",",i2,")",2x,"=",es23.15," Zws(",i2,",",i2,")",2x,"=",es23.15," Rws(",i2,",",i2,")",2x,"=",es23.15," Zwc(",i2,",",i2,")",2x,"=",es23.15)
1015 format("Rwc(",i1,",",i2,")",3x,"=",es23.15," Zws(",i1,",",i2,")",3x,"=",es23.15," Rws(",i1,",",i2,")",3x,"=",es23.15," Zwc(",i1,",",i2,")",3x,"=",es23.15)

1020 format("Vns(",i3,",",i1,")",2x,"=",es23.15," Bns(",i3,",",i1,")",2x,"=",es23.15," Vnc(",i3,",",i1,")",2x,"=",es23.15," Bnc(",i3,",",i1,")",2x,"=",es23.15)
1021 format("Vns(",i2,",",i1,")",3x,"=",es23.15," Bns(",i2,",",i1,")",3x,"=",es23.15," Vnc(",i2,",",i1,")",3x,"=",es23.15," Bnc(",i2,",",i1,")",3x,"=",es23.15)
1022 format("Vns(",i1,",",i1,")",4x,"=",es23.15," Bns(",i1,",",i1,")",4x,"=",es23.15," Vnc(",i1,",",i1,")",4x,"=",es23.15," Bnc(",i1,",",i1,")",4x,"=",es23.15)
1023 format("Vns(",i3,",",i2,")",1x,"=",es23.15," Bns(",i3,",",i2,")",1x,"=",es23.15," Vnc(",i3,",",i2,")",1x,"=",es23.15," Bnc(",i3,",",i2,")",1x,"=",es23.15)
1024 format("Vns(",i2,",",i2,")",2x,"=",es23.15," Bns(",i2,",",i2,")",2x,"=",es23.15," Vnc(",i2,",",i2,")",2x,"=",es23.15," Bnc(",i2,",",i2,")",2x,"=",es23.15)
1025 format("Vns(",i1,",",i2,")",3x,"=",es23.15," Bns(",i1,",",i2,")",3x,"=",es23.15," Vnc(",i1,",",i2,")",3x,"=",es23.15," Bnc(",i1,",",i2,")",3x,"=",es23.15)

  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing numericlist ;")') cput-cpus, myid
  endif

  write(iunit,'("&numericlist")')
  write(iunit,'(" Linitialize = ",i9            )') Linitialize
  write(iunit,'(" LautoinitBn = ",i9            )') LautoinitBn
  write(iunit,'(" Lzerovac    = ",i9            )') Lzerovac
  write(iunit,'(" Ndiscrete   = ",i9            )') Ndiscrete
  write(iunit,'(" Nquad       = ",i9            )') Nquad
  write(iunit,'(" iMpol       = ",i9            )') iMpol
  write(iunit,'(" iNtor       = ",i9            )') iNtor
  write(iunit,'(" Lsparse     = ",i9            )') Lsparse
  write(iunit,'(" Lsvdiota    = ",i9            )') Lsvdiota
  write(iunit,'(" imethod     = ",i9            )') imethod
  write(iunit,'(" iorder      = ",i9            )') iorder
  write(iunit,'(" iprecon     = ",i9            )') iprecon
  write(iunit,'(" iotatol     = ",es23.15       )') iotatol
  write(iunit,'(" Lextrap     = ",i9            )') Lextrap
  write(iunit,'(" Mregular    = ",i9            )') Mregular
  write(iunit,'(" Lrzaxis     = ",i9            )') Lrzaxis
  write(iunit,'(" Ntoraxis    = ",i9            )') Ntoraxis
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing locallist ;")') cput-cpus, myid
  endif

  write(iunit,'("&locallist")')
  write(iunit,'(" LBeltrami   = ",i9            )') LBeltrami
  write(iunit,'(" Linitgues   = ",i9            )') Linitgues
  write(iunit,'(" Lmatsolver  = ",i9            )') Lmatsolver
  write(iunit,'(" NiterGMRES  = ",i9            )') NiterGMRES
  write(iunit,'(" LGMRESprec  = ",i9            )') LGMRESprec
  write(iunit,'(" epsGMRES    = ",es23.15       )') epsGMRES
  write(iunit,'(" epsILU      = ",es23.15       )') epsILU

  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing globallist ;")') cput-cpus, myid
  endif

  write(iunit,'("&globallist")')
  write(iunit,'(" Lfindzero   = ",i9            )') Lfindzero
  write(iunit,'(" escale      = ",es23.15       )') escale
  write(iunit,'(" opsilon     = ",es23.15       )') opsilon
  write(iunit,'(" pcondense   = ",es23.15       )') pcondense
  write(iunit,'(" epsilon     = ",es23.15       )') epsilon
  write(iunit,'(" wpoloidal   = ",es23.15       )') wpoloidal
  write(iunit,'(" upsilon     = ",es23.15       )') upsilon
  write(iunit,'(" forcetol    = ",es23.15       )') forcetol
  write(iunit,'(" c05xmax     = ",es23.15       )') c05xmax
  write(iunit,'(" c05xtol     = ",es23.15       )') c05xtol
  write(iunit,'(" c05factor   = ",es23.15       )') c05factor
  write(iunit,'(" LreadGF     = ",L9            )') LreadGF
  write(iunit,'(" mfreeits    = ",i9            )') mfreeits
  write(iunit,'(" gBntol      = ",es23.15       )') gBntol
  write(iunit,'(" gBnbld      = ",es23.15       )') gBnbld
  write(iunit,'(" vcasingeps  = ",es23.15       )') vcasingeps
  write(iunit,'(" vcasingtol  = ",es23.15       )') vcasingtol
  write(iunit,'(" vcasingits  = ",i9            )') vcasingits
  write(iunit,'(" vcasingper  = ",i9            )') vcasingper
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing diagnosticslist ;")') cput-cpus, myid
  endif

  write(iunit,'("&diagnosticslist")')
  write(iunit,'(" odetol      = ",es23.15       )') odetol
  write(iunit,'(" nPpts       = ",i9            )') nPpts
  write(iunit,'(" Ppts        = ",es23.15       )') Ppts
  write(iunit,'(" nPtrj       = ",256i6         )') nPtrj(1:Mvol)
  write(iunit,'(" LHevalues   = ",L9            )') LHevalues
  write(iunit,'(" LHevectors  = ",L9            )') LHevectors
  write(iunit,'(" LHmatrix    = ",L9            )') LHmatrix
  write(iunit,'(" Lperturbed  = ",i9            )') Lperturbed
  write(iunit,'(" dpp         = ",i9            )') dpp
  write(iunit,'(" dqq         = ",i9            )') dqq
  write(iunit,'(" dRZ         = ",es23.15       )') dRZ
  write(iunit,'(" Lcheck      = ",i9            )') Lcheck
  write(iunit,'(" Ltiming     = ",L9            )') Ltiming
  write(iunit,'("/")')

  if( Wwrtend ) then ; cput = GETTIME ; write(ounit,'("wrtend : ",f10.2," : myid=",i3," ; writing screenlist ;")') cput-cpus, myid
  endif

  write(iunit,'("&screenlist")')
  if( Wreadin           ) write(iunit,'(" Wreadin = ",L1                )') Wreadin
  if( Wwrtend           ) write(iunit,'(" Wwrtend = ",L1                )') Wwrtend
  if( Wmacros           ) write(iunit,'(" Wmacros = ",L1                )') Wmacros
  write(iunit,'("/")')

  do imn = 1, mn ; write(iunit,'(2i6,1024es23.15)') im(imn), in(imn)/Nfp, ( iRbc(imn,vvol), iZbs(imn,vvol), iRbs(imn,vvol), iZbc(imn,vvol), vvol = 1, Nvol )
  enddo

  close(iunit)


  RETURN(wrtend)



end subroutine wrtend



subroutine IsMyVolume(vvol)

LOCALS

INTEGER, intent(in) :: vvol



IsMyVolumeValue = -1 ! Error value - Problem with vvol / id
if( myid.ne.modulo(vvol-1,ncpu) ) then
  IsMyVolumeValue = 0
else
  IsMyVolumeValue = 1
endif

end subroutine IsMyVolume



subroutine WhichCpuID(vvol, cpu_id)

LOCALS

INTEGER            :: vvol, cpu_id



cpu_id = modulo(vvol-1,ncpu)

end subroutine WhichCpuID



end module allglobal



module fftw_interface ! JAB; 25 Jul 17

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  TYPE(C_PTR)                            :: planf        !< FFTW-related (?)
  TYPE(C_PTR)                            :: planb        !< FFTW-related (?)
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxin(:,:,:)  !< FFTW-related (?)
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxout(:,:,:) !< FFTW-related (?)

end module fftw_interface


