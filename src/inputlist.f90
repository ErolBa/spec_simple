
module inputlist

    implicit none

    integer, parameter :: MNvol = 256 !< The maximum value of \c Nvol is \c MNvol=256.
    integer, parameter :: MMpol = 128 !< The maximum value of \c Mpol is \c MNpol=64.
    integer, parameter :: MNtor = 128 !< The maximum value of \c Ntor is \c MNtor=64.

    integer :: Igeometry = 3 !< selects Cartesian, cylindrical or toroidal geometry;
    integer :: Istellsym = 1 !< stellarator symmetry is enforced if \c Istellsym==1
    integer :: Lfreebound = 0 !< compute vacuum field surrounding plasma
    real(8) :: phiedge = 1.0 !< total enclosed toroidal magnetic flux;
    real(8) :: curtor = 0.0 !< total enclosed (toroidal) plasma current;
    real(8) :: curpol = 0.0 !< total enclosed (poloidal) linking current;
    real(8) :: gamma = 0.0 !< adiabatic index; cannot set \f$|\gamma| = 1\f$
    integer :: Nfp = 1 !< field periodicity
    integer :: Nvol = 1 !< number of volumes
    integer :: Mpol = 0 !< number of poloidal Fourier harmonics
    integer :: Ntor = 0 !< number of toroidal Fourier harmonics
    integer :: Lrad(1:MNvol + 1) = 4 !< Chebyshev resolution in each volume
    integer :: Lconstraint = -1 !< selects constraints; primarily used in ma02aa() and mp00ac().
    real(8) :: tflux(1:MNvol + 1) = 0.0 !< toroidal flux, \f$\psi_t\f$, enclosed by each interface
    real(8) :: pflux(1:MNvol + 1) = 0.0 !< poloidal flux, \f$\psi_p\f$, enclosed by each interface
    real(8) :: helicity(1:MNvol) = 0.0 !< helicity, \f${\cal K}\f$, in each volume, \f${\cal V}_i\f$
    real(8) :: pscale = 0.0 !< pressure scale factor
    real(8) :: pressure(1:MNvol + 1) = 0.0 !< pressure in each volume
    integer :: Ladiabatic = 0 !< logical flag
    real(8) :: adiabatic(1:MNvol + 1) = 0.0 !< adiabatic constants in each volume
    real(8) :: mu(1:MNvol + 1) = 0.0 !< helicity-multiplier, \f$\mu\f$, in each volume
    real(8) :: Ivolume(1:MNvol + 1) = 0.0 !< Toroidal current constraint normalized by \f$\mu_0\f$ (\f$I_{volume} = \mu_0\cdot [A]\f$), in each volume.
    real(8) :: Isurf(1:MNvol) = 0.0 !< Toroidal current normalized by \f$\mu_0\f$ at each interface (cumulative). This is the sum of all pressure driven currents.
    integer :: pl(0:MNvol) = 0 !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: ql(0:MNvol) = 0 !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: pr(0:MNvol) = 0 !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: qr(0:MNvol) = 0 !< "inside" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    real(8) :: iota(0:MNvol) = 0.0 !< rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on inner side of each interface
    integer :: lp(0:MNvol) = 0 !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: lq(0:MNvol) = 0 !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: rp(0:MNvol) = 0 !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    integer :: rq(0:MNvol) = 0 !< "outer" interface rotational-transform is \f$\mbox{$\,\iota\!\!$-} = (p_l+\gamma p_r)/(q_l+\gamma q_r)\f$,
    real(8) :: oita(0:MNvol) = 0.0 !< rotational-transform, \f$\mbox{$\,\iota\!\!$-}\f$, on outer side of each interface
    real(8) :: mupftol = 1.0e-14 !< accuracy to which \f$\mu\f$ and \f$\Delta\psi_p\f$ are required
    integer :: mupfits = 8 !< an upper limit on the transform/helicity constraint iterations;
    real(8) :: rpol = 1.0 !< poloidal extent of slab (effective radius)
    real(8) :: rtor = 1.0 !< toroidal extent of slab (effective radius)
    integer :: Lreflect = 0 !< =1 reflect the upper and lower bound in slab, =0 do not reflect

    real(8) :: Rac(0:MNtor) = 0.0 !<     stellarator symmetric coordinate axis;
    real(8) :: Zas(0:MNtor) = 0.0 !<     stellarator symmetric coordinate axis;
    real(8) :: Ras(0:MNtor) = 0.0 !< non-stellarator symmetric coordinate axis;
    real(8) :: Zac(0:MNtor) = 0.0 !< non-stellarator symmetric coordinate axis;

    real(8) :: Rbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric boundary components;
    real(8) :: Zbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric boundary components;
    real(8) :: Rbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric boundary components;
    real(8) :: Zbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric boundary components;

    real(8) :: Rwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric boundary components of wall;
    real(8) :: Zws(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric boundary components of wall;
    real(8) :: Rws(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric boundary components of wall;
    real(8) :: Zwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric boundary components of wall;

    real(8) :: Vns(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric normal field at boundary; vacuum component;
    real(8) :: Bns(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !<     stellarator symmetric normal field at boundary; plasma component;
    real(8) :: Vnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric normal field at boundary; vacuum component;
    real(8) :: Bnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0 !< non-stellarator symmetric normal field at boundary; plasma component;

    integer :: Linitialize = 0 !< Used to initialize geometry using a regularization / extrapolation method
    integer :: LautoinitBn = 1 !< Used to initialize \f$B_{ns}\f$ using an initial fixed-boundary calculation
    integer :: Lzerovac = 0 !< Used to adjust vacuum field to cancel plasma field on computational boundary
    integer :: Ndiscrete = 2 !< resolution of the real space grid on which fast Fourier transforms are performed is given by \c Ndiscrete*Mpol*4
    integer :: Nquad = -1 !< Resolution of the Gaussian quadrature
    integer :: iMpol = -4 !< Fourier resolution of straight-fieldline angle on interfaces
    integer :: iNtor = -4 !< Fourier resolution of straight-fieldline angle on interfaces;
    integer :: Lsparse = 0 !< controls method used to solve for rotational-transform on interfaces
    integer :: Lsvdiota = 0 !< controls method used to solve for rotational-transform on interfaces;
    integer :: imethod = 3 !< controls iterative solution to sparse matrix
    integer :: iorder = 2 !< controls real-space grid resolution for constructing the straight-fieldline angle;
    integer :: iprecon = 0 !< controls iterative solution to sparse matrix arising in real-space transformation
    real(8) :: iotatol = -1.0 !< tolerance required for iterative construction of straight-fieldline angle;
    integer :: Lextrap = 0 !< geometry of innermost interface is defined by extrapolation
    integer :: Mregular = -1 !< maximum regularization factor
    integer :: Lrzaxis = 1 !< controls the guess of geometry axis in the innermost volume or initialization of interfaces
    integer :: Ntoraxis = 3 !< the number of \f$n\f$ harmonics used in the Jacobian \f$m=1\f$ harmonic elimination method;

    integer :: LBeltrami = 4 !< Control flag for solution of Beltrami equation
    integer :: Linitgues = 1 !< controls how initial guess for Beltrami field is constructed
    integer :: Lposdef = 0 !< redundant;
    real(8) :: maxrndgues = 1.0 !< the maximum random number of the Beltrami field if \c Linitgues = 3
    integer :: Lmatsolver = 3 !< 1 for LU factorization, 2 for GMRES, 3 for GMRES matrix-free
    integer :: NiterGMRES = 200 !< number of max iteration for GMRES
    real(8) :: epsGMRES = 1e-14 !< the precision of GMRES
    integer :: LGMRESprec = 1 !< type of preconditioner for GMRES, 1 for ILU sparse matrix
    real(8) :: epsILU = 1e-12 !< the precision of incomplete LU factorization for preconditioning

    integer :: Lfindzero = 0 !< use Newton methods to find zero of force-balance, which is computed by dforce()
    real(8) :: escale = 0.0 !< controls the weight factor, \c BBweight, in the force-imbalance harmonics
    real(8) :: opsilon = 1.0 !< weighting of force-imbalance
    real(8) :: pcondense = 2.0 !< spectral condensation parameter
    real(8) :: epsilon = 0.0 !< weighting of spectral-width constraint
    real(8) :: wpoloidal = 1.0 !< "star-like" poloidal angle constraint radial exponential factor
    real(8) :: upsilon = 1.0 !< weighting of "star-like" poloidal angle constraint
    real(8) :: forcetol = 1.0e-10 !< required tolerance in force-balance error; only used as an initial check
    real(8) :: c05xmax = 1.0e-06 !< required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$
    real(8) :: c05xtol = 1.0e-12 !< required tolerance in position, \f${\bf x} \equiv \{ R_{i,v}, Z_{i,v}\}\f$
    real(8) :: c05factor = 1.0e-02 !< used to control initial step size in
    logical :: LreadGF = .true. !< read \f$\nabla_{\bf x} {\bf F}\f$ from file \c ext.GF
    integer :: mfreeits = 0 !< maximum allowed free-boundary iterations
    real(8) :: bnstol = 1.0e-06 !< redundant;
    real(8) :: bnsblend = 0.666 !< redundant;
    real(8) :: gBntol = 1.0e-06 !< required tolerance in free-boundary iterations
    real(8) :: gBnbld = 0.666 !< normal blend
    real(8) :: vcasingeps = 1.e-12 !< regularization of Biot-Savart; see bnorml(), casing()
    real(8) :: vcasingtol = 1.e-08 !< accuracy on virtual casing integral; see bnorml(), casing()
    integer :: vcasingits = 8 !< minimum number of calls to adaptive virtual casing routine; see casing()
    integer :: vcasingper = 1 !< periods of integragion  in adaptive virtual casing routine; see casing()
    integer :: mcasingcal = 8 !< minimum number of calls to adaptive virtual casing routine; see casing(); redundant;

    real(8) :: odetol = 1.0e-07 !< o.d.e. integration tolerance for all field line tracing routines
    real(8) :: absreq = 1.0e-08 !< redundant
    real(8) :: relreq = 1.0e-08 !< redundant
    real(8) :: absacc = 1.0e-04 !< redundant
    real(8) :: epsr = 1.0e-08 !< redundant
    integer :: nPpts = 0 !< number of toroidal transits used (per trajectory) in following field lines
    real(8) :: Ppts = 0.0 !< stands for Poincare plot theta start. Chose at which angle (normalized over \f$\pi\f$) the Poincare field-line tracing start.
    integer :: nPtrj(1:MNvol + 1) = -1 !< number of trajectories in each annulus to be followed in constructing Poincaré plot
    logical :: LHevalues = .false. !< to compute eigenvalues of \f$\nabla {\bf F}\f$
    logical :: LHevectors = .false. !< to compute eigenvectors (and also eigenvalues) of \f$\nabla {\bf F}\f$
    logical :: LHmatrix = .false. !< to compute and write to file the elements of \f$\nabla {\bf F}\f$
    integer :: Lperturbed = 0 !< to compute linear, perturbed equilibrium
    integer :: dpp = -1 !< perturbed harmonic
    integer :: dqq = -1 !< perturbed harmonic
    integer :: Lerrortype = 0 !< the type of error output for Lcheck=1
    integer :: Ngrid = -1 !< the number of points to output in the grid, -1 for Lrad(vvol)
    real(8) :: dRZ = 1e-5 !< difference in geometry for finite difference estimate (debug only)
    integer :: Lcheck = 0 !< implement various checks
    logical :: Ltiming = .false. !< to check timing
    logical :: Ltransform = .false. !< to evaluate iota and straight field line coordinates
    real(8) :: fudge = 1.0e-00 !< redundant
    real(8) :: scaling = 1.0e-00 !< redundant

    logical :: Wmanual = .false.
    logical :: Wrzaxis = .false.
    logical :: Wpackxi = .false.
    logical :: Wvolume = .false.
    logical :: Wcoords = .false.
    logical :: Wbasefn = .false.
    logical :: Wmemory = .false.
    logical :: Wmetrix = .false.
    logical :: Wma00aa = .false.
    logical :: Wmatrix = .false.
    logical :: Wmp00ac = .false.
    logical :: Wma02aa = .false.
    logical :: Wpackab = .false.
    logical :: Wtr00ab = .false.
    logical :: Wcurent = .false.
    logical :: Wlforce = .false.
    logical :: Wintghs = .false.
    logical :: Wmtrxhs = .false.
    logical :: Wlbpol = .false.
    logical :: Wbrcast = .false.
    logical :: Wdfp100 = .false.
    logical :: Wdfp200 = .false.
    logical :: Wdforce = .false.
    logical :: Wnewton = .false.
    logical :: Wcasing = .false.
    logical :: Wjo00aa = .false.
    logical :: Wpp00aa = .false.
    logical :: Wpp00ab = .false.
    logical :: Wbfield = .false.
    logical :: Wstzxyz = .false.
    logical :: Whesian = .false.
    logical :: Wra00aa = .false.
    logical :: Wnumrec = .false.
    logical :: Wdcuhre = .false.
    logical :: Wminpack = .false.
    logical :: Wiqpack = .false.
    logical :: Wrksuite = .false.
    logical :: Wi1mach = .false.
    logical :: Wd1mach = .false.
    logical :: Wilut = .false.
    logical :: Witers = .false.
    logical :: Wsphdf5 = .false.
    logical :: Wpreset = .false.
    logical :: Wglobal = .false.
    logical :: Wxspech = .false.
    logical :: Wbuild_vector_potential = .false. !< \todo: what is this?
    logical :: Wreadin = .false. !< write screen output of readin()
    logical :: Wwrtend = .false. !< write screen output of wrtend()
    logical :: Wmacros = .false. !< write screen output from expanded macros

    namelist /physicslist/ &
        Igeometry, &
        Istellsym, &
        Lfreebound, &
        phiedge, &
        curtor, &
        curpol, &
        gamma, &
        Nfp, &
        Nvol, &
        Mpol, &
        Ntor, &
        Lrad, &
        Lconstraint, &
        tflux, &
        pflux, &
        helicity, &
        pscale, &
        pressure, &
        Ladiabatic, &
        adiabatic, &
        mu, &
        Ivolume, &
        Isurf, &
        pl, &
        ql, &
        pr, &
        qr, &
        iota, &
        lp, &
        lq, &
        rp, &
        rq, &
        oita, &
        mupftol, &
        mupfits, &
        rpol, &
        rtor, &
        Lreflect, &
        Rac, &
        Zas, &
        Ras, &
        Zac, &
        Rbc, &
        Zbs, &
        Rbs, &
        Zbc, &
        Rwc, &
        Zws, &
        Rws, &
        Zwc, &
        Vns, &
        Bns, &
        Vnc, &
        Bnc

    namelist /numericlist/ &
        Linitialize, &
        LautoinitBn, &
        Lzerovac, &
        Ndiscrete, &
        Nquad, &
        iMpol, &
        iNtor, &
        Lsparse, &
        Lsvdiota, &
        imethod, &
        iorder, &
        iprecon, &
        iotatol, &
        Lextrap, &
        Mregular, &
        Lrzaxis, &
        Ntoraxis

    namelist /locallist/ &
        LBeltrami, &
        Linitgues, &
        maxrndgues, &
        maxrndgues, &
        Lmatsolver, &
        NiterGMRES, &
        epsGMRES, &
        LGMRESprec, &
        epsILU, &
        Lposdef

    namelist /globallist/ &
        Lfindzero, &
        escale, &
        opsilon, &
        pcondense, &
        epsilon, &
        wpoloidal, &
        upsilon, &
        forcetol, &
        c05xmax, &
        c05xtol, &
        c05factor, &
        LreadGF, &
        mfreeits, &
        bnstol, &
        bnsblend, &
        gBntol, &
        gBnbld, &
        vcasingeps, &
        vcasingtol, &
        vcasingits, &
        vcasingper, &
        mcasingcal

    namelist /diagnosticslist/ &
        odetol, &
        absreq, &
        relreq, &
        absacc, &
        epsr, &
        nPpts, &
        Ppts, &
        nPtrj, &
        LHevalues, &
        LHevectors, &
        LHmatrix, &
        Lperturbed, &
        dpp, &
        dqq, &
        Lerrortype, &
        Ngrid, &
        Lcheck, &
        dRZ, &
        Ltiming, &
        Ltransform, &
        fudge, &
        scaling

    namelist /screenlist/ &
        Wmanual, &
        Wrzaxis, &
        Wpackxi, &
        Wvolume, &
        Wcoords, &
        Wbasefn, &
        Wmemory, &
        Wmetrix, &
        Wma00aa, &
        Wmatrix, &
        Wmp00ac, &
        Wma02aa, &
        Wpackab, &
        Wtr00ab, &
        Wcurent, &
        Wlforce, &
        Wintghs, &
        Wmtrxhs, &
        Wlbpol, &
        Wbrcast, &
        Wdfp100, &
        Wdfp200, &
        Wdforce, &
        Wnewton, &
        Wcasing, &
        Wjo00aa, &
        Wpp00aa, &
        Wpp00ab, &
        Wbfield, &
        Wstzxyz, &
        Whesian, &
        Wra00aa, &
        Wnumrec, &
        Wdcuhre, &
        Wminpack, &
        Wiqpack, &
        Wrksuite, &
        Wi1mach, &
        Wd1mach, &
        Wilut, &
        Witers, &
        Wsphdf5, &
        Wpreset, &
        Wglobal, &
        Wxspech, &
        Wbuild_vector_potential, &
        Wreadin, &
        Wwrtend, &
        Wmacros

contains

    subroutine initialize_inputs

        implicit none

        Igeometry = 3
        Istellsym = 1
        Lfreebound = 0
        phiedge = 1.0
        curtor = 0.0
        curpol = 0.0
        gamma = 0.0
        Nfp = 1
        Nvol = 1
        Mpol = 0
        Ntor = 0
        Lrad(1:MNvol + 1) = 4
        Lconstraint = -1
        tflux(1:MNvol + 1) = 0.0
        pflux(1:MNvol + 1) = 0.0
        helicity(1:MNvol) = 0.0
        pscale = 0.0
        pressure(1:MNvol + 1) = 0.0
        Ladiabatic = 0
        adiabatic(1:MNvol + 1) = 0.0
        mu(1:MNvol + 1) = 0.0
        Ivolume(1:MNvol + 1) = 0.0
        Isurf(1:MNvol) = 0.0
        pl(0:MNvol) = 0
        ql(0:MNvol) = 0
        pr(0:MNvol) = 0
        qr(0:MNvol) = 0
        iota(0:MNvol) = 0.0
        lp(0:MNvol) = 0
        lq(0:MNvol) = 0
        rp(0:MNvol) = 0
        rq(0:MNvol) = 0
        oita(0:MNvol) = 0.0
        rpol = 1.0
        rtor = 1.0

        Rac(0:MNtor) = 0.0
        Zas(0:MNtor) = 0.0
        Ras(0:MNtor) = 0.0
        Zac(0:MNtor) = 0.0

        Rbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Zbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Rbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Zbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

        Rwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Zws(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Rws(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Zwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

        Vns(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Bns(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Vnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
        Bnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

        mupftol = 1.0e-14
        mupfits = 8

        Lreflect = 0

        Linitialize = 0
        LautoinitBn = 1
        Lzerovac = 0
        Ndiscrete = 2
        Nquad = -1
        iMpol = -4
        iNtor = -4
        Lsparse = 0
        Lsvdiota = 0
        imethod = 3
        iorder = 2
        iprecon = 0
        iotatol = -1.0
        Lextrap = 0
        Mregular = -1
        Lrzaxis = 1
        Ntoraxis = 3

        LBeltrami = 4
        Linitgues = 1
        Lposdef = 0
        maxrndgues = 1.0
        Lmatsolver = 1
        NiterGMRES = 200
        epsGMRES = 1e-14
        LGMRESprec = 1
        epsILU = 1e-12

        Lfindzero = 0
        escale = 0.0
        opsilon = 1.0
        pcondense = 2.0
        epsilon = 0.0
        wpoloidal = 1.0
        upsilon = 1.0
        forcetol = 1.0e-10
        c05xmax = 1.0e-06
        c05xtol = 1.0e-12
        c05factor = 1.0e-02
        LreadGF = .true.
        mfreeits = 0
        bnstol = 1.0e-06
        bnsblend = 0.666
        gBntol = 1.0e-06
        gBnbld = 0.666
        vcasingeps = 1.e-12
        vcasingtol = 1.e-08
        vcasingits = 8
        vcasingper = 1
        mcasingcal = 8

        odetol = 1.0e-07
        absreq = 1.0e-08
        relreq = 1.0e-08
        absacc = 1.0e-04
        epsr = 1.0e-08
        nPpts = 0
        Ppts = 0.0
        nPtrj(1:MNvol + 1) = -1
        LHevalues = .false.
        LHevectors = .false.
        LHmatrix = .false.
        Lperturbed = 0
        dpp = -1
        dqq = -1
        Lerrortype = 0
        Ngrid = -1
        dRZ = 1e-5
        Lcheck = 0
        Ltiming = .false.
        Ltransform = .false.
        fudge = 1.0e-00
        scaling = 1.0e-00

    end subroutine initialize_inputs

end module inputlist
