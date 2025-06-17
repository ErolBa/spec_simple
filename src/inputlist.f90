
module inputlist

    implicit none

    integer, parameter :: MNvol = 256
    integer, parameter :: MMpol = 128
    integer, parameter :: MNtor = 128

    integer :: Igeometry = 3
    integer :: Istellsym = 1
    integer :: Lfreebound = 0
    real(8) :: phiedge = 1.0
    real(8) :: curtor = 0.0
    real(8) :: curpol = 0.0
    real(8) :: gamma = 0.0
    integer :: Nfp = 1
    integer :: Nvol = 1
    integer :: Mpol = 0
    integer :: Ntor = 0
    integer :: Lrad(1:MNvol + 1) = 4
    integer :: Lconstraint = -1
    real(8) :: tflux(1:MNvol + 1) = 0.0
    real(8) :: pflux(1:MNvol + 1) = 0.0
    real(8) :: helicity(1:MNvol) = 0.0
    real(8) :: pscale = 0.0
    real(8) :: pressure(1:MNvol + 1) = 0.0
    integer :: Ladiabatic = 0
    real(8) :: adiabatic(1:MNvol + 1) = 0.0
    real(8) :: mu(1:MNvol + 1) = 0.0
    real(8) :: Ivolume(1:MNvol + 1) = 0.0
    real(8) :: Isurf(1:MNvol) = 0.0
    integer :: pl(0:MNvol) = 0
    integer :: ql(0:MNvol) = 0
    integer :: pr(0:MNvol) = 0
    integer :: qr(0:MNvol) = 0
    real(8) :: iota(0:MNvol) = 0.0
    integer :: lp(0:MNvol) = 0
    integer :: lq(0:MNvol) = 0
    integer :: rp(0:MNvol) = 0
    integer :: rq(0:MNvol) = 0
    real(8) :: oita(0:MNvol) = 0.0
    real(8) :: mupftol = 1.0e-14
    integer :: mupfits = 8
    real(8) :: rpol = 1.0
    real(8) :: rtor = 1.0
    integer :: Lreflect = 0

    real(8) :: Rac(0:MNtor) = 0.0
    real(8) :: Zas(0:MNtor) = 0.0
    real(8) :: Ras(0:MNtor) = 0.0
    real(8) :: Zac(0:MNtor) = 0.0

    real(8) :: Rbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Zbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Rbs(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Zbc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

    real(8) :: Rwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Zws(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Rws(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Zwc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

    real(8) :: Vns(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Bns(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Vnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0
    real(8) :: Bnc(-MNtor:MNtor, -MMpol:MMpol) = 0.0

    integer :: Linitialize = 0
    integer :: LautoinitBn = 1
    integer :: Lzerovac = 0
    integer :: Ndiscrete = 2
    integer :: Nquad = -1
    integer :: iMpol = -4
    integer :: iNtor = -4
    integer :: Lsparse = 0
    integer :: Lsvdiota = 0
    integer :: imethod = 3
    integer :: iorder = 2
    integer :: iprecon = 0
    real(8) :: iotatol = -1.0
    integer :: Lextrap = 0
    integer :: Mregular = -1
    integer :: Lrzaxis = 1
    integer :: Ntoraxis = 3

    integer :: LBeltrami = 4
    integer :: Linitgues = 1
    integer :: Lposdef = 0
    real(8) :: maxrndgues = 1.0
    integer :: Lmatsolver = 3
    integer :: NiterGMRES = 200
    real(8) :: epsGMRES = 1e-14
    integer :: LGMRESprec = 1
    real(8) :: epsILU = 1e-12

    integer :: Lfindzero = 0
    real(8) :: escale = 0.0
    real(8) :: opsilon = 1.0
    real(8) :: pcondense = 2.0
    real(8) :: epsilon = 0.0
    real(8) :: wpoloidal = 1.0
    real(8) :: upsilon = 1.0
    real(8) :: forcetol = 1.0e-10
    real(8) :: c05xmax = 1.0e-06
    real(8) :: c05xtol = 1.0e-12
    real(8) :: c05factor = 1.0e-02
    logical :: LreadGF = .true.
    integer :: mfreeits = 0
    real(8) :: bnstol = 1.0e-06
    real(8) :: bnsblend = 0.666
    real(8) :: gBntol = 1.0e-06
    real(8) :: gBnbld = 0.666
    real(8) :: vcasingeps = 1.e-12
    real(8) :: vcasingtol = 1.e-08
    integer :: vcasingits = 8
    integer :: vcasingper = 1
    integer :: mcasingcal = 8

    real(8) :: odetol = 1.0e-07
    real(8) :: absreq = 1.0e-08
    real(8) :: relreq = 1.0e-08
    real(8) :: absacc = 1.0e-04
    real(8) :: epsr = 1.0e-08
    integer :: nPpts = 0
    real(8) :: Ppts = 0.0
    integer :: nPtrj(1:MNvol + 1) = -1
    logical :: LHevalues = .false.
    logical :: LHevectors = .false.
    logical :: LHmatrix = .false.
    integer :: Lperturbed = 0
    integer :: dpp = -1
    integer :: dqq = -1
    integer :: Lerrortype = 0
    integer :: Ngrid = -1
    real(8) :: dRZ = 1e-5
    integer :: Lcheck = 0
    logical :: Ltiming = .false.
    logical :: Ltransform = .false.
    real(8) :: fudge = 1.0e-00
    real(8) :: scaling = 1.0e-00

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
    logical :: Wbuild_vector_potential = .false.
    logical :: Wreadin = .false.
    logical :: Wwrtend = .false.
    logical :: Wmacros = .false.

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
