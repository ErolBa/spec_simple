
subroutine allocate_Beltrami_matrices(vvol, LcomputeDerivatives)

    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer, intent(in) :: vvol
    logical, intent(in) :: LcomputeDerivatives
    integer :: NN

    NN = NAdof(vvol)

    if (allocated(dMA)) deallocate (dMA)
    allocate (dMA(0:NN, 0:NN), stat=astat)
    dMA(0:NN, 0:NN) = zero
    if (allocated(dMD)) deallocate (dMD)
    allocate (dMD(0:NN, 0:NN), stat=astat)
    dMD(0:NN, 0:NN) = zero

    if (allocated(dMB)) deallocate (dMB)
    allocate (dMB(0:NN, 0:2), stat=astat)
    dMB(0:NN, 0:2) = zero
    if (allocated(dMG)) deallocate (dMG)
    allocate (dMG(0:NN), stat=astat)
    dMG(0:NN) = zero

    if (allocated(solution)) deallocate (solution)
    allocate (solution(1:NN, -1:2), stat=astat)
    solution(1:NN, -1:2) = zero

    if (allocated(MBpsi)) deallocate (MBpsi)
    allocate (MBpsi(1:NN), stat=astat)
    MBpsi(1:NN) = zero

end subroutine allocate_Beltrami_matrices

subroutine deallocate_Beltrami_matrices(LcomputeDerivatives)

    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LcomputeDerivatives

    deallocate (dMA, stat=astat)
    deallocate (dMD, stat=astat)

    deallocate (dMB, stat=astat)

    deallocate (dMG, stat=astat)

    deallocate (solution, stat=astat)

    deallocate (MBpsi, stat=astat)

end subroutine deallocate_Beltrami_matrices

subroutine allocate_geometry_matrices(vvol, LcomputeDerivatives)

    use constants, only: zero

    use fileunits

    use inputlist, only: Wmemory, Wmacros, Mpol, Lrad

    use allglobal

    use cputiming

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    integer :: vvol

    logical, intent(in) :: LcomputeDerivatives

    integer :: ll, lldof, jjdof, iidof

    ll = Lrad(vvol)

    if (Lcoordinatesingularity) then
        lldof = (Lrad(vvol) - mod(Lrad(vvol), 2))/2

        iidof = mn
        jjdof = mn

    else
        lldof = Lrad(vvol)

        iidof = mn
        jjdof = mn

    end if

    if (allocated(guvijsave)) deallocate (guvijsave)
    allocate (guvijsave(1:Ntz, 1:3, 1:3, 1:Iquad(vvol)), stat=astat)
    guvijsave(1:Ntz, 1:3, 1:3, 1:Iquad(vvol)) = zero

    if (allocated(DToocc)) deallocate (DToocc)
    allocate (DToocc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DToocc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(TTssss)) deallocate (TTssss)
    allocate (TTssss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TTssss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(TDstsc)) deallocate (TDstsc)
    allocate (TDstsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TDstsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(TDszsc)) deallocate (TDszsc)
    allocate (TDszsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    TDszsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(DDttcc)) deallocate (DDttcc)
    allocate (DDttcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDttcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(DDtzcc)) deallocate (DDtzcc)
    allocate (DDtzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDtzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
    if (allocated(DDzzcc)) deallocate (DDzzcc)
    allocate (DDzzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
    DDzzcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

    if (allocated(Tss)) deallocate (Tss)
    allocate (Tss(0:lldof, 1:mn), stat=astat)
    Tss(0:lldof, 1:mn) = zero
    if (allocated(Dtc)) deallocate (Dtc)
    allocate (Dtc(0:lldof, 1:mn), stat=astat)
    Dtc(0:lldof, 1:mn) = zero
    if (allocated(Dzc)) deallocate (Dzc)
    allocate (Dzc(0:lldof, 1:mn), stat=astat)
    Dzc(0:lldof, 1:mn) = zero
    if (allocated(Ttc)) deallocate (Ttc)
    allocate (Ttc(0:lldof, 1:mn), stat=astat)
    Ttc(0:lldof, 1:mn) = zero
    if (allocated(Tzc)) deallocate (Tzc)
    allocate (Tzc(0:lldof, 1:mn), stat=astat)
    Tzc(0:lldof, 1:mn) = zero

    if (NOTstellsym) then

        if (allocated(DToocs)) deallocate (DToocs)
        allocate (DToocs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DToocs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DToosc)) deallocate (DToosc)
        allocate (DToosc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DToosc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DTooss)) deallocate (DTooss)
        allocate (DTooss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DTooss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(TTsscc)) deallocate (TTsscc)
        allocate (TTsscc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsscc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TTsscs)) deallocate (TTsscs)
        allocate (TTsscs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsscs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TTsssc)) deallocate (TTsssc)
        allocate (TTsssc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TTsssc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(TDstcc)) deallocate (TDstcc)
        allocate (TDstcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TDstcs)) deallocate (TDstcs)
        allocate (TDstcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TDstss)) deallocate (TDstss)
        allocate (TDstss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDstss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(TDszcc)) deallocate (TDszcc)
        allocate (TDszcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszcc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TDszcs)) deallocate (TDszcs)
        allocate (TDszcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(TDszss)) deallocate (TDszss)
        allocate (TDszss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        TDszss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(DDttcs)) deallocate (DDttcs)
        allocate (DDttcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDttsc)) deallocate (DDttsc)
        allocate (DDttsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDttss)) deallocate (DDttss)
        allocate (DDttss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDttss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(DDtzcs)) deallocate (DDtzcs)
        allocate (DDtzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDtzsc)) deallocate (DDtzsc)
        allocate (DDtzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDtzss)) deallocate (DDtzss)
        allocate (DDtzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDtzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(DDzzcs)) deallocate (DDzzcs)
        allocate (DDzzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzcs(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDzzsc)) deallocate (DDzzsc)
        allocate (DDzzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzsc(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero
        if (allocated(DDzzss)) deallocate (DDzzss)
        allocate (DDzzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof), stat=astat)
        DDzzss(0:lldof, 0:lldof, 1:iidof, 1:jjdof) = zero

        if (allocated(Tsc)) deallocate (Tsc)
        allocate (Tsc(0:lldof, 1:mn), stat=astat)
        Tsc(0:lldof, 1:mn) = zero
        if (allocated(Dts)) deallocate (Dts)
        allocate (Dts(0:lldof, 1:mn), stat=astat)
        Dts(0:lldof, 1:mn) = zero
        if (allocated(Dzs)) deallocate (Dzs)
        allocate (Dzs(0:lldof, 1:mn), stat=astat)
        Dzs(0:lldof, 1:mn) = zero
        if (allocated(Tts)) deallocate (Tts)
        allocate (Tts(0:lldof, 1:mn), stat=astat)
        Tts(0:lldof, 1:mn) = zero
        if (allocated(Tzs)) deallocate (Tzs)
        allocate (Tzs(0:lldof, 1:mn), stat=astat)
        Tzs(0:lldof, 1:mn) = zero

    end if

end subroutine allocate_geometry_matrices

subroutine deallocate_geometry_matrices(LcomputeDerivatives)

    use constants, only: zero

    use fileunits

    use inputlist, only: Wmemory, Wmacros

    use allglobal

    use cputiming

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LcomputeDerivatives

    Lsavedguvij = .false.
    deallocate (guvijsave, stat=astat)

    deallocate (DToocc, stat=astat)
    deallocate (TTssss, stat=astat)
    deallocate (TDstsc, stat=astat)
    deallocate (TDszsc, stat=astat)
    deallocate (DDttcc, stat=astat)
    deallocate (DDtzcc, stat=astat)
    deallocate (DDzzcc, stat=astat)

    deallocate (Tss, stat=astat)
    deallocate (Dtc, stat=astat)
    deallocate (Dzc, stat=astat)
    deallocate (Ttc, stat=astat)
    deallocate (Tzc, stat=astat)

    if (NOTstellsym) then

        deallocate (DToocs, stat=astat)
        deallocate (DToosc, stat=astat)
        deallocate (DTooss, stat=astat)

        deallocate (TTsscc, stat=astat)
        deallocate (TTsscs, stat=astat)
        deallocate (TTsssc, stat=astat)

        deallocate (TDstcc, stat=astat)
        deallocate (TDstcs, stat=astat)
        deallocate (TDstss, stat=astat)

        deallocate (TDszcc, stat=astat)
        deallocate (TDszcs, stat=astat)
        deallocate (TDszss, stat=astat)

        deallocate (DDttcs, stat=astat)
        deallocate (DDttsc, stat=astat)
        deallocate (DDttss, stat=astat)

        deallocate (DDtzcs, stat=astat)
        deallocate (DDtzsc, stat=astat)
        deallocate (DDtzss, stat=astat)

        deallocate (DDzzcs, stat=astat)
        deallocate (DDzzsc, stat=astat)
        deallocate (DDzzss, stat=astat)

        deallocate (Tsc, stat=astat)
        deallocate (Dts, stat=astat)
        deallocate (Dzs, stat=astat)
        deallocate (Tts, stat=astat)
        deallocate (Tzs, stat=astat)

    end if

end subroutine deallocate_geometry_matrices
