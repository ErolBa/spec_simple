
subroutine packxi(NGdof, position, Mvol, mn, iRbc, iZbs, iRbs, iZbc, packorunpack, LComputeDerivatives, LComputeAxis)

    use constants, only: zero

    use numerical, only:

    use fileunits, only: ounit

    use inputlist, only: Wpackxi, Igeometry, Ntor, Nvol, Lfindzero

    use cputiming, only: Tpackxi

    use allglobal, only: ncpu, myid, cpus, im, in, MPI_COMM_SPEC, &
                         YESstellsym, NOTstellsym, &
                         ajk, Nt, Nz, Ntz, iRij, iZij, tRij, tZij, &
                         ijreal, ijimag, jireal, jiimag, efmn, ofmn, cfmn, sfmn, evmn, odmn, comn, simn, &
                         psifactor, Rscale

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    logical, intent(in) :: LComputeDerivatives
    logical, intent(in) :: LComputeAxis

    integer, intent(in) :: NGdof, Mvol, mn
    real(8) :: position(0:NGdof), iRbc(1:mn, 0:Mvol), iZbs(1:mn, 0:Mvol), iRbs(1:mn, 0:Mvol), iZbc(1:mn, 0:Mvol)
    character :: packorunpack

    integer :: lvol, jj, kk, irz, issym, idof, ifail, ivol

    idof = 0

    do lvol = 1, Mvol - 1

        do jj = 1, mn

            do irz = 0, 1

                if (Igeometry < 3 .and. irz == 1) cycle

                do issym = 0, 1

                    if (YESstellsym .and. issym == 1) cycle

                    if (issym == 0 .and. irz == 1 .and. jj == 1) cycle
                    if (issym == 1 .and. irz == 0 .and. jj == 1) cycle

                    idof = idof + 1

                    select case (packorunpack)

                    case ('P')

                        if (irz == 0 .and. issym == 0) position(idof) = iRbc(jj, lvol)/psifactor(jj, lvol)
                        if (irz == 1 .and. issym == 0) position(idof) = iZbs(jj, lvol)/psifactor(jj, lvol)
                        if (irz == 0 .and. issym == 1) position(idof) = iRbs(jj, lvol)/psifactor(jj, lvol)
                        if (irz == 1 .and. issym == 1) position(idof) = iZbc(jj, lvol)/psifactor(jj, lvol)

                    case ('U')

                        if (irz == 0 .and. issym == 0) iRbc(jj, lvol) = position(idof)*psifactor(jj, lvol)
                        if (irz == 1 .and. issym == 0) iZbs(jj, lvol) = position(idof)*psifactor(jj, lvol)
                        if (irz == 0 .and. issym == 1) iRbs(jj, lvol) = position(idof)*psifactor(jj, lvol)
                        if (irz == 1 .and. issym == 1) iZbc(jj, lvol) = position(idof)*psifactor(jj, lvol)

                    end select

                end do

            end do

        end do

    end do

    if (YESstellsym) then
        ; ; iZbs(1, 0:Mvol) = zero
        ; ; iRbs(1:mn, 0:Mvol) = zero
        ; ; iZbc(1:mn, 0:Mvol) = zero
    else
        ; ; iZbs(1, 0:Mvol) = zero
        ; ; iRbs(1, 0:Mvol) = zero
        ; 
    end if

    select case (packorunpack)

    case ('P')

    case ('U')

        ivol = 1

    end select

end subroutine packxi

