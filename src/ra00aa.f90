
subroutine ra00aa(writeorread)

    use constants, only: zero

    use numerical, only:

    use fileunits, only: ounit, aunit

    use inputlist, only: Wmacros, Wra00aa, Nfp, Mpol, Ntor, Lrad

    use cputiming, only: Tra00aa

    use allglobal, only: myid, ncpu, cpus, ext, Mvol, mn, im, in, Ate, Aze, Ato, Azo

    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    character, intent(in) :: writeorread

    logical :: exist

    integer :: vvol, oldMvol, oldMpol, oldNtor, oldmn, oldNfp, oldLrad, ii, jj, minLrad, llmodnp, ideriv, sumLrad
    integer, allocatable :: oldim(:), oldin(:)
    real(8), allocatable :: oldAte(:), oldAze(:), oldAto(:), oldAzo(:)
    real(8), allocatable :: allAte(:, :), allAze(:, :), allAto(:, :), allAzo(:, :)

    ideriv = 0

    select case (writeorread)

    case ('W')

        sumLrad = sum(Lrad(1:Mvol) + 1)

        allocate (allAte(1:sumLrad, 1:mn))
        allocate (allAze(1:sumLrad, 1:mn))
        allocate (allAto(1:sumLrad, 1:mn))
        allocate (allAzo(1:sumLrad, 1:mn))

        sumLrad = 1
        do vvol = 1, Mvol

            do ii = 1, mn
                allAte(sumLrad:sumLrad + Lrad(vvol), ii) = Ate(vvol, ideriv, ii)%s(0:Lrad(vvol))
                allAze(sumLrad:sumLrad + Lrad(vvol), ii) = Aze(vvol, ideriv, ii)%s(0:Lrad(vvol))
                allAto(sumLrad:sumLrad + Lrad(vvol), ii) = Ato(vvol, ideriv, ii)%s(0:Lrad(vvol))
                allAzo(sumLrad:sumLrad + Lrad(vvol), ii) = Azo(vvol, ideriv, ii)%s(0:Lrad(vvol))
            end do
            sumLrad = sumLrad + Lrad(vvol) + 1
        end do

        sumLrad = sum(Lrad(1:Mvol) + 1)

        deallocate (allAte)
        deallocate (allAze)
        deallocate (allAto)
        deallocate (allAzo)

    case ('R')

        if (myid == 0) then

            inquire (file="."//trim(ext)//".sp.A", exist=exist)

            if (.not. exist) then; write (ounit, '("ra00aa : ",f10.2," : myid=",i3," ; error ; .ext.sp.A does not exist ;")') cput - cpus, myid; goto 9998
            end if

            open (aunit, file="."//trim(ext)//".sp.A", status="old", form="unformatted", iostat=ios)

            if (ios /= 0) then; write (ounit, '("ra00aa : ",f10.2," : myid=",i3," ; error ; opening .ext.sp.A ;")') cput - cpus, myid; goto 9997
            end if

            read (aunit, iostat=ios) oldMvol, oldMpol, oldNtor, oldmn, oldNfp

            if (ios /= 0) then
                write (ounit, '("ra00aa : ",f10.2," : myid=",i3," ; error ; reading oldMvol, oldMpol, oldNtor, oldmn, oldNfp;")') cput - cpus, myid
                goto 9997
            end if

            if (oldNfp /= Nfp) then; write (ounit, '("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Nfp ; ")') cput - cpus, myid; goto 9997
            end if
            if (oldMvol /= Mvol) then; write (ounit, '("ra00aa : ",f10.2," : myid=",i3," ; error ; inconsistent Mvol ;")') cput - cpus, myid; goto 9997
            end if

            if (allocated(oldim)) deallocate (oldim)
            allocate (oldim(1:oldmn), stat=astat)
            oldim(1:oldmn) = 0
            if (allocated(oldin)) deallocate (oldin)
            allocate (oldin(1:oldmn), stat=astat)
            oldin(1:oldmn) = 0

            read (aunit, iostat=ios) oldim(1:oldmn)
            read (aunit, iostat=ios) oldin(1:oldmn)

            do vvol = 1, oldMvol

                read (aunit, iostat=ios) oldLrad

                minLrad = min(oldLrad, Lrad(vvol))

                if (allocated(oldAte)) deallocate (oldAte)
                allocate (oldAte(0:oldLrad), stat=astat)
                oldAte(0:oldLrad) = zero
                if (allocated(oldAze)) deallocate (oldAze)
                allocate (oldAze(0:oldLrad), stat=astat)
                oldAze(0:oldLrad) = zero
                if (allocated(oldAto)) deallocate (oldAto)
                allocate (oldAto(0:oldLrad), stat=astat)
                oldAto(0:oldLrad) = zero
                if (allocated(oldAzo)) deallocate (oldAzo)
                allocate (oldAzo(0:oldLrad), stat=astat)
                oldAzo(0:oldLrad) = zero

                do jj = 1, oldmn

                    read (aunit, iostat=ios) oldAte(0:oldLrad)
                    read (aunit, iostat=ios) oldAze(0:oldLrad)
                    read (aunit, iostat=ios) oldAto(0:oldLrad)
                    read (aunit, iostat=ios) oldAzo(0:oldLrad)

                    do ii = 1, mn
                        if (im(ii) == oldim(jj) .and. in(ii) == oldin(jj)) then; Ate(vvol, ideriv, ii)%s(0:minLrad) = oldAte(0:minLrad)
                            ; ; Aze(vvol, ideriv, ii)%s(0:minLrad) = oldAze(0:minLrad)
                            ; ; Ato(vvol, ideriv, ii)%s(0:minLrad) = oldAto(0:minLrad)
                            ; ; Azo(vvol, ideriv, ii)%s(0:minLrad) = oldAzo(0:minLrad)
                        end if
                    end do

                end do

                deallocate (oldAte, stat=astat)
                deallocate (oldAze, stat=astat)
                deallocate (oldAto, stat=astat)
                deallocate (oldAzo, stat=astat)

            end do

            deallocate (oldim, stat=astat)
            deallocate (oldin, stat=astat)

9997        continue

            close (aunit)

9998        continue

        end if

        do vvol = 1, Mvol

            llmodnp = 0

            do ii = 1, mn

            end do
            do ii = 1, mn

            end do

        end do

    case default

        if (.true.) then
            write (6, '("ra00aa :      fatal : myid=",i3," ; .true. ; invalid writeorread flag supplied on input;")') myid

            stop "ra00aa : .true. : invalid writeorread flag supplied on input ;"
        end if

    end select

end subroutine ra00aa

