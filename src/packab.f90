
subroutine packab(packorunpack, lvol, NN, solution, ideriv)

    use constants, only: zero

    use numerical, only: small

    use fileunits, only: ounit

    use inputlist, only: Wpackab, Lrad

    use cputiming, only: Tpackab

    use allglobal, only: myid, ncpu, cpus, MPI_COMM_SPEC, &
                         mn, im, in, Ate, Aze, Ato, Azo, YESstellsym, NOTstellsym, &
                         TT, &
                         Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                         Lmavalue, Lmbvalue, Lmcvalue, Lmdvalue, Lmevalue, Lmfvalue, Lmgvalue, Lmhvalue

    use mpi
    implicit none
    integer :: ierr, astat, ios, nthreads, ithread
    real(8) :: cput, cpui, cpuo = 0

    character, intent(in) :: packorunpack
    integer, intent(in) :: lvol, NN, ideriv
    real(8) :: solution(1:NN)

    integer :: ii, ll, id, llrad

    llrad = Lrad(lvol) ! shorthand;

    select case (packorunpack)

    case ('U')

        if (YESstellsym) then

            do ii = 1, mn
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ate(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Aze(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; ; Ato(lvol, ideriv, ii)%s(ll) = zero
                    ; ; ; Azo(lvol, ideriv, ii)%s(ll) = zero
                end do ! end of do ll;
            end do ! end of do ii;

        else ! NOTstellsym;

            ; ii = 1
            do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                else; Ate(lvol, ideriv, ii)%s(ll) = zero
                end if
                ; ; id = Aze(lvol, 0, ii)%i(ll); 
                if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                else; Aze(lvol, ideriv, ii)%s(ll) = zero
                end if
                ; ; ; Ato(lvol, ideriv, ii)%s(ll) = zero ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
                ; ; ; Azo(lvol, ideriv, ii)%s(ll) = zero ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
            end do
            do ii = 2, mn
                do ll = 0, llrad; id = Ate(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ate(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ate(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Aze(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Aze(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Aze(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Ato(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Ato(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Ato(lvol, ideriv, ii)%s(ll) = zero
                    end if
                    ; ; id = Azo(lvol, 0, ii)%i(ll); 
                    if (id /= 0) then; Azo(lvol, ideriv, ii)%s(ll) = solution(id)
                    else; Azo(lvol, ideriv, ii)%s(ll) = zero
                    end if
                end do ! end of do ll;
            end do ! end of do ii;

        end if ! end of if( YESstellsym );

    case ('P')

        solution = zero

        do ii = 1, mn
            do ll = 0, llrad
                ; ; id = Ate(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Ate(lvol, ideriv, ii)%s(ll)
                ; ; id = Aze(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Aze(lvol, ideriv, ii)%s(ll)
                if (ii > 1 .and. NOTstellsym) then; id = Ato(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Ato(lvol, ideriv, ii)%s(ll)
                    ; ; id = Azo(lvol, 0, ii)%i(ll); if (id /= 0) solution(id) = Azo(lvol, ideriv, ii)%s(ll)
                end if
            end do ! end of do ll;
        end do ! end of do ii;

    end select ! end of select case( packorunpack );

end subroutine packab

