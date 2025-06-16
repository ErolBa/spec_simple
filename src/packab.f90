
subroutine packab( packorunpack, lvol, NN, solution, ideriv )



  use constants, only : zero

  use numerical, only : small

  use fileunits, only : ounit

  use inputlist, only : Wpackab, Lrad

  use cputiming, only : Tpackab

  use allglobal, only : myid, ncpu, cpus, MPI_COMM_SPEC, &
                        mn, im, in, Ate, Aze, Ato, Azo, YESstellsym, NOTstellsym, &
                        TT, YESMatrixFree, &
                        Lma, Lmb, Lmc, Lmd, Lme, Lmf, Lmg, Lmh, &
                        Lmavalue, Lmbvalue, Lmcvalue, Lmdvalue, Lmevalue, Lmfvalue, Lmgvalue, Lmhvalue



  LOCALS

  CHARACTER, intent(in) :: packorunpack
  INTEGER  , intent(in) :: lvol, NN, ideriv
  REAL                  :: solution(1:NN)

  INTEGER               :: ii, ll, id, llrad

  BEGIN(packab)



  llrad = Lrad(lvol) ! shorthand;



  select case( packorunpack )

  case( 'U' )

   if( YESstellsym ) then

    do ii = 1, mn
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ;                           ; Ato(lvol,ideriv,ii)%s(ll) = zero
      ;               ;                           ; Azo(lvol,ideriv,ii)%s(ll) = zero
     enddo ! end of do ll;
    enddo ! end of do ii;

    if (YESMatrixFree) then
      do ii = 1, mn
        ;                  ; id = Lma(lvol,ii)
        if (id/=0) then; Lmavalue(lvol,ii) = solution(id)
        else           ; Lmavalue(lvol,ii) = zero
        endif

        ;                  ; id = Lmb(lvol,ii)
        if (id/=0) then; Lmbvalue(lvol,ii) = solution(id)
        else           ; Lmbvalue(lvol,ii) = zero
        endif

        if( ii.gt.1 ) then
          ;                ; id = Lme(lvol,ii)
          if (id/=0) then; Lmevalue(lvol,ii) = solution(id)
          else           ; Lmevalue(lvol,ii) = zero
          endif

        else               ; id = Lmg(lvol,ii)
          if (id/=0) then; Lmgvalue(lvol,ii) = solution(id)
          else           ; Lmgvalue(lvol,ii) = zero
          endif

          ;                ; id = Lmh(lvol,ii)
          if (id/=0) then; Lmhvalue(lvol,ii) = solution(id)
          else           ; Lmhvalue(lvol,ii) = zero
          endif
        endif
      enddo ! ii
    endif ! YESMatrixFree
   else ! NOTstellsym;

    ;  ii = 1
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ;                           ; Ato(lvol,ideriv,ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
      ;               ;                           ; Azo(lvol,ideriv,ii)%s(ll) = zero         ! sin( m \t - n \z ) = 0 for (m,n)=(0,0);
     enddo
    do ii = 2, mn
     do ll = 0, llrad ; id = Ate(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ate(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ate(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Aze(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Aze(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Aze(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Ato(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Ato(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Ato(lvol,ideriv,ii)%s(ll) = zero
       endif
      ;               ; id = Azo(lvol,0,ii)%i(ll) ;
       if (id/=0) then; Azo(lvol,ideriv,ii)%s(ll) = solution(id)
       else           ; Azo(lvol,ideriv,ii)%s(ll) = zero
       endif
     enddo ! end of do ll;
    enddo ! end of do ii;

    if (YESMatrixFree) then
      do ii = 1, mn
        ;                  ; id = Lma(lvol,ii)
        if (id/=0) then; Lmavalue(lvol,ii) = solution(id)
        else           ; Lmavalue(lvol,ii) = zero
        endif

        ;                  ; id = Lmb(lvol,ii)
        if (id/=0) then; Lmbvalue(lvol,ii) = solution(id)
        else           ; Lmbvalue(lvol,ii) = zero
        endif

        if( ii.gt.1 ) then ; id = Lmc(lvol,ii)
          if (id/=0) then; Lmcvalue(lvol,ii) = solution(id)
          else           ; Lmcvalue(lvol,ii) = zero
          endif

          ;                ; id = Lmd(lvol,ii)
          if (id/=0) then; Lmdvalue(lvol,ii) = solution(id)
          else           ; Lmdvalue(lvol,ii) = zero
          endif

          ;                ; id = Lme(lvol,ii)
          if (id/=0) then; Lmevalue(lvol,ii) = solution(id)
          else           ; Lmevalue(lvol,ii) = zero
          endif

          ;                ; id = Lmf(lvol,ii)
          if (id/=0) then; Lmfvalue(lvol,ii) = solution(id)
          else           ; Lmfvalue(lvol,ii) = zero
          endif

        else               ; id = Lmg(lvol,ii)
          if (id/=0) then; Lmgvalue(lvol,ii) = solution(id)
          else           ; Lmgvalue(lvol,ii) = zero
          endif

          ;                ; id = Lmh(lvol,ii)
          if (id/=0) then; Lmhvalue(lvol,ii) = solution(id)
          else           ; Lmhvalue(lvol,ii) = zero
          endif
        endif
      enddo ! ii
    endif ! YESMatrixFree
   endif ! end of if( YESstellsym );

  case( 'P' )


   solution = zero

   do ii = 1, mn
    do ll = 0, llrad
     ;                                    ; id = Ate(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Ate(lvol,ideriv,ii)%s(ll)
     ;                                    ; id = Aze(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Aze(lvol,ideriv,ii)%s(ll)
     if( ii.gt.1 .and. NOTstellsym ) then ; id = Ato(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Ato(lvol,ideriv,ii)%s(ll)
      ;                                   ; id = Azo(lvol,0,ii)%i(ll) ; if (id/=0) solution(id) = Azo(lvol,ideriv,ii)%s(ll)
     endif
    enddo ! end of do ll;
   enddo ! end of do ii;

  end select ! end of select case( packorunpack );


  RETURN(packab)



end subroutine packab


