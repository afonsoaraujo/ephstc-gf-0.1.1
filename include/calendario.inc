! -------------------------------------------------------------------------
!  Calendar: Trata objetos de um calendario
!
!   Nelson Luis da Costa Dias
!   21-abr-1988  (em PASCAL)
!   19-jun-1991  (em C)
!
!   Afonso Augusto Araujo
!   22-fev-2002  (em Fortran)
! ------------------------------------------------------------------------- 

!!DEC$ IF .not. DEFINED (CALENDARIO_INC)
!!DEC$ DEFINE CALENDARIO_INC

      integer        ::   NNDias(2,12)
      integer        ::   NDiasAno(2)
      character(4)   ::   NomMes(2,12)

      integer        ::   bb
      integer        ::   mm
      
      data ((NNDias(bb,mm),mm=1,12),bb=1,2)                           &
                  / 31,28,31,30,31,30,31,31,30,31,30,31,              & 
                    31,29,31,30,31,30,31,31,30,31,30,31 /


      data ((NomMes(bb,mm),mm=1,12),bb=1,2)                           &
                  / 'JAN','FEV','MAR','ABR','MAI','JUN',              &
                    'JUL','AGO','SET','OUT','NOV','DEZ',              &
                    'JAN','FEB','MAR','APR','MAY','JUN',              &
                    'JUL','AUG','SEP','OCT','NOV','DEC' /

      
    
      data NDiasAno / 365,366 /

!!DEC$ ENDIF
