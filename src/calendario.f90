!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                   MODULE CALENDARIO                  ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module Calendario
! -------------------------------------------------------------------------
!   Calendar: Trata objetos de um calendario
!
!   Nelson Luis da Costa Dias
!   21-abr-1988  (em PASCAL)
!   19-jun-1991  (em C)
!   Afonso Augusto Magalhães de Araujo
!   22-fev-2002  (compilado em FORTRAN90)
!   Afonso Augusto Magalhãaes de Araujo
!   15-out-2004 (módulo em Fortran90)   
! -------------------------------------------------------------------------

!declaração de variáveis 
!   include 'calendario.inc'
   
!contains

!
! -------------------------------------------------------------------------
!   --> Bissexto: Diz se um ano e' bissexto ou nao
!
!          Entrada:
!             Ano      -- Numero do ano de 0 em diante
!
!          Saida:
!             Bissexto -- Um valor logico ( TRUE ou FALSE )
! -------------------------------------------------------------------------
!
   SUBROUTINE BISSEXTO( Ano, Bis )

   implicit none
   
   include 'calendario.inc'

   integer   ::   Ano
   integer   ::   Bis

   IF ( (MOD(Ano,4)   .eq. 0) .and. (MOD(Ano,100) .ne. 0) .or.        &
      (MOD(Ano,400) .eq. 0) ) THEN
       Bis = 2
   ELSE
       Bis = 1
   ENDIF

   RETURN   

   end subroutine BISSEXTO

!
!-----------------------------------------------------------------------
!  --> DataCor: Dado o fato de o ano ser ou nao bissexto e um par mes e
!               dia, retorna um dia corrido entre 1 e 366
!
!         Entrada:
!            Bis    -- Se o ano e' ou nao bissexto
!            Mes    -- Mes do ano
!            Dia    -- Dia do mes
!
!         Saida:
!            DiAno  -- Dia do ano
!------------------------------------------------------------------------
   SUBROUTINE DATACOR( Bis, Mes, Dia, DiAno)

   implicit none

   integer   ::   Bis
   integer   ::   Mes
   integer   ::   Dia
   integer   ::   DiAno

   integer    ::  k
   
   include 'calendario.inc'

       
   DiAno = Dia
   k = 1 ;
!      write (6, '(3x,a,i4,a)') 'k = ', k,','
!      write (6, '(3x,a,i4,a)') 'Mes = ', Mes,','
   do while  ( k .lt. Mes )
      DiAno = DiAno + NNDias(Bis, k)
      k = k + 1
   end do

   RETURN

  END subroutine DATACOR


!end module Calendario

