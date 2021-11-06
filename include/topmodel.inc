
! constantes
!---------------------------------------------------------------------
   integer    ::   NACMX
   parameter(NACMX=40)

   integer    ::   NCHMX
   parameter(NCHMX=10)
   
! parametros  (Obs.: colocar no NAMELIST)
!---------------------------------------------------------------------   
   real      ::   CHVDT
   real      ::   RVDT
   real      ::   Q0

   common /topmodel001/ CHVDT,RVDT,Q0
   
   character(128)   ::   itfile    ! arquivo contendo IT ordenado
   
   common /topmodel002/ itfile
   
   integer   ::   IN
   
   common /topmodel003/ IN
   
