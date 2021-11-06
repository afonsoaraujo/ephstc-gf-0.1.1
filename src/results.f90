!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                     MODULE RESULTS                   ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!module results
   
subroutine Results (Rnflx,Lhflx,Shflx,Ghflx)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Analise results and generate statistics
!
!
!-----------------------------------------------------------------------      
   implicit none

   include 'globcst.inc'
   include 'stations.inc'
   include 'phycst.inc'
   include 'soilcst.inc'
   
   real      ::   Rnflx(NSTEPS)
   real      ::   Lhflx(NSTEPS)
   real      ::   Shflx(NSTEPS)
   real      ::   Ghflx(NSTEPS)    

   integer   ::   NVARS,NCOLS
   parameter(NVARS=4, NCOLS=27)
   
   real      ::   VARCal(NVARS,NSTEPS)
   real      ::   VARObs(NVARS,NSTEPS)
   real      ::   buffer1(NCOLS)
   integer   ::   it,ncol,icol
   real      ::   E,RMSE,MBE,CORR
   integer   ::   aux,i,j,k,l
   
   real,allocatable      :: EBObs(:)
!   real,allocatable      :: EBCal(:)
   integer   ::   istatus
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Início do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   
   aux  = len_trim(flx1doutfile)
   do k = 1, nflxstations
      write(flx1doutfile(aux+1:aux+12),'(i3)') flxstation(k)
      open(20+k,file=flx1doutfile,form='formatted',status='unknown')
   enddo

! --------------------------------------------------------------------   
!  Allocate memory - Aloca memória
! --------------------------------------------------------------------
   if (.not. allocated(EBObs)) ALLOCATE(EBObs (nsteps),STAT=istatus)
   CALL check_alloc_status(istatus, "ephst:EBObs")
   EBObs = 0
   
!   if (.not. allocated(EBCal)) ALLOCATE(EBCal (nx,ny),STAT=istatus)
!   CALL check_alloc_status(istatus, "ephst:EBCal")
!   EBCal = 0   

! Verifica balanco de energia observado
! --------------------------------------------------------------------
   do it = 1, NSTEPS
      EBObs(it) = Rnflx(it) - Lhflx(it) - Shflx(it) - Ghflx(it)
      write(60,'(5f18.8)') EBObs(it),Rnflx(it),Lhflx(it),Shflx(it),Ghflx(it)
   enddo


!  Constroi o vetor com as variaveis observadas
!-----------------------------------------------------------------------
   do it = 1, NSTEPS
      VARObs(1,it) = Rnflx(it)
      VARObs(2,it) = Lhflx(it)
      VARObs(3,it) = Shflx(it)
      VARObs(4,it) = Ghflx(it)
   enddo

   write(140,'(4f18.8)') emissg,emissa,albdveg(10)
   do k = 1, nflxstations
      i=flxpoints(k,1)
      j=flxpoints(k,2)

!  Constroi o vetor temp com as variaveis a serem analisadas
!-----------------------------------------------------------------------
      do it = 1, NSTEPS
         icol = 7              ! coluna do arquivo de saida/leitura onde esta avariavel...

!      write(6,'(a,i4,3f18.8)') 'Estou aqui!',icol,VARObs(1,it),VARObs(2,it),VARObs(3,it)
! Le variaveis calculadas do arquivo de saida do EPHST
!-----------------------------------------------------------------------
         read(20+k,*) (buffer1(ncol), ncol=1,NCOLS)
         do l = 1, NVARS
            VARCal(l,it) = buffer1(icol)
	    icol = icol + 1
         enddo
         write(110+k,'(8f18.8)') VARCal(1,it),VARObs(1,it),VARCal(2,it),VARObs(2,it), &
	                         VARCal(3,it),VARObs(3,it),VARCal(4,it),VARObs(4,it)
      enddo
   
!      write(6,'(a)')  '   '
      write(140,'(a,2i4)') 'ponto ==>', i,j
      do l = 1, NVARS
         call Results_Nash (VARCal(l,:),VARObs(l,:),E,RMSE,MBE,CORR)
!        write(6,'(3f18.8)') E,RMSE,MBE
      enddo   ! enddo NSTEPS

   enddo   ! enddo nflxstations
   
   do k = 1, nflxstations
      close (20+k)
      close (110+k)
   enddo

return

end subroutine Results


subroutine Results_Nash (VARCal,VARObs,E,RMSE,MBE,CORR)

   implicit none
   include 'globcst.inc'   

   real      ::   VARObs(NSTEPS)
   real      ::   VARCal(NSTEPS)

   real      ::   SUMVARObs,SUMVARCal
   real      ::   SUMVARObs2,SUMVARCal2,SUMVARObsCal
   real      ::   SUMN,SUMD,SUM1,SUM2
   real      ::   VARObsM,VARCalM
   real      ::   E,RMSE,MBE,CORR

   real      ::   VARQ,VARE
   integer   ::   IT

!
!  OBJECTIVE FUNCTION CALCULATIONS
!-----------------------------------------------------------------------
   SUMVARObs = 0.
   SUMVARCal = 0.
   SUMVARCal2 = 0.
   SUMVARObs2 = 0.
   SUMVARObsCAl = 0.
   SUMN  = 0.
   SUMD  = 0.
   SUM1  = 0.
   SUM2  = 0.

   do IT=1,NSTEPS
      SUMVARObs = SUMVARObs + VARObs(IT)
      SUMVARCal = SUMVARCal + VARCal(IT) 
      SUMVARObs2 = SUMVARObs2 +  VARObs(IT)*VARObs(IT)
      SUMVARCal2 = SUMVARCal2 +  VARCal(IT)*VARCal(IT)
      SUMVARObsCal = SUMVARObsCal + VARObs(IT)*VARCal(IT)
   enddo
   VARObsM = SUMVARObs/NSTEPS
   VARCalM = SUMVArCal/NSTEPS

   do IT=1,NSTEPS
      SUMN = SUMN + (VARObs(IT)-VARCal(IT)) * (VARObs(IT)-VARCal(IT))
      SUMD = SUMD + (VARObs(IT)-VARObsM) * (VARObs(IT)-VARObsM)
      SUM1 = SUM1 + (VARCal(IT)-VARObs(IT)) * (VARCal(IT)-VARObs(IT))
      SUM2 = SUM2 + (VARCal(IT)-VARObs(IT))
   enddo
      
   E    = 1-(SUMN/SUMD)
   RMSE = sqrt(SUM1/NSTEPS)
   MBE  = SUM2/NSTEPS
   CORR = (SUMVARObsCal - (1/NSTEPS)*(SUMVARCal * SUMVARObs)) /       &
          sqrt( (SUMVARCal2 - (1/NSTEPS)*(SUMVARCal*SUMVARCal)) *     &
	        (SUMVARObs2 - (1/NSTEPS)*(SUMVARObs*SUMVARObs)))


   write(140,'(5f18.8)') E,RMSE,MBE,VARObsM,CORR


return

end subroutine Results_Nash
