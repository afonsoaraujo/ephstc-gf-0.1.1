!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                   MODULE AGGREGATE                   ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module aggregate


!/*****************************************************************************
!  Aggregate()
!  
!  Calculate the average values for the different fluxes and state variables
!  over the basin.  Only the runoff is calculated as a total volume instead
!  of an average.  In the current implementation the local radiation
!  elements are not stored for the entire area.  Therefore these components
!  are aggregated in AggregateRadiation() inside MassEnergyBalance().
!
!  
!*****************************************************************************/
subroutine Aggregate(nx,ny,maskmap,rnflx,shflx,lhflx,ghflx,precip,    &
                     RunOff,SatFlow,OutFlow,wetsfc,wetdp2,wetdp3,     &
	             wetcnp,Wini,TTRnflxAvg,TTLhflxAvg,TTShflxAvg,    &
	             TTGhflxAvg,TTPrecipAvg,TTEvapAvg,TTRunOffAvg,    &
		     TTSatFlowAvg,TTOutflowAvg,TTWetAvg,TTWiniAvg,    &
		     TTRnflxVar,TTLhflxVar,TTShflxVar,TTGhflxVar,     &
		     TTPrecipVar,TTEvapVar,TTRunOffVar,TTSatFlowVar,  &
		     TTOutflowVar,TTWetVar,TTWiniVar,TArea)

   implicit none
   include 'globcst.inc'
   include 'phycst.inc'
   include 'soilcst.inc'
   include 'grid.inc'
   include 'alloc.inc'
   
   integer   ::   nx,ny
   integer   ::   maskmap   (nx,ny)
   real(4)   ::   rnflx     (nx,ny)
   real(4)   ::   shflx     (nx,ny)
   real(4)   ::   lhflx     (nx,ny)
   real(4)   ::   ghflx     (nx,ny)
   real(4)   ::   precip    (nx,ny)     ! Precipitation rate at the surface   
   real(4)   ::   RunOff    (nx,ny)     ! Total RunOff 
   real(4)   ::   SatFlow   (nx,ny)
   real(4)   ::   OutFlow   (nx,ny)
   real(4)   ::   wetsfc    (nx,ny)   
   real(4)   ::   wetdp2    (nx,ny)   
   real(4)   ::   wetdp3    (nx,ny)   
   real(4)   ::   wetcnp    (nx,ny)
   real(4)   ::   Wini      (nx,ny)
   real      ::   TTRnflxAvg,  TTRnflxVar
   real      ::   TTShflxAvg,  TTShflxVar
   real      ::   TTLhflxAvg,  TTLhflxVar
   real      ::   TTGhflxAvg,  TTGhflxVar
   real      ::   TTPrecipAvg, TTPrecipVar
   real      ::   TTEvapAvg,   TTEvapVar
   real      ::   TTRunOffAvg, TTRunOffVar
   real      ::   TTSatFlowAvg,TTSatFlowVar
   real      ::   TTOutFlowAvg,TTOutFlowVar
   real      ::   TTWetAvg,    TTWetVar
   real      ::   TTWiniAvg,   TTWiniVar
   real(4)   ::   TArea
   
   integer   ::    NPixels			!/* Number of pixels in the basin */
   integer   ::    i,j,k 			!/* counters */
   real(4),allocatable  ::  Buff01(:)

!---------------------------------------------------------------------
!  Variáveis de depuração - Debug variables
!---------------------------------------------------------------------
   integer   ::   istatus            
!
!#####################################################################
! 
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!
!#####################################################################
!
!   write(6,'(a)')    '/-----------------------------------------------/   '
!   write(6,'(1x,a)') '/---- Entrando na rotina Aggregate--------------/'
!
! --------------------------------------------------------------------   
!  Allocate memory - Aloca memória
! --------------------------------------------------------------------
   if (.not. allocated(Buff01)) ALLOCATE(Buff01 (nx*ny),STAT=istatus)
   CALL check_alloc_status(istatus, "ephst:Buf01")
   Buff01 = 0
! --------------------------------------------------------------------   

   NPixels    = 0
   TTRnflxAvg    = 0.
   TTShflxAvg    = 0.
   TTLhflxAvg    = 0.
   TTGhflxAvg    = 0.
   TTPrecipAvg   = 0.
   TTRunOffAvg   = 0.
   TTSatFlowAvg  = 0.
   TTOutflowAvg  = 0.
   TTWetAvg      = 0.
   TTWiniAvg     = 0.
   TArea         = 0.
   Buff01        = 0.
   
!	/* aggregate the fluxes data */
!---------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx
      if (maskmap(i,j) == 1) then
	 if (isnan(lhflx(i,j)) == .false.) then
	    NPixels = NPixels + 1
	    TTRnflxAvg = TTRnflxAvg + rnflx(i,j)
	    Buff01(NPixels) = rnflx(i,j)
	 endif
      endif
   enddo
   enddo
   TTRnflxAvg = TTRnflxAvg / NPixels
   call stat_var( NPixels, Buff01, TTRnflxAvg, TTRnflxVar )	
	 
   Buff01  = 0.
   NPixels = 0
   do j = 1, ny
   do i = 1, nx	
      if (maskmap(i,j) == 1) then
	 if (isnan(lhflx(i,j)) == .false.) then
	    NPixels = NPixels + 1
	    TTLhflxAvg = TTLhflxAvg + lhflx(i,j)
	    Buff01(NPixels) = lhflx(i,j)
	 endif
      endif
   enddo
   enddo
   TTLhflxAvg = TTLhflxAvg / NPixels
   call stat_var( NPixels, Buff01, TTLhflxAvg, TTLhflxVar )

   Buff01  = 0.
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then   
	 if (isnan(lhflx(i,j)) == .false.) then
            TTShflxAvg = TTShflxAvg +  shflx(i,j)
	    Buff01(NPixels) = shflx(i,j)
	 endif
      endif
   enddo
   enddo
   TTShflxAvg = TTShflxAvg / NPixels   
   call stat_var( NPixels, Buff01, TTShflxAvg, TTShflxVar )	
  
   Buff01  = 0.   
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
	 if (isnan(lhflx(i,j)) == .false.) then
            TTGhflxAvg = TTGhflxAvg +  ghflx(i,j)
	    Buff01(NPixels) = ghflx(i,j)
	 endif
      endif
   enddo
   enddo 
   TTGhflxAvg = TTGhflxAvg / NPixels
   call stat_var( NPixels, Buff01, TTGhflxAvg, TTGhflxVar )	  
	    
!	/* aggregate precipitation data */
!---------------------------------------------------------------------
   NPixels = 0
   Buff01  = 0.
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
         NPixels = NPixels + 1
         TTPrecipAvg = TTPrecipAvg + precip(i,j)*(dtsfc/1000)
	 Buff01(NPixels) = precip(i,j)*(dtsfc/1000)
      endif 
   enddo
   enddo 
   TTPrecipAvg = TTPrecipAvg / NPixels
   call stat_var( NPixels, Buff01, TTPrecipAvg, TTPrecipVar )	  

!	/* aggregate evaporation data */
!---------------------------------------------------------------------
   NPixels       = 0
   Buff01  = 0.   
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
      	 if (isnan(lhflx(i,j)) == .false.) then
	    NPixels = NPixels + 1	 
            TTEvapAvg = TTEvapAvg + (Lhflx(i,j)/lathv)*(dtsfc/1000)
	    Buff01(NPixels) = (Lhflx(i,j)/lathv)*(dtsfc/1000)
	 endif
      endif 
   enddo
   enddo 
   TTEvapAvg  = TTEvapAvg  / NPixels
   call stat_var( NPixels, Buff01, TTEvapAvg, TTEvapVar )	  

!	/* aggregate RunOff data */
!---------------------------------------------------------------------
   NPixels = 0
   Buff01  = 0.
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
         NPixels = NPixels + 1	 
	 TTRunOffAvg   = TTRunOffAvg  + RunOff(i,j)
	 Buff01(NPixels) = RunOff(i,j)
      endif
   enddo
   enddo
   TTRunOffAvg   = TTRunOffAvg   / NPixels   
   call stat_var( NPixels, Buff01, TTRunOffAvg, TTRunOffVar )	  

!	/* aggregate subflowdata */
!---------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
         TTSatFlowAvg = TTSatFlowAvg + (-SatFlow(i,j))
	 Buff01(NPixels) = (-SatFlow(i,j))
      endif
   enddo
   enddo
   TTSatFlowAvg   = TTSatFlowAvg   / NPixels   
   call stat_var( NPixels, Buff01, TTSatFlowAvg, TTSatFlowVar )	  

!	/* aggregate outflow */
!---------------------------------------------------------------------	 
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
         TTOutflowAvg = TTOutflowAvg + OutFlow(i,j)
	 Buff01(NPixels) = OutFlow(i,j)
      endif
   enddo
   enddo
   TTOutFlowAvg   = TTOutFlowAvg   / NPixels   
   call stat_var( NPixels, Buff01, TTOutFlowAvg, TTOutFlowVar )	  	 

!	/* aggregate soil moisture data */
!---------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
	 TTWetAvg  = TTWetAvg  + (wetsfc(i,j)*(d1)+wetdp2(i,j)*(d2-d1)&
	             )/2 + wetdp3(i,j)*(d3-d2) + wetcnp(i,j)*dr
	 Buff01(NPixels) = (wetsfc(i,j)*(d1)+wetdp2(i,j)*(d2-d1)      &
	             )/2 + wetdp3(i,j)*(d3-d2) + wetcnp(i,j)*dr
      endif
   enddo
   enddo
   TTWetAvg   = TTWetAvg   / NPixels   
   call stat_var( NPixels, Buff01, TTWetAvg, TTWetVar )	  	 

!	/* aggregate initial soil moisture data */
!---------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx	 
      if (maskmap(i,j) == 1) then
         TTWiniAvg = TTWiniAvg + Wini(i,j)
	 Buff01(NPixels) = Wini(i,j)
      endif
   enddo
   enddo
   TTWiniAvg   = TTWiniAvg   / NPixels   
   call stat_var( NPixels, Buff01, TTWiniAvg, TTWiniVar )	  	 
  
!  /* Area media da bacia - average basin area
! --------------------------------------------------------------------   
   TArea = NPixels * (dx*dy)
!   write(6,'(a,f18.4,6f18.10)') 'TArea= ',TArea,TTPrecipAvg,TTEvapAvg, &
!                            TTRunOffAvg,TTSatFlowAvg,TTOutflowAvg,    &
!			    OutFlow(411,47)
! --------------------------------------------------------------------

end subroutine Aggregate


!
! -----------------------------------------------------------------------
!   --> stat_f4_2: calculates the average and variance of an array
! -----------------------------------------------------------------------
!
   SUBROUTINE stat_var( n, x, xavg, xvar )

      integer ::  n          ! # of elements
      real    ::  x(n)          ! data array
      real    ::  xavg          ! average 
      real    ::  xvar          ! variance

      integer ::  int 
      integer ::  i        ! counter

      real    ::  fn
      real    ::  desv
      real    ::  s1
      real    ::  s2c

!      write (6, '(a)')    ' Entrando em stat! '
!      write (6, '(3x,a,i10,a)')    'n = ', n,','
!      write (6, '(3x,a,f15.3,a)') 'x = ', x(2),','


      s2c = 0.0
! -------------------------------------------------------------------------
!   loop to calculate s2c
! -------------------------------------------------------------------------
!
      do  j = 1, n
         desv = x(j) - xavg
         s2c = s2c + (desv*desv)
      enddo

!
! -------------------------------------------------------------------------
!     obtains standard deviation
! -------------------------------------------------------------------------
!
      xvar = s2c / n   
!      write (6, '(3x,a,f15.3,a)') 's2c = ', s2c,','
!      write (6, '(3x,a,f15.3,a)') 'xvar = ', xvar,','

      return

      END SUBROUTINE stat_var  ! stat_var
