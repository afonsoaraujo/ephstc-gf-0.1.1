!##################################################################
!##################################################################
!######                                                      ######
!######                     MODULE OUTPUT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!module output

!interfaces
!interface IOutput
!   module procedure 
!   module procedure 
!end interface   

!contains

! --------------------------------------------------------------------
!  EphstcOut()
!
! --------------------------------------------------------------------
subroutine EphstcOut(it,nx,ny,maskmap,topomap,slopemap,aspectmap,     &
                    fgradmap,soilmap,sdepthmap,vegmap,ndvimap,laimap, &
		    vegfrac,tair,qvair,pres,relh,wspd,                &
                    rsi,rsr,precip,rnflx,shflx,lhflx,ghflx,Eg,Etr,Er, &
		    qvsatts,wetsfc,wetdp2,wetdp3,wetcnp,Wini,tsfc,    &
		    tdeep,Runoff,RunOffei,RunOffes,RunOffCnp,         &
                    SatFlow,Outflow,wtdepthmap,dirmap,Q,QOBS,RnflxObs,&
		    LhflxObs,ShflxObs,GhflxObs,TTRnflxAvg,TTLhflxAvg, &
		    TTShflxAvg,TTGhflxAvg,TTPrecipAvg,TTEvapAvg,      &
		    TTRunOffAvg,TTSatFlowAvg,TTOutflowAvg,TTWetAvg,   &
		    TTWiniAvg,TTRnflxVar,TTLhflxVar,TTShflxVar,       &
		    TTGhflxVar,TTPrecipVar,TTEvapVar,TTRunOffVar,     &
		    TTSatFlowVar,TTOutflowVar,TTWetVar,TTWiniVar,TArea)
!
!#######################################################################
!
!  Declaração de variáveis
!
!#######################################################################
!  
   implicit none 
   include 'globcst.inc'
   include 'grid.inc'
   include 'calendario.inc'
   
   integer   ::   nx,ny,it
   integer   ::   maskmap(nx,ny)
   real      ::   topomap(nx,ny)
   real      ::   slopemap(nx,ny)
   real      ::   aspectmap(nx,ny)
   real      ::   fgradmap(nx,ny)
   real      ::   soilmap(nx,ny)
   real      ::   sdepthmap(nx,ny)
   real      ::   vegmap(nx,ny)
   real      ::   ndvimap(nx,ny)
   real      ::   laimap(nx,ny)
   real      ::   vegfrac(nx,ny)

! Variáveis 2D
!-----------------------------------------------------------------------   
   real      ::   tair(nx,ny)
   real      ::   qvair(nx,ny)
   real      ::   pres(nx,ny)
   real      ::   relh(nx,ny)
   real      ::   wspd(nx,ny)
   real      ::   rsi(nx,ny)
   real      ::   rsr(nx,ny)
   real      ::   precip(nx,ny)
   real      ::   rnflx(nx,ny)
   real      ::   shflx(nx,ny)
   real      ::   lhflx(nx,ny)
   real      ::   Ghflx(nx,ny)
   real      ::   Eg(nx,ny)
   real      ::   Etr(nx,ny)
   real      ::   Er(nx,ny)
   real      ::   qvsatts(nx,ny)
   real      ::   wetsfc(nx,ny)
   real      ::   wetdp2(nx,ny)
   real      ::   wetdp3(nx,ny)
   real      ::   wetcnp(nx,ny)
   real      ::   Wini(nx,ny)
   real      ::   tsfc(nx,ny)
   real      ::   tdeep(nx,ny)
   real      ::   Runoff(nx,ny)
   real      ::   RunOffei(nx,ny)
   real      ::   Runoffes(nx,ny)
   real      ::   RunoffCnp(nx,ny)      
   real      ::   OutFlow(nx,ny)
   real      ::   SatFlow(nx,ny)
   real      ::   wtdepthmap(nx,ny)
   integer   ::   dirmap(nx,ny,NDIRS)   
   real      ::   Q
   real      ::   QOBS
   real      ::   RnflxObs(nx,ny)
   real      ::   LhflxObs(nx,ny)
   real      ::   ShflxObs(nx,ny)
   real      ::   GhflxObs(nx,ny)
   
! Variáveis 1D (Agregadas)
!-----------------------------------------------------------------------   
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
   real      ::   TArea
   real      ::   TPrecip
   real      ::   TEvap
   real      ::   TRunOff
   real      ::   TWet
   real      ::   TWini
   
! Direção do escoamento
!-----------------------------------------------------------------------
   real,allocatable      :: DirU(:,:)
   real,allocatable      :: DirV(:,:)   
   integer,allocatable   :: DirM(:,:)
   real , allocatable    :: DirMax(:,:)
   
! Balanço de Massa
!-----------------------------------------------------------------------   
   real,allocatable      :: MassBalance(:,:)
   real,allocatable      :: EBObs(:,:)
   real,allocatable      :: DeltaW(:,:)
      
! Diversos   
!-----------------------------------------------------------------------
   integer(4)   ::   i,j,aux,aux1,istatus
   integer(4)   ::   k,irec, sfcoutopt
   integer(4)   ::   medopt
   LOGICAL      ::   firstcall        ! First call flag of this subroutine
   character(138)  ::  sfcoutfl
   SAVE firstcall
   DATA firstcall/.true./   
!
!#######################################################################
!
!  Include:
!
!#######################################################################
!
   include 'phycst.inc'
   include 'soilcst.inc'
   include 'stations.inc'
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   namelist /ephstc_output/ flx2dout,flx1dout,w2dout,w1dout,mbout,    &
                            avgfout,topmodelout,                      &
			    nflxstations,nWstations,flx2doutfile,     &
                            flx1doutfile,ephstchkfile,mbfile,avgffile,&
			    topmodelfile
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(a)')    '/-----------------------------------------------/'
   write(6,'(1x,a)') '/---- Escrevendo as saidas do modelo!-----------/'   
!
! --------------------------------------------------------------------   
!  Lê lista de nomes - Read namelists
! --------------------------------------------------------------------
   REWIND (UNIT=5)
   READ (5, ephstc_output)
! --------------------------------------------------------------------   
!  Allocate memory - Aloca memória
! --------------------------------------------------------------------
   if (.not. allocated(MassBalance)) ALLOCATE(MassBalance (nx,ny),STAT=istatus)
   CALL check_alloc_status(istatus, "ephstc:MassBalance")
   MassBalance = 0
   
   if (.not. allocated(EBObs)) ALLOCATE(EBObs (nx,ny),STAT=istatus)
   CALL check_alloc_status(istatus, "ephstc:EBObs")
   EBObs = 0
   
   if (.not. allocated(DeltaW)) ALLOCATE(DeltaW (nx,ny),STAT=istatus)
   CALL check_alloc_status(istatus, "ephstc:DeltaW")
   DeltaW = 0

! --------------------------------------------------------------------   
!  Primeira chamada
! --------------------------------------------------------------------      
   if (firstcall) then
      if (flx1dout == 1) then   ! fluxos verticais 1d
         aux  = len_trim(flx1doutfile)
         do k = 1, nflxstations
            write(flx1doutfile(aux+1:aux+12),'(i3.3)') flxstation(k)
            open(20+k,file=flx1doutfile,form='formatted',status='unknown')
         enddo
      endif
      
      if (mbout == 1) then    ! Balanco de massa       
         aux1 = len_trim(mbfile)
         do k = 1, nflxstations
            write(mbfile(aux1+1:aux1+12),'(i3,4i2.2)') flxstation(k), &
                      hora, diam,mes,(ano-1900)                      
            open(40+k,file=mbfile,form='formatted',status='unknown')
	 enddo   
      endif
      
      if (avgfout == 1)  then  ! Geracao de campos medios
         open(12,file=avgffile,form='formatted',status='unknown') 
      endif
      
      if (topmodelout == 1) then  ! Arquivo para inicialiacao externa do TOPMODEL 
            open(13,file=topmodelfile,form='formatted',status='unknown')
      endif
           
      firstcall = .false.
   end if
! --------------------------------------------------------------------   
!  Balanço de massa em cada célula de grade - Mass balance
!  Repare que as unidades estão todas em [m].
! --------------------------------------------------------------------      
   if (mbout == 1) then
      do j = 1, ny
      do i = 1, nx
         DeltaW(i,j) = ((wetsfc(i,j)*d1 + wetdp2(i,j)*(d2-d1))/2 +    &
	                 wetdp3(i,j)*(d3-d2) + wetcnp(i,j)*dr) -      &
			 Wini(i,j)
			 
! Unidades de volume/area [m3/m2], ou seja em cada ponto de grade
          MassBalance(i,j) = (Wini(i,j) + precip(i,j) - Eg(i,j) -     &
	                      Etr(i,j) - Er(i,j) + DeltaW(i,j)) *     &
			     (dtsfc/rhow) - RunOff(i,j) 
	 			     
      enddo
      enddo
   
      do k = 1, nflxstations                       
         i=flxpoints(k,1)
         j=flxpoints(k,2)
         write(40+k,'(11f18.12)') wetdp2(i,j)*d2,wetdp3(i,j)*(d3-d2), &
	                          wetcnp(i,j)*dr, precip(i,j)*(dtsfc  &
				  /rhow),Eg(i,j)*(dtsfc/rhow),Etr(i,j)&
				  *(dtsfc/rhow),Er(i,j)*(dtsfc/rhow), &
				  Wini(i,j),DeltaW(i,j),              &
				  MassBalance(i,j)
      enddo
   endif
! --------------------------------------------------------------------   
!  Saídas 1D de variáveis 2D - formato Ascii
! --------------------------------------------------------------------
   if (flx1dout == 1) then
   
      do k = 1, nflxstations
         i=flxpoints(k,1)
         j=flxpoints(k,2)

         write(20+k,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2,32f16.5)')&
                 ano,mes,diam,diaj,hora,mint,                            &
                 rnflx(i,j),lhflx(i,j),shflx(i,j),Ghflx(i,j),            &
                 rsi(i,j),rsr(i,j),(lathv*Eg(i,j)),(lathv*Etr(i,j)),     &
                 (lathv*Er(i,j)),wetsfc(i,j),wetdp2(i,j),wetdp3(i,j),    &
                 tair(i,j),tsfc(i,j),tdeep(i,j),                         &
		 precip(i,j)*(-1)*(dtsfc),                               &
                 RunOff(i,j),RunOffes(i,j),                              &
                 RunOffei(i,j)*(dtroutsfc/1000),                         &		 
                 RunOffCnp(i,j)*(dtroutsfc/1000),                        &
                 SatFlow(i,j),wtdepthmap(i,j),qvair(i,j),qvsatts(i,j),   &
		 wetcnp(i,j),OutFlow(i,j)
!-----------------------------------------------------------------------
! GERA FORCANTES 1D 
!-----------------------------------------------------------------------
!         write(50+k,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2,22f16.5)')&
!                 ano,mes,diam,diaj,hora,mint,                            &
!		 tair(i,j),pres(i,j),relh(i,j),wspd(i,j),rsi(i,j),       &
!		 rsr(i,j),precip(i,j)
!-----------------------------------------------------------------------		 
      enddo

   endif

! --------------------------------------------------------------------   
!  Saída 1D de variáveis Agregadas - formato Ascii
! --------------------------------------------------------------------
!   avgfout = 1
   if (avgfout == 1) then

!  unidades convertidas em aggregate para ===> [m]
      write(12,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2,21f19.12)') &
              ano,mes,diam,diaj,hora,mint,TTRnflxAvg,sqrt(TTRnflxVar),&
	      TTLhflxAvg,sqrt(TTLhflxVar),TTShflxAvg,sqrt(TTShflxVar),&
	      TTGhflxAvg,sqrt(TTGhflxVar),TTPrecipAvg,                &
	      sqrt(TTPrecipVar),TTEvapAvg,sqrt(TTEvapVar),TTRunOffAvg,&
	      sqrt(TTRunOffVar),TTWetAvg,sqrt(TTWetVar),TTWiniAvg,    &
	      sqrt(TTWiniVar),TTOutFlowAvg*(TArea/dtroutsfc),         &
	      OutFlow(411,47)*(TArea/dtroutsfc)
	      
   endif
   

!  Saída 1D para inicializacao externa do TOPMODEL - formato Ascii
!  unidades em [m/h]
! --------------------------------------------------------------------
   if (topmodelout == 1) then  
!      write(6,'(a,2f18.8,i4)') '1- ',QOBS,TArea,dtroutsfc
      if (isnan(TTEvapAvg) .eqv. .false.) then             
         write(13,'(3f18.8)') TTEvapAvg, TTPrecipAvg, QOBS *          &
	                      dtroutsfc/TArea
      else 
         write(13,'(3f18.8)') 0.0000, 0.0000, QOBS * dtroutsfc/       &
	                       TArea
      endif   
   endif

! --------------------------------------------------------------------   
!  Saída 2D - formato Grads 
! --------------------------------------------------------------------
   if (flx2dout == 1) then

!     Calcula componentes meridional e zonal da direção do escoamento
! --------------------------------------------------------------------   
!      call FluxDir (nx,ny,maskmap,dirmap,SatFlow,DirU,DirV,DirM)
!      DirMax = DirM
!      do j=1,ny
!      do i=1,nx
!         write(6,'(a,f18.8,i8)') 'DirMax= ', Dirmax(i,j), DirM(i,j)
!      enddo
!      enddo  
! --------------------------------------------------------------------   
      aux  = len_trim(flx2doutfile)
      write(flx2doutfile(aux+1:aux+10),'(i4.4,3i2.2)') ano,mes, &
            diam,hora
      
! --------------------------------------------------------------------
!     Abre arquivos de saída e verificação
! --------------------------------------------------------------------
      open(11,file=flx2doutfile,form='unformatted',access='direct',   &
           recl=nx*ny,status='unknown')

! --------------------------------------------------------------------
!     Conversões de unidades para visualização...
! --------------------------------------------------------------------
!      precip(i,j)   = precip(i,j)*(-1)*dtsfc    ! [kg/s.m2] ==> [mm/h]
!      RunOff(i,j)   = RunOff(i,j)*dtroutsfc     ! [kg/s.m2] ==> [mm/h]
!      RunOffes(i,j) = RunOffes(i,j)*dtroutsfc   ! [kg/s.m2] ==> [mm/h]
!      RunOffei(i,j) = RunOffei(i,j)*dtroutsfc   ! [kg/s.m2] ==> [mm/h]
!      RunOffCnp(i,j)= RunOffCnp(i,j)*dtroutsfc  ! [kg/s.m2] ==> [mm/h]
!      SatFlow(i,j)  = SatFlow(i,j)*(-1)*1000    ! [m]       ==> [mm]
! --------------------------------------------------------------------
      irec = 1
      write(11,rec=irec) ((rnflx(i,j),i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((precip(i,j)*dtsfc,i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((lhflx(i,j),i=1,nx),j=1,ny)            
      irec=irec+1
      write(11,rec=irec) ((shflx(i,j),i=1,nx),j=1,ny)            
      irec=irec+1
      write(11,rec=irec) ((Ghflx(i,j),i=1,nx),j=1,ny)            
      irec=irec+1
      write(11,rec=irec) ((wetsfc(i,j),i=1,nx),j=1,ny)            
      irec=irec+1
      write(11,rec=irec) ((wetdp2(i,j),i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((wetdp3(i,j),i=1,nx),j=1,ny)      
      irec=irec+1
      write(11,rec=irec) ((wtdepthmap(i,j),i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((RunOff(i,j)*dtroutsfc,i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((RunOffes(i,j)*dtroutsfc,i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((RunOffei(i,j)*dtroutsfc,i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((SatFlow(i,j)*1000,i=1,nx),j=1,ny)
      irec=irec+1
      write(11,rec=irec) ((OutFlow(i,j)*1000,i=1,nx),j=1,ny)
      irec=irec+1
   endif
 

! --------------------------------------------------------------------   
!  Gera campos médios sobre a bacia
! --------------------------------------------------------------------
!   medopt = 0
!   if (medopt == 1) then
!      
!      evapm   = 0.0
!      precipm = 0.0
!      k = 0
!      write(6,'(a,2i4,f16.8)') 'it = ', it,dtforc,QOBS
!      do j = 1, ny
!      do i = 1, nx
!         if (maskmap(i,j) == 1) then
!	     ArB = ArB +1
!	  endif1
!
!         if (maskmap(i,j) == 1 .and. lhflx(i,j) .gt. 0.0) then
!             
!            if (isnan(lhflx(i,j)) == .false.) then
!               evapm   = evapm + (lhflx(i,j)/lathv)*(3600/1000)
!               precipm = precipm + precip(i,j)*(3600/1000)
!               k = k + 1
!            endif
!            write(6,'(4i6,4f18.8)') maskmap(i,j),it,i,j,lhflx(i,j),precip(i,j),evapm,precipm
!         endif
!      enddo
!      enddo
!      
!      ArB = ArB*(dx*dy)
!      write(6,'(a,f18.6)') 'Ab= ', ArB
!      evapm = evapm/k
!      precipm = precipm/k
!      
!      if (isnan(evapm) == .false.) then             
!         write(13,'(3f18.8)') evapm, precipm, QOBS*dtroutsfc/ArB
!      else 
!         write(13,'(3f18.8)') 0.0000, 0.0000, QOBS*dtroutsfc/ArB
!      endif
!      
!   endif
   
! --------------------------------------------------------------------   
!  Fecha arquivos de saida
! --------------------------------------------------------------------
   if (it == nsteps) then
      if (flx1dout == 1) then
         do k = 1, nflxstations   
	    close (20+k)
	    close (40+k)
	    close (110+k)
	 enddo
      endif
   endif
          
end subroutine EphstcOut


! --------------------------------------------------------------------
!  FluxDir()
!
! --------------------------------------------------------------------
subroutine FluxDir (nx,ny,maskmap,dirmap,SatFlow,DirU,DirV,DirM)
!
!#######################################################################
!
!  Declaração de variáveis
!
!#######################################################################
!  
   implicit none 

   include 'grid.inc'
   
   integer   ::   nx,ny,it
   integer   ::   maskmap(nx,ny)
   
   integer   ::   dirmap(nx,ny,NDIRS)
   real      ::   SatFlow(nx,ny)
   real      ::   DirU   (nx,ny)
   real      ::   DirV   (nx,ny)
   integer   ::   DirM   (nx,ny)   
   
   integer   ::   i,j,k
   integer   ::   dirmaxx
!
!#######################################################################
!
!  Variáveis de depuração - Debug variables
!
!#######################################################################
!
   integer   ::   istatus
!
!#######################################################################
!
!  Include:
!
!#######################################################################
!
   include 'globcst.inc'
   include 'phycst.inc'
   include 'stations.inc'
!
!#######################################################################
! 
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!
!#######################################################################
!
   write(6,'(a)') '   '
   write(6,'(1x,a)') ' Entrando na rotina FluxDir...'
!
!#######################################################################
!
!  Aloca memória - Memory alloc
!
!#######################################################################
!
!  ALLOCATE(DirMax    (nx,ny),STAT=istatus)
!  CALL check_alloc_status(istatus, "ephstc:DirMax")
!  DirMax = 0
!-----------------------------------------------------------------------
!  /*  */
!-----------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx
   
      if (maskmap(i,j) .ne. 1 ) then
      
         DirU(i,j)   = -9999
         DirV(i,j)   = -9999
         DirM(i,j)   = -9999
      else

         dirmaxx = 0
         do k = 1, NDIRS
            if (dirmap(i,j,k) > dirmaxx) then
               DirM(i,j) = k
               dirmaxx = dirmap(i,j,k)	       
            endif
         enddo
         
         if (DirM(i,j) == 1) then
            DirU(i,j) = 0.0
            DirV(i,j) = abs(SatFlow(i,j))
         else if (DirM(i,j) == 2) then
            DirU(i,j) = abs(SatFlow(i,j) * sin(0.7854))
            DirV(i,j) = abs(SatFlow(i,j) * cos(0.7854))
         else if (DirM(i,j) == 3) then
            DirU(i,j) = abs(SatFlow(i,j))
            DirV(i,j) = 0.0
         else if (DirM(i,j) == 4) then
            DirU(i,j) =  abs(SatFlow(i,j) * cos(0.7854))
            DirV(i,j) = -abs(SatFlow(i,j) * sin(0.7854))
         else if (DirM(i,j) == 5) then
            DirU(i,j) = 0.0
            DirV(i,j) = -abs(SatFlow(i,j))
         else if (DirM(i,j) == 6) then
            DirU(i,j) = -abs(SatFlow(i,j) * cos(0.7854))
            DirV(i,j) = -abs(SatFlow(i,j) * sin(0.7854))
         else if (DirM(i,j) == 7) then
            DirU(i,j) = 0.0
            DirV(i,j) = -abs(SatFlow(i,j))
         else if (DirM(i,j) == 8) then
            DirU(i,j) =  abs(SatFlow(i,j) * cos(0.7854))
            DirV(i,j) = -abs(SatFlow(i,j) * sin(0.7854))
         endif
     
      endif   
         
   enddo
   enddo

return

end subroutine FluxDir


!end module output
