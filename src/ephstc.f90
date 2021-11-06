!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                    PROGRAM EPHSTC                    ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
PROGRAM EPHSTC
!
!#######################################################################
!
!  Declaração de variáveis - Variable declarations:
!
!#######################################################################
!   
!use init          ! módulo para inicialização das variáveis
!use calendario    ! funções de calendário
!use derivs
!use route         ! módulo para propagação do escoamento lateralmente
!use output        ! rotinas para geração de saídas do modelo
!
!#######################################################################
!
!  Include files
!
!#######################################################################
!
   implicit none

   include 'globcst.inc'     ! Global constants that control model
                             ! execution
   include 'soilcst.inc'
   include 'calendario.inc' ! Define as variáveis de calendário
   include 'grid.inc'
   include 'alloc.inc'      ! Memory allocation declaration and interfaces
   include 'topmodel.inc'
!
!#######################################################################
!
!  Estruturas de dados - Data structures
!
!#######################################################################
!   
   integer   ::   nx,ny,nzsoil

   integer,allocatable   ::   maskmap(:,:)   
   real,   allocatable   ::   topomap(:,:)
   real,   allocatable   ::   aspectmap(:,:)
   real,   allocatable   ::   slopemap(:,:)
   real,   allocatable   ::   lnatbmap(:,:)
   integer,allocatable   ::   dirmap(:,:,:)
   integer,allocatable   ::   totaldirmap(:,:)
   real,   allocatable   ::   fgradmap(:,:)
   real,   allocatable   ::   hidromap(:,:)

   integer,allocatable   ::   soilmap(:,:)
   real,   allocatable   ::   sdepthmap(:,:)
   integer,allocatable   ::   vegmap(:,:)
   real,   allocatable   ::   vegfrac(:,:)
   real,   allocatable   ::   ndvimap(:,:)
   real,   allocatable   ::   laimap(:,:)
   real,   allocatable   ::   wtdepthmap(:,:)
   real,   allocatable   ::   wtlevelmap(:,:)
   real,   allocatable   ::   tair(:,:)
   real,   allocatable   ::   pres(:,:)
   real,   allocatable   ::   relh(:,:)
   real,   allocatable   ::   wspd(:,:)
   real,   allocatable   ::   rsi(:,:)
   real,   allocatable   ::   rsr(:,:)
   real,   allocatable   ::   precip(:,:)
  
   real,   allocatable   ::   tsfc(:,:)
   real,   allocatable   ::   tdeep(:,:)
   real,   allocatable   ::   wetsfc(:,:)
   real,   allocatable   ::   wetdp2(:,:)
   real,   allocatable   ::   wetdp3(:,:)   
   real,   allocatable   ::   wetcnp(:,:)
   real,   allocatable   ::   snowdpth(:,:)
   real,   allocatable   ::   qvsatta(:,:)
   real,   allocatable   ::   qvsatts(:,:)
   real,   allocatable   ::   qvsfc(:,:)
   real,   allocatable   ::   qvair(:,:)
   real,   allocatable   ::   rhoa(:,:)
   
   real,   allocatable   ::   cdha(:,:)
   real,   allocatable   ::   cdqa(:,:)   
   real,   allocatable   ::   f34(:,:)
   
   real,   allocatable   ::   ct(:,:)
   real,   allocatable   ::   rnflx(:,:)
   real,   allocatable   ::   shflx(:,:)
   real,   allocatable   ::   lhflx(:,:)
   real,   allocatable   ::   ghflx(:,:)
   real,   allocatable   ::   Eg(:,:)
   real,   allocatable   ::   Etr(:,:)
   real,   allocatable   ::   Er(:,:)
   real,   allocatable   ::   RunOff(:,:)
   real,   allocatable   ::   RunOffei(:,:)   
   real,   allocatable   ::   RunOffes(:,:)   
   real,   allocatable   ::   RunOffCnp(:,:)
   real,   allocatable   ::   SatFlow(:,:)
   real,   allocatable   ::   OutFlow(:,:)
   real,   allocatable   ::   Wini(:,:)
!-----------------------------------------------------------------------
! Variáveis Agregadas - Aggregate variables
!-----------------------------------------------------------------------
   real                  ::   TTRnflxAvg,  TTRnflxVar
   real                  ::   TTShflxAvg,  TTShflxVar
   real                  ::   TTLhflxAvg,  TTLhflxVar
   real                  ::   TTGhflxAvg,  TTGhflxVar
   real                  ::   TTPrecipAvg, TTPrecipVar
   real                  ::   TTEvapAvg,   TTEvapVar
   real                  ::   TTRunOffAvg, TTRunOffVar
   real                  ::   TTSatFlowAvg,TTSatFlowVar
   real                  ::   TTOutFlowAvg,TTOutFlowVar
   real                  ::   TTWetAvg,    TTWetVar
   real                  ::   TTWiniAvg,   TTWiniVar
   real                  ::   TArea
   real                  ::   TPrecip
   real                  ::   TEvap
   real                  ::   TRunOff
   real                  ::   TWet
   real                  ::   TWini
!-----------------------------------------------------------------------
! Variáveis para a atualização da linha d'água
!-----------------------------------------------------------------------
   real                  ::   NRootLayers
   real,   allocatable   ::   LPorosity(:)
   real,   allocatable   ::   LFCap(:)
   real,   allocatable   ::   RootDepth(:)
   real,   allocatable   ::   Moisture(:)
!-----------------------------------------------------------------------
   real(4)               ::   CalcWTableDepth
!   real(4)               ::   HeadSlopeAspect
! Variáveis do topmodel   
!-----------------------------------------------------------------------
   real,   allocatable   ::   Q(:)
   integer               ::   NAC,IA,NSC
   real                  ::   AC(NACMX)
   real                  ::   ST(NACMX)
   real                  ::   SUMAC
   real                  ::   TL,NCH
   real                  ::   ACH(NACMX)
   real                  ::   D(NACMX)
   real                  ::   AR(NCHMX)
   integer               ::   ND,NR   
! Variáveis observadas
!-----------------------------------------------------------------------   
   real,   allocatable   ::   RnflxObs(:)
   real,   allocatable   ::   LhflxObs(:)
   real,   allocatable   ::   ShflxObs(:)
   real,   allocatable   ::   GhflxObs(:)
   real,   allocatable   ::   QOBS(:)   
!
!#######################################################################
!
!  Variáveis gerais - 
!
!#######################################################################
!
   integer               ::   forcflopt
   character(128)        ::   tairfile,presfile,relhfile,wspdfile,    &
                              rsifile,rsrfile,precipfile
!
!#######################################################################
!
!  Variáveis locais - Local variables
!
!#######################################################################
!   
   integer   ::   i,j,k,tint,tintn,Layer
   integer   ::   bis,aux1,aux2,ii,jj
!
!#######################################################################
!
!  Variáveis de depuração - Debug variables
!
!#######################################################################
!
!   logical   ::   STAT_ASSOC
   integer   ::   istatus
!
!#######################################################################
!
!  Interfaces - Interfaces:
!
!#######################################################################
!

!
!#######################################################################
!
!  Listas de nomes - Namelists:
!
!#######################################################################
!
   namelist /ephstc_setup/   runname,initime,tstartsfc,tstopsfc,dtbig,  &
                             dtsfc,dtforc,dtroutsfc
   namelist /ephstc_mapgrid/ nx,ny,nzsoil,DX,DY,Xorig,Yorig,OffSetX,    &
                             OffSetY,nx_forc,ny_forc,DX_Forc,DY_Forc,   &
                             NDIRS	
   namelist /ephstc_forcfl/  forcflopt,tairfile,presfile,relhfile,      &
                             wspdfile,rsifile,rsrfile,precipfile               	
   namelist /ephstc_runopt/  ephstmode,epstopt,epelopt
   namelist /ephstc_soilebm/ moist
   namelist /ephstc_isba3l/  qpw
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Início do código executável
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!#######################################################################
!
!  Leitura das listas de nomes - Namelists reading:
!
!#######################################################################
!
!  Lê NAMELISTS 
   rewind(UNIT=5)
   read (5, ephstc_setup)
   read (5, ephstc_mapgrid)
   read (5, ephstc_forcfl)
   read (5, ephstc_runopt)   
   read (5, ephstc_soilebm)

   neq=6
   nsteps = (tstopsfc - tstartsfc   ) + 1
!
!#######################################################################
!
!  Aloca memória - Memory alloc
!
!#######################################################################
!
  ALLOCATE(maskmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:maskmap")
  maskmap = 0
  
  ALLOCATE(topomap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:topomap")
  topomap = 0
  
  ALLOCATE(slopemap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:slopmap")
  slopemap = 0
  
  ALLOCATE(aspectmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:aspectmap")
  aspectmap = 0
  
  ALLOCATE(lnatbmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:lnatbmap")
  lnatbmap = 0.

  ALLOCATE(fgradmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:fgradmap")
  fgradmap = 0.

  ALLOCATE(dirmap    (nx,ny,NDIRS),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:dirmap")
  dirmap = 0
  
  ALLOCATE(totaldirmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:totaldirmap")
  totaldirmap = 0  

  ALLOCATE(hidromap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:hidromap")
  hidromap = 0

  ALLOCATE(soilmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:soilmap")
  soilmap = 0
  
  ALLOCATE(sdepthmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:sdepthmap")
  sdepthmap = 0
  
  ALLOCATE(vegmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:vegmap")
  vegmap = 0

  ALLOCATE(vegfrac    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:vegfrac")
  vegfrac = 0

  ALLOCATE(ndvimap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:ndvimap")
  ndvimap = 0
  
  ALLOCATE(laimap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:laimap")
  laimap = 0
  
  ALLOCATE(wtdepthmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wtdepthmap")
  wtdepthmap = 0
  
  ALLOCATE(wtlevelmap    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wtlevelmap")
  wtlevelmap = 0  
  
  ALLOCATE(tair    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:tair")
  tair = 0
  
  ALLOCATE(wspd    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wspd")
  wspd = 0

  ALLOCATE(pres    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:pres")
  pres = 0
  
  ALLOCATE(relh    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:relh")
  relh = 0
    
  ALLOCATE(precip    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:precip")
  precip = 0
  
  ALLOCATE(rsi    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:rsi")
  rsi = 0
  
  ALLOCATE(rsr    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:rsr")
  rsr = 0

  ALLOCATE(tsfc    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:tsfc")
  tsfc = 0
  
  ALLOCATE(tdeep    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:tdeep")
  tdeep = 0
  
  ALLOCATE(wetsfc    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wetsfc")
  wetsfc = 0
  
  ALLOCATE(wetdp2    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wetdp2")
  wetdp2 = 0
  
  ALLOCATE(wetdp3    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wetdp3")
  wetdp3 = 0
  
  ALLOCATE(wetcnp    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:wetcnp")
  wetcnp = 0
  
  ALLOCATE(snowdpth    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:snowdpth")
  snowdpth = 0  

  ALLOCATE(qvsatta    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:qvsatta")
  qvsatta = 0
  
  ALLOCATE(qvsatts    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:qvsatts")
  qvsatts = 0
  
  ALLOCATE(qvsfc    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:qvsfc")
  qvsfc = 0
  
  ALLOCATE(qvair    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:qvair")
  qvair = 0  
  
  ALLOCATE(rhoa    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:rhoa")
  rhoa = 0
  
  ALLOCATE(cdqa    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:cdqa")
  cdqa = 0
  
  ALLOCATE(cdha    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:cdha")
  cdha = 0    
  
  ALLOCATE(f34    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:f34")
  f34 = 0  
  
  ALLOCATE(ct    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:ct")
  ct = 0
  
  ALLOCATE(rnflx    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:rnflx")
  rnflx = 0
  
  ALLOCATE(shflx    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:shflx")
  shflx = 0
  
  ALLOCATE(lhflx    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:lhflx")
  lhflx = 0
  
  ALLOCATE(ghflx    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:ghflx")
  ghflx = 0
  
  ALLOCATE (Eg    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Eg")
  Eg = 0  
  
  ALLOCATE(Etr    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Etr")
  Etr = 0
  
  ALLOCATE(Er    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Er")
  Er = 0  
  
  ALLOCATE(RunOff    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RunOff")
  RunOff = 0

  ALLOCATE(RunOffei    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RunOffei")
  RunOffei = 0
  
  ALLOCATE(RunOffes    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RunOffes")
  RunOffes = 0
  
  ALLOCATE(RunOffCnp    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RunOffCnp")
  RunOffCnp = 0
  
  ALLOCATE(SatFlow    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:SatFlow")
  SatFlow = 0    
  
  ALLOCATE(OutFlow    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:OutFlow")
  OutFlow = 0
  
  ALLOCATE(Q    (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Q")
  Q = 0      
  
  ALLOCATE(Wini    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Wini")
  Wini = 0      
!---------------------------------------------------------------------    
  if (.not. allocated(LPorosity)) ALLOCATE(LPorosity  (nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:LPorosity")
  LPorosity = 0
  
  if (.not. allocated(LFCap)) ALLOCATE(LFCap  (nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:LFCap")
  LFCap = 0 
  
  if (.not. allocated(RootDepth)) ALLOCATE(RootDepth (nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RootDepth")
  RootDepth = 0
  
  if (.not. allocated(Moisture)) ALLOCATE(Moisture (nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:Moisture")
  Moisture = 0
!---------------------------------------------------------------------    
  if (.not. allocated(RnflxObs)) ALLOCATE(RnflxObs (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:RnflxObs")
  RnflxObs = 0      

  if (.not. allocated(GhflxObs)) ALLOCATE(GhflxObs (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:GhflxObs")
  GhflxObs = 0  
  
  if (.not. allocated(ShflxObs)) ALLOCATE(ShflxObs (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:ShflxObs")
  ShflxObs = 0  
  
  if (.not. allocated(LhflxObs)) ALLOCATE(LhflxObs (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:LhflxObs")
  LhflxObs = 0  
  
  if (.not. allocated(QOBS)) ALLOCATE(QOBS (nsteps),STAT=istatus)
  CALL check_alloc_status(istatus, "ephstc:QOBS")
  QOBS = 0  
!
!#######################################################################
!
!  Processo de inicialização - Inicialization processes:
!
!#######################################################################
!  
!-----------------------------------------------------------------------   
!  Inicializa Calendário
!-----------------------------------------------------------------------   
   call InitCalend(initime,ano,mes,diam,diaj,hora,mint)
!-----------------------------------------------------------------------   
!  Inicializa MapGrid
!-----------------------------------------------------------------------   
!   call InitMapGrid(nx,ny,nzsoil,dx,dy,Xorig,Yorig,OffSetX,OffSetY    &
!                    nx_forc,ny_forc,DX_Forc,DY_Forc)

!-----------------------------------------------------------------------   
!  Inicializa os forçantes do modelo (ReadForc)
!-----------------------------------------------------------------------   
   call ReadForcantes(nx,ny,nx_forc,ny_forc,tairfile,presfile,        &
                      relhfile,wspdfile,rsifile,rsrfile,precipfile,   &
                      tair,pres,relh,wspd,rsi,rsr,precip) 
!-----------------------------------------------------------------------   
!  Inicializa o estado do modelo (Model State)
!-----------------------------------------------------------------------   
   call InitModelState (nx,ny,nzsoil,tair,pres,relh,tsfc,tdeep,       &
                        wetsfc,wetdp2,wetdp3,wetcnp,snowdpth,         &
                        qvsatta,qvsatts,qvsfc,rhoa,Wini)
!-----------------------------------------------------------------------    
!  Inicializa Mapas do terreno (TerrainMaps)
!-----------------------------------------------------------------------   
   call InitTerrainMaps(nx,ny,nzsoil,maskmap,topomap,aspectmap,       &
                        slopemap,lnatbmap,fgradmap,dirmap,totaldirmap,&
                        hidromap,soilmap,sdepthmap,vegmap,vegfrac,    &
                        ndvimap,laimap,wetsfc,wetdp2,wetdp3,          &
                        wtdepthmap,wtlevelmap)
!-----------------------------------------------------------------------
!  Inicializa Mapa de indice topografico			
!-----------------------------------------------------------------------	
   if (epelopt == 2) then
      call InitTopmData (TArea,NAC,IA,AC,ST,SUMAC,TL,NCH,ACH,D,ND,NR, &
                         AR)
   endif
!-----------------------------------------------------------------------      ! Inicializa parâmetros do Svat
! **implementar processo de calibracao aqui!
!-----------------------------------------------------------------------   
   call InitSvatpar (nx,ny,cdha,cdqa)
!-----------------------------------------------------------------------   
! Inicializa dados observados
!-----------------------------------------------------------------------   
!   call InitObsrvData (QOBS,RnflxOBS,ShflxObs,LhflxOBS,GhflxObs)
   
   write(6,'(a)') ' Processo de inicialização concluído! '

!#######################################################################
!
!  Processo de integração - Integration processes:
!
!#######################################################################
!-----------------------------------------------------------------------   
!  Integração do  ephstc - ephstc integration
!-----------------------------------------------------------------------   
!Informações para depuração
   write(6,'(a)') '  '
   write(6,'(a)') ' Iniciando processo de integração temporal! '
   write(6,'(a,1x,i6)') ' tstartsfc= ', tstartsfc
   write(6,'(a,1x,i6)') ' tstopsfc = ', tstopsfc
   
   nsfcst = INT(dtbig/dtsfc)
   do tint = tstartsfc, tstopsfc
   
!Informações para depuração
      write(6,'(a)') '#################################################'
      write(6,'(a,i4)') ' Tempo de integracao => ', tint
!-----------------------------------------------------------------------   
!  Ephstc - esquema de parametrização hidrológica da superfície terrestre
!-----------------------------------------------------------------------
      if (epstopt == 0) then
         call Calcqvair(nx,ny,pres,tair,relh,qvair)
         call CalcRnflx(nx,ny,vegmap,tair,tsfc,rsi,rsr,rnflx,qvair) 	 
         
         call soilebm (nx,ny,nzsoil,soilmap,vegmap,laimap,vegfrac,    &
                       tsfc,tdeep,wetsfc,wetdp2,wetcnp,snowdpth,      &
                       qvsfc,wspd,pres,rhoa,precip,                   &
                       tair,qvair,cdha,cdqa,rsi,                      &
                       rnflx,shflx,lhflx,ghflx,ct,Eg,Etr,Er,          &
                       qvsatts,qvsatta,f34)
         
      else if (epstopt == 1 .or. epstopt == 2) then



         call Calcqvair(nx,ny,pres,tair,relh,qvair)
         call CalcRnflx(nx,ny,vegmap,tair,tsfc,rsi,rsr,rnflx,qvair)
	 
         call Depst3l2d(nx,ny,nzsoil,soilmap,vegmap,laimap,vegfrac,   &
                         tsfc,tdeep,wetsfc,wetdp2,wetdp3,wetcnp,      &
                         snowdpth,qvsfc,wspd,pres,rhoa,precip,        &
                         tair,qvair,cdha,cdqa,rsi,                    &
                         rnflx,shflx,lhflx,ghflx,ct,Eg,Etr,Er,        &
                         qvsatts,qvsatta,f34,                         &
                         SatFlow,RunOffei,RunOffCnp,wtdepthmap,       &
                         sdepthmap)
         
	 if (epstopt == 1) then
	    do j = 1, ny
	    do i = 1, nx
	       RunOff(i,j) = RunOffei(i,j)*(dtsfc/1000)
	    enddo
	    enddo
	 endif

      endif
!-----------------------------------------------------------------------   
!     CalcWtDepth - Atualiza a linha d'água em função do movimento 
!                   vertival da água na célula
!-----------------------------------------------------------------------      
!---------------------------------------------------------------------
!  Se a redistribuição lateral do esocamento estiver ativa a profundidade
!  da linha d'água deve ser atualizada corretamente!
!
!  /* Calculate the depth of the water table based on the soil moisture 
!     profile and adjust the soil moisture profile, to assure that the 
!     soil moisture is never more than the maximum allowed soil moisture 
!     amount,i.e. the porosity.  A negative water table depth means that
!     the water is ponding on the surface.  This amount of water becomes 
!     surface Runoff */
!---------------------------------------------------------------------
      if (epstopt == 2 .or. epelopt == 1) then  ! Rever essas opções

         do j = 1, ny
         do i = 1, nx

            do k =1, nzsoil
               LPorosity(k) = porosity(soilmap(i,j))      ! constante ao longo de z
               LFCap(k) = wfc(soilmap(i,j))               ! constante ao longo de z
            enddo

            call CalcRootDepth(nzsoil,vegmap(i,j),NRootlayers,        &
                               RootDepth,sdepthmap(i,j),i,j) 
            Moisture(1) = wetsfc(i,j)
            Moisture(2) = wetdp2(i,j)
            Moisture(3) = wetdp3(i,j)

            wtdepthmap(i,j) = CalcWTableDepth(NRootLayers,sdepthmap(i,j),&
                                              RootDepth,LFCap,           &
                                              LPorosity,Moisture,i,j)

            wetsfc(i,j) = Moisture(1)
            wetdp2(i,j) = Moisture(2)
            wetdp3(i,j) = Moisture(3)

            RunOffes(i,j) = 0.  
            if (wtdepthmap(i,j) < 0.0) then

               RunOffes(i,j) = (-(wtdepthmap(i,j)))
               wtdepthmap(i,j) = 0.0
	       
            endif

         enddo
         enddo

      endif
!-----------------------------------------------------------------------
!     Epel - Esquema de parametrização do escoamento lateral
!     Propagação lateral do escoamento - Distribute subsurface water
!     laterally
!-----------------------------------------------------------------------
      if (epelopt == 1 .or. epelopt == 2) then

!  Calcula/atualiza os gradientes do escoamento
! --------------------------------------------------------------------  
         call HeadSlopeAspect(nx,ny,maskmap,topomap,wtlevelmap,       &
                              slopemap,aspectmap,fgradmap,dirmap,     &
                              totaldirmap)

!  Redistribui o escoamento lateralmente
! --------------------------------------------------------------------  
         call RoutSubSfc (dtroutsfc,nx,ny,nzsoil,maskmap,topomap,     &
	                  aspectmap,slopemap,fgradmap,dirmap,soilmap, &
                          vegmap,sdepthmap,wtdepthmap,wtlevelmap,     &
                          OutFlow,SatFlow)

!  Agrega dados - Aggregate data
! --------------------------------------------------------------------        !         call  Aggregate(nx,ny,maskmap,rnflx,shflx,lhflx,ghflx,precip,&
!                     RunOff,SatFlow,OutFlow,wetsfc,wetdp2,wetdp3,     &
!	             wetcnp,Wini,TTRnflxAvg,TTLhflxAvg,TTShflxAvg,    &
!	             TTGhflxAvg,TTPrecipAvg,TTEvapAvg,TTRunOffAvg,    &
!		     TTSatFlowAvg,TTOutflowAvg,TTWetAvg,TTWiniAvg,    &
!		     TTRnflxVar,TTLhflxVar,TTShflxVar,TTGhflxVar,     &
!		     TTPrecipVar,TTEvapVar,TTRunOffVar,TTSatFlowVar,  &
!		     TTOutflowVar,TTWetVar,TTWiniVar,TArea)
!			 
!         TPrecip = TPrecip + TTPrecipAvg
!	 TEvap   = TEvap   + TTEvapAvg
!	 TRunOff = TRunoff + TTRunOffAvg
!	 TWet    = TWet    + TTWetAvg
!	 TWini   = TWini   + TTWiniAvg

!  Propaga o escoamento superficialmente (overland and channel)
! --------------------------------------------------------------------
          if (epelopt == 2) then

             call RoutSfcTOPM (nx,ny,tint,dtroutsfc,maskmap,lnatbmap, &
	                       RunOff,TArea,NAC,IA,AC,ST,SUMAC,TL,NCH, &
			       ACH,D,ND,NR,AR,TTLhflxAvg,TTPrecipAvg, &
			       TTOutFlowAvg,TTRunOffAvg,Q)

	  else if (epelopt == 3) then  !ainda nao implementado...
	  
!            RoutSfcDHSVM1 (dtroutsfc,nx,ny,maskmap,dirmap,totaldirmap,RunOff) ! Não
!               implementado ainda...

          endif

      endif
      
! --------------------------------------------------------------------        !  Agrega dados - Aggregate data
! --------------------------------------------------------------------        !      if (epelopt .ne. 1 .and. epelopt .ne. 2) then
!         call  Aggregate(nx,ny,maskmap,rnflx,shflx,lhflx,ghflx,precip,&
!                     RunOff,SatFlow,OutFlow,wetsfc,wetdp2,wetdp3,     &
!	             wetcnp,Wini,TTRnflxAvg,TTLhflxAvg,TTShflxAvg,    &
!	             TTGhflxAvg,TTPrecipAvg,TTEvapAvg,TTRunOffAvg,    &
!		     TTSatFlowAvg,TTOutflowAvg,TTWetAvg,TTWiniAvg,    &
!		     TTRnflxVar,TTLhflxVar,TTShflxVar,TTGhflxVar,     &
!		     TTPrecipVar,TTEvapVar,TTRunOffVar,TTSatFlowVar,  &
!		     TTOutflowVar,TTWetVar,TTWiniVar,TArea)
!	 
!                                               !Unidades
!         TPrecip = TPrecip + TTPrecipAvg       ![mm/h]
!	 TEvap   = TEvap   + TTEvapAvg         !
!	 TRunOff = TRunoff + TTRunOffAvg       !
!	 TWet    = TWet    + TTWetAvg          !
!	 TWini   = TWini   + TTWiniAvg         !
!      endif
      
!-----------------------------------------------------------------------
!     Saídas do modelo - Model outputs
!-----------------------------------------------------------------------
      call EphstcOut(tint,nx,ny,maskmap,topomap,slopemap,aspectmap,   &
                    fgradmap,soilmap,sdepthmap,vegmap,ndvimap,laimap, &
		    vegfrac,tair,qvair,pres,relh,wspd,                &
                    rsi,rsr,precip,rnflx,shflx,lhflx,ghflx,Eg,Etr,Er, &
		    qvsatts,wetsfc,wetdp2,wetdp3,wetcnp,Wini,tsfc,    &
		    tdeep,Runoff,RunOffei,RunOffes,RunOffCnp,         &
                    SatFlow,OutFlow,wtdepthmap,dirmap,Q,QOBS,RnflxObs,&
		    LhflxObs,ShflxObs,GhflxObs,TTRnflxAvg,TTLhflxAvg, &
		    TTShflxAvg,TTGhflxAvg,TTPrecipAvg,TTEvapAvg,      &
		    TTRunOffAvg,TTSatFlowAvg,TTOutflowAvg,TTWetAvg,   &
		    TTWiniAvg,TTRnflxVar,TTLhflxVar,TTShflxVar,       &
		    TTGhflxVar,TTPrecipVar,TTEvapVar,TTRunOffVar,     &
		    TTSatFlowVar,TTOutflowVar,TTWetVar,TTWiniVar,TArea)
		    
!-----------------------------------------------------------------------
!     Prepara para novo passo de tempo - prepare the new time step
!-----------------------------------------------------------------------
      call NextStep(ano,mes,diam,diaj,hora,mint,dtforc,forcflopt,     &
                    tairfile,presfile,relhfile,wspdfile,rsifile,      &
		    rsrfile,precipfile,nx,ny,maskmap,topomap,         &
		    wtdepthmap,wtlevelmap,slopemap,aspectmap,fgradmap)

!-----------------------------------------------------------------------                    
!     Lê novo forçante do modelo - read a new forcing
!-----------------------------------------------------------------------
      if (tint .lt. nsteps-1) then
         call ReadForcantes(nx,ny,nx_forc,ny_forc,tairfile,           &
	                    presfile,relhfile,wspdfile,rsifile,       &
			    rsrfile,precipfile,tair,pres,relh,wspd,   &
			    rsi,rsr,precip)
      endif

   enddo   ! end of time step (dtbig) integration

!-----------------------------------------------------------------------      !     Gera estatíticas para análise - 
!     Generate statistics for analysis
!-----------------------------------------------------------------------
!   if (flx1dout == 1) then   
!      call Results (RnflxObs,LhflxObs,ShflxObs,GhflxObs)
!   endif
   
   write(6,'(a)') ' Programa encerrado com sucesso!'   
   
end program EPHSTC
