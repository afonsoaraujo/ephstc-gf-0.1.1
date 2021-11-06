!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                      MODULE INIT                     ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module init
 
!use readdata
!use calc

!   include 'globcst.inc'
!   include 'phycst.inc'
!   include 'soilcst.inc'
!   include 'calendario.inc'

!contains

!
! -------------------------------------------------------------
!  InitCalend()
!
! -------------------------------------------------------------
subroutine InitCalend(initime,ano,mes,diam,diaj,hora,mint)

   implicit none
   
   character(128)   ::   initime
   integer          ::   ano
   integer          ::   mes
   integer          ::   diam
   integer          ::   diaj
   integer          ::   hora
   integer          ::   mint
   
   integer          ::   bis
   
   include 'calendario.inc'


! -------10/02/2015
! Ajustes feitos devido a troca de compilador - ifc para gfortran
! 
   
!  Inicializa variáveis de calendario
!   ano   = inum(initime(1:4))
!   mes   = inum(initime(6:7))
!   diam  = inum(initime(9:10))
   write (initime(1:4),  '(i4)')   ano
   write (initime(6:7),  '(i2.2)') mes
   write (initime(9:10), '(i2.2)') diam

   call bissexto(ano,bis)
   call datacor( bis,mes,diam, diaj)
!   hora  = inum(initime(12:13))
!   mint   = inum(initime(15:16))

   write (initime(12:13), '(i2.2)') hora
   write (initime(15:16), '(i2.2)') mint

   write(6,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2)') ano,  &
              mes, diam, diaj, hora, mint   

!--------


end subroutine InitCalend

!
! -------------------------------------------------------------
!  InitMapGrid()
!
! -------------------------------------------------------------
subroutine InitMapGrid(nx,ny,nzsoil,x,y,dx,dy,dxy,Xorig,Yorig,        &
                       OffSetX,OffSetY,nx_forc,ny_forc,               &
                       DX_Forc,DY_Forc)
! futuramente essa subrotina deverá ser usada para definição de 
! grades em diversas projeções e para conversão de projeções
!
!#######################################################################
!
!  Declaração de variáveis
!
!#######################################################################
!
   implicit none
!
   integer          ::   nx,ny,nzsoil
   real(4)          ::   dx,dy,dxy,X,Y,Xorig,Yorig,OffSetX,OffSety
   integer          ::   nx_forc,ny_forc
   real(4)          ::   DX_forc,DY_forc
  
   integer          ::   i,j
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   namelist /ephstc_mapgrid/ nx,ny,nzsoil,DX,DY,Xorig,Yorig,OffSetX,    &
                             OffSetY,DX_Forc,DY_Forc
! 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(1x,a)') ' Inicializando a grade do modelo!'

!  Leirura das NAMELISTS
   REWIND (UNIT=5)
   READ (5, ephstc_mapgrid)

   do j = 1, ny
      Y= Yorig + j
      do i = 1, nx
         X = Xorig + i
      enddo
   enddo
  
   dxy = DX*DX + DY*DY
   dxy = SQRT(dxy)
   
!   write(6,'(2i6,f10.2)') nx, ny, dxy
   
end subroutine InitMapGrid

!
! -------------------------------------------------------------
!  InitForcantes()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...), bem 
!     como outras características do solo.
!     
!     Verificar a necessidade futura de passar TopoMap para 
!     permitir o uso de Mask na leitura dos arquivos???
! -------------------------------------------------------------
subroutine ReadForcantes (nx,ny,nx_forc,ny_forc,tairfile,presfile,    &
                          relhfile,wspdfile,rsifile,rsrfile,          &
                          precipfile,tair,pres,&
                          relh,wspd,rsi,rsr,precip)         

   implicit none
   include 'globcst.inc'
   
   integer          ::   nx,ny,nx_forc,ny_forc
   character(128)   ::   tairfile,presfile,relhfile,wspdfile,rsifile, &
                         rsrfile,precipfile
   real             ::   tair(nx,ny),   tairforc  (nx_forc,ny_forc)
   real             ::   pres(nx,ny),   presforc  (nx_forc,ny_forc)
   real             ::   relh(nx,ny),   relhforc  (nx_forc,ny_forc)
   real             ::   wspd(nx,ny),   wspdforc  (nx_forc,ny_forc)
   real             ::   rsi(nx,ny),    rsiforc   (nx_forc,ny_forc)
   real             ::   rsr(nx,ny),    rsrforc   (nx_forc,ny_forc)
   real             ::   precip(nx,ny), precipforc(nx_forc,ny_forc)
   real             ::   precipag, precipag_ac
   real             ::   precipfag, precipfag_ac
   
   integer          ::   i,j,k,ii,jj

   logical          ::   firstcall
   SAVE             ::   firstcall
   DATA                  firstcall /.true./
   character(128)   ::   precipagfl
!
!#######################################################################
!
!     Namelists:
!
!#######################################################################
!
!   namelist /ephstc_forcfl/  tairfile,presfile,relhfile,wspdfile,        &
!                            rsifile,rsrfile,precipfile      
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   if (firstcall) then
      firstcall = .false.
      write(6,*) 'Primeira chamada à ReadForc'
   endif
   

! Leitura dos Forçantes do Modelo
   write(6,'(1x,a,i6)') ' Lendo forçantes do modelo...!', ano

!   write(6,'(1x,a)') ' Lendo Temperatura do Ar na superfície...!'
!   write(6,'(1x,a)') tairfile
   call readbinf(1,1,tairfile,NX_Forc,NY_Forc,tairforc)
!   write(6,'(1x,a)') ' Temperatura do Ar na superfície lida!'   
   
!   write(6,'(1x,a)') ' Lendo Pressão do Ar na superfície...!'
!   write(6,'(1x,a)') AirPressfl
   call readbinf(1,2,presfile,NX_Forc,NY_Forc,presforc)
!   write(6,'(1x,a)') ' Temperatura do Ar na superfície lida!'   

!   write(6,'(1x,a)') ' Lendo Umidade Relativa do Ar na superfície...!'
!   write(6,'(1x,a)') Relhfl
   call readbinf(1,3,relhfile,NX_Forc,NY_Forc,relhforc)

!   write(6,'(1x,a)') ' Lendo Velocidade do vento na superfície...!'
!   write(6,'(1x,a)') Wspdfl
   call readbinf(1,4,wspdfile,NX_Forc,NY_Forc,wspdforc)

!   write(6,'(1x,a)') ' Lendo Radiação solar insidente...!'
!   write(6,'(1x,a)') Rsifl
   call readbinf(1,5,rsifile,NX_Forc,NY_Forc,rsiforc)

!   write(6,'(1x,a)') ' Lendo Radiação solar refletida...!'
!   write(6,'(1x,a)') Rsrfl
   call readbinf(1,6,rsrfile,NX_Forc,NY_Forc,rsrforc)

!   write(6,'(1x,a)') ' Lendo Precipitação...!'
!   write(6,'(1x,a)') Precipfl
   call readbinf(1,7,precipfile,NX_Forc,NY_Forc,precipforc)

! "Interpolação" para o domínio do modelo
! No momento o processo de "interpolação consiste em simplesmente
! considerar o campo homogêneo.....

!   write(6,'(4i4)') nx_forc,ny_forc,NX_forc,NY_forc
!   do j=1,NY_Forc
!   do i=1,NX_Forc
!      write(6,'(2i4,3f16.4)') i,j,tair(i,j),rsiforc(i,j),rsrforc(i,j)
!   enddo
!   enddo

   precipag  = 0.0
   precipfag = 0.0

!   do jj = 1, ny_forc
!   do ii = 1, nx_forc
!      precipfag = precipfag + precipforc(ii,jj)
!   enddo
!   enddo

   jj = 1
   do j = 1, ny
      ii = 1
      do i = 1, nx
         tair(i,j)   = tairforc(ii,jj)         ! Conversão para Kelvin
         pres(i,j)   = presforc(ii,jj)         ! Conversão de mb para Pa
         relh(i,j)   = relhforc(ii,jj) / 100
         wspd(i,j)   = wspdforc(ii,jj) 
         rsi(i,j)    = rsiforc(ii,jj)
         rsr(i,j)    = rsrforc(ii,jj)
         precip(i,j) = precipforc(ii,jj)       ! Coversao mm/h ==> kg/s.m2
	 
!	 write(6,'(2f12.4,4i4)') tair(i,j),tairforc(ii,jj),i,j,ii, jj
	 
!	 write(6,'(2f12.4,4i4)') rsi(i,j),rsiforc(ii,jj),i,j,ii, jj

         if (i > ii*(int(nx/nx_forc)+1)) then
            ii = ii + 1
            if (ii > nx) then
               ii = nx
            endif
         endif
         
      enddo
      
      if (j > jj*(int(ny/ny_forc))) then
         jj = jj + 1
         if (jj > ny) then
            jj = ny
         endif
      endif
   enddo
    
   return

end subroutine ReadForcantes

!
! -------------------------------------------------------------
!  InitModelState()
!
! esta rotina deverá permitir a inicialização e reinicialização
! dos estados iniciais do modelo.Inicialmente serão inicializa-
! das as seguintes variáveis:
! 1. WaterTableDepth - será inicializada a partir da Rotina 
!                      CalcWtDepth alimentada com as umidades
!                      iniciais do solo Ws,W2 e W3.
! -------------------------------------------------------------
subroutine InitModelState (nx,ny,nzsoil,tair,pres,relh,tsfc,tdeep,    &
                           wetsfc,wetdp2,wetdp3,wetcanp,snowdepth,    &
                           qvsatta,qvsatts,qvsfc,rhoa,Wini)
   implicit none
   include 'soilcst.inc'

   integer          ::   nx,ny,nzsoil
   real             ::   tair(nx,ny)
   real             ::   pres(nx,ny)
   real             ::   relh(nx,ny)
   real             ::   tsfc(nx,ny)
   real             ::   tdeep(nx,ny)
   real             ::   wetsfc(nx,ny)
   real             ::   wetdp2(nx,ny)
   real             ::   wetdp3(nx,ny)
   real             ::   wetcanp(nx,ny)
   real             ::   snowdepth(nx,ny)
   real             ::   qvsatta(nx,ny)
   real             ::   qvsatts(nx,ny)
   real             ::   qvsfc(nx,ny)
   real             ::   rhoa(nx,ny)
   real             ::   Wini(nx,ny)

   integer          ::   i,j
   real             ::   aw,bw
!
!#######################################################################
!
!  Include:
!
!#######################################################################
!
   include 'globcst.inc'
   include 'phycst.inc'
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   namelist /ephstc_initmodelstate/ smopt,smxyfile,smzfile,tsfcini,    &
                                   tdp2ini,wsini,w2ini,w3ini,wcini,   &
				   sdini
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   write(6,'(a)') '   '
   write(6,'(1x,a)') 'Definindo o estado inicial do modelo/Initializing Model States...!'
!---------------------------------------------------------------------
!  Lê namelsits
!---------------------------------------------------------------------      
   REWIND(UNIT=5)
   READ (5, ephstc_initmodelstate)

!---------------------------------------------------------------------
!  Inicializa variáveis para integração....
!---------------------------------------------------------------------
   do j=1,ny
   do i=1,nx
      tsfc(i,j)      = tsfcini              !Tg
      tdeep(i,j)     = tdp2ini              !T2 
      wetsfc(i,j)    = wsini                !Wg - devem ser substituidos por valores
      wetdp2(i,j)    = w2ini                !W2
      wetdp3(i,j)    = w3ini                !W3
      if (epstopt == 0) then
         wetdp3(i,j)    = 0.0
      endif
      wetcanp(i,j)   = wcini                !Wr
      snowdepth(i,j) = sdini                ! snowdepth
      Wini(i,j)      = (wetsfc(i,j)*d1 + wetdp2(i,j)*(d2-d1))/2 +     &
                        wetdp3(i,j)*(d3-d2) + wetcanp(i,j)*dr
   enddo
   enddo
	 
!---------------------------------------------------------------------
!  Cálculo da umidade específica de saturação para AirTemp
!---------------------------------------------------------------------
   do j=1, NY
   do i=1, NX
   
      if ( tair(i,j) .ge. 273.15 ) then
         aw = 17.27
         bw = 35.5
      else
         aw = 21.875
         bw = 7.5
      end if
         
      qvsatta(i,j) = (0.622/pres(i,j)) * 611 * exp(aw * (tair(i,j) -  &
         273.15) / (tair(i,j) - bw))
      
      if (tair(i,j) .ge. 273.15) then
         aw = 17.27
         bw = 35.5
      else
         aw = 21.875
         bw = 7.5
      end if

      qvsatts(i,j) = (0.622 / pres(i,j)) * 611 * exp(aw * (tair(i,j) -   &
         273.15) /  (tair(i,j) - bw))
      
!      Umidade específica do ar
       qvsfc(i,j) = relh(i,j) * qvsatta(i,j)
	 
!      write(6,'(3f15.8)') qvsatta(i,j), qvsatts(i,j),qva(i,j)

      rhoa(i,j) = roa
   enddo
   enddo      

   return

end subroutine InitModelState
!
! -------------------------------------------------------------
!  InitTerrainMaps()
! -------------------------------------------------------------
subroutine InitTerrainMaps(nx,ny,nzsoil,maskmap,topomap,aspectmap,    &
                           slopemap,lnatbmap,fgradmap,dirmap,totaldirmap,&
                           hidromap,soilmap,sdepthmap,vegmap,vegfrac, &
                           ndvimap,laimap,wetsfc,wetdp2,wetdp3,       &
                           wtdepthmap,wtlevelmap)

   implicit none
   include 'grid.inc'
   
   integer          ::   nx,ny,nzsoil

   integer          ::   maskmap(nx,ny)
   real             ::   topomap(nx,ny)
   real             ::   aspectmap(nx,ny)
   real             ::   slopemap(nx,ny)
   real             ::   lnatbmap(nx,ny)
   real             ::   fgradmap(nx,ny)
   integer          ::   dirmap(nx,ny,NDIRS) 
   integer          ::   totaldirmap(nx,ny)             
   integer          ::   hidromap(nx,ny)
   integer          ::   soilmap(nx,ny)
   real             ::   sdepthmap(nx,ny)
   integer          ::   vegmap(nx,ny)
   real             ::   vegfrac(nx,ny)
   real             ::   ndvimap(nx,ny)
   real             ::   laimap(nx,ny)
   real             ::   wetsfc(nx,ny)
   real             ::   wetdp2(nx,ny)
   real             ::   wetdp3(nx,ny)
   real             ::   wtdepthmap(nx,ny)
   real             ::   wtlevelmap(nx,ny)
!
!#######################################################################
!  Includes
!#######################################################################
!   
   include 'soilcst.inc'
   include 'globcst.inc'
!
!#######################################################################
!  Namelists:
!#######################################################################
!
   namelist /ephstc_sfcfl/   maskfile,topofile,aspectfile,slopefile,   &
	                    lnatbfile,hidrofile,soilfile,sdepthfile,  &
			    vegfile,ndvifile
   namelist /ephstc_initopt/ slopeopt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executavel
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(a)') '   '
   write(6,'(1x,a)') 'Inicializando mapas do terreno - Initializing terrain maps...!'

   REWIND(UNIT=5)
   READ (5, ephstc_sfcfl)
   READ (5, ephstc_initopt)
   
   call InitMaskMap (nx,ny,maskfile,maskmap)

   call InitSoilMaps (nx,ny,nzsoil,soilfile,sdepthfile,soilmap,       &
                      sdepthmap)

   call InitVegMaps (nx,ny,vegfile,ndvifile,vegmap,vegfrac,ndvimap,   &
                     laimap)
                     
   call InitTopoMap (nx,ny,maskmap,topomap)                     
						  
   call InitWtMaps (nx,ny,nzsoil,sdepthmap,topomap,vegmap,soilmap,    &
                    wetsfc,wetdp2,wetdp3,wtlevelmap,wtdepthmap)
		    
   call InitSlopeAspectMaps (nx,ny,maskmap,wtlevelmap,wtdepthmap,     &
                             topomap,aspectmap,slopemap,fgradmap,     &
                             dirmap,totaldirmap)
			     
   call InitLnaTbMap (nx,ny,maskmap,lnatbmap)			     
                             
!   call InitHidroMap (nx,ny,maskmap,hidromap)                             
                      
end subroutine InitTerrainMaps

!
! -------------------------------------------------------------
!  InitMaskMap()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...).
! -------------------------------------------------------------
subroutine InitMaskMap(nx,ny,maskfile,maskmap)

   integer          ::   nx,ny
   character(128)   ::   maskfile
   integer          ::   maskmap(nx,ny)
   
   integer          ::   i,j             ! Counters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(a)') '   Lendo Máscara da bacia...'
   write(6,'(1x,a)') maskfile

!  Abre arquivo com máscara
   open(1,file=maskfile,status='unknown')
   
!  Le arquivo 
   do j = 1, NY
      read(1,*)(maskmap(i,j),     i=1,NX)
   enddo
   
   do j = 1, NY
   do i = 1, NX
      if (maskmap(i,j) == -9999) then
         maskmap(i,j) = 0
      endif
   enddo
   enddo

   write(6,'(a)') ' OK, Máscara da bacia lida!'
   close(1)

   return

end subroutine InitMaskMap

!
! -------------------------------------------------------------
!  InitSoilMaps()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...), bem 
!     como outras características do solo.
!     
!     Verificar a necessidade futura de passar TopoMap para 
!     permitir o uso de Mask na leitura dos arquivos???
! -------------------------------------------------------------
subroutine InitSoilMaps (nx,ny,nzsoil,soilfile,sdepthfile,soilmap,    &
                         sdepthmap)
   implicit none		       

   integer          ::   nx,ny,nzsoil
   character(128)   ::   soilfile,sdepthfile
   integer          ::   soilmap(nx,ny)
   real             ::   sdepthmap(nx,ny)

   integer          ::   soilmapopt 
   integer          ::   i,j,k                     ! Counters

!  debug auxiliary variables
   integer          :: STAT_ALLOC   
   
   write(6,'(a)') 'Lendo Mapas de solos  ' 
   write(6,'(a)') soilfile
   write(6,'(a)'),sdepthfile

!  Open Soil datasets
   open(1,file=soilfile, status='unknown')
   open(2,file=sdepthfile,status='unknown')
   
!  Read Soil Maps
   do j = 1, NY
      read(1,*)(soilmap(i,j), i=1,NX)
   enddo
   
   soilmapopt = 2
   if (soilmapopt == 2) then
!  Converte os códigos de textura de solo de GRIJALVA para o EPHSTC[ARPS]

   do j = 1, NY
   do i = 1, NX
      
      if ( soilmap(i,j) == 1 ) then 
   
         soilmap(i,j) = 6
   
      else if ( soilmap(i,j) == 2 ) then
   
         soilmap(i,j) = 1
   
      else if ( soilmap(i,j) == 3 ) then
   
         soilmap(i,j) = 11
	 
      else if ( soilmap(i,j) == 0 ) then	 

         soilmap(i,j) = 13
	 
      endif

   enddo
   enddo
   
   endif
   
   
   do j = 1, NY
      read(2,*)(sdepthmap(i,j), i=1,NX)
   enddo
   
   write(6,'(a)') ' OK, Mapas de solos lidos !'
   close(1)
   close(2)

   return

end subroutine InitSoilMaps

! -------------------------------------------------------------
!  InitVegMaps()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...), bem 
!     como outras características do solo.
!     
!     Verificar a necessidade futura de passar TopoMap para 
!     permitir o uso de Mask na leitura dos arquivos???
! -------------------------------------------------------------
subroutine InitVegMaps(nx,ny,vegfile,ndvifile,vegmap,vegfrac,ndvimap, &
                       laimap)

   implicit none

   integer          ::   nx,ny
   character(128)   ::   vegfile, ndvifile
   integer          ::   vegmap(nx,ny)
   real             ::   vegfrac(nx,ny)
   real             ::   ndvimap(nx,ny)
   real             ::   laimap(nx,ny)
   
   integer          ::   vegmapopt  ! opcao de tabela de classificacao.

   integer          ::   i,j                     ! Counters
   real,external    ::   CalcLaiNDVI
   
!  debug auxiliary variables
   integer          :: STAT_ALLOC   
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!   
   include 'soilcst.inc'   
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executavel
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   write(6,'(a)') 'Lendo Mapas de vegetação...!  ' 
   write(6,'(a)') vegfile
   write(6,'(a)') ndvifile

!  Open Vegetation datasets
   open(1,file=vegfile, status='unknown')
   open(2,file=ndvifile,status='unknown')
   
!  Read Veg Maps
   do j = 1, NY
      read(1,*)(vegmap(i,j),  i=1,NX)
      read(2,*)(ndvimap(i,j), i=1,NX)
   enddo

   vegmapopt = 2 
   if (vegmapopt == 1) then
!  Converte os códigos de vegetação do SGP97 para o EPHSTC[ARPS]

   do j = 1, NY
   do i = 1, NX
      if ( vegmap(i,j) == 2 ) then
      
         vegmap(i,j) = 13
   
      else if ( vegmap(i,j) == 6 ) then
   
         vegmap(i,j) = 3
   
      else if ( vegmap(i,j) == 13 ) then 
   
         vegmap(i,j) = 4   
   
      else if ( vegmap(i,j) == 0 .OR.  vegmap(i,j) == 7 .OR.     &
      vegmap(i,j) == 8 ) then
   
         vegmap(i,j) = 5
   
      else if ( vegmap(i,j) == 1 .OR. vegmap(i,j) == 3  .OR.     &
      vegmap(i,j) == 4  .OR. vegmap(i,j) == 5  .OR.              &
      vegmap(i,j) == 10 .OR. vegmap(i,j) == 11 .OR.              &
      vegmap(i,j) == 12 ) then
   
         vegmap(i,j) = 10
   
      else if ( vegmap(i,j) == 9 .OR. vegmap(i,j) == -9999 ) then
      
         vegmap(i,j) = 14
	 
      else if ( vegmap(i,j) == 13 ) then
      
      vegmap(i,j) = 12
      
      endif
      
   enddo
   enddo
      
   else if (vegmapopt == 2) then
!  Converte os códigos de vegetação de GRIJALVA para o EPHSTC[ARPS]

   do j = 1, NY
   do i = 1, NX
      if ( vegmap(i,j) == 16 ) then

         vegmap(i,j) = 1
   
      else if ( vegmap(i,j) == 23 ) then
   
         vegmap(i,j) = 3
   
      else if ( vegmap(i,j) == 7 .OR.  vegmap(i,j) == 27 .OR.     &
      vegmap(i,j) == 28 .OR. vegmap(i,j) == 32) then
   
         vegmap(i,j) = 4
      
      else if ( vegmap(i,j) == 8 ) then
   
         vegmap(i,j) = 5
   
      else if ( vegmap(i,j) == 20 .OR.  vegmap(i,j) == 18 .OR.     &
      vegmap(i,j) == 34) then
   
         vegmap(i,j) = 6
	 
      else if ( vegmap(i,j) == 19 .OR. vegmap(i,j) == 21 .OR.      &
                vegmap(i,j) == 26 .OR. vegmap(i,j) == 29 .OR.      &
		vegmap(i,j) == 30 .OR. vegmap(i,j) == 31 .OR.      &
		vegmap(i,j) == 36) then
		
         vegmap(i,j) = 7

      else if ( vegmap(i,j) == 22 ) then
   
         vegmap(i,j) = 8

      else if ( vegmap(i,j) == 42 ) then
   
         vegmap(i,j) = 9

      else if ( vegmap(i,j) ==  4 .OR. vegmap(i,j) ==  9 .OR.      &
                vegmap(i,j) == 10 .OR. vegmap(i,j) == 11 .OR.      &
		vegmap(i,j) == 12 .OR. vegmap(i,j) == 13 .OR.      &
		vegmap(i,j) == 14 .OR. vegmap(i,j) == 15 .OR.      &
		vegmap(i,j) == 17 .OR. vegmap(i,j) == 24 .OR.      &
		vegmap(i,j) == 25 .OR. vegmap(i,j) == 33 .OR.      &
		vegmap(i,j) == 35 .OR. vegmap(i,j) == 37 .OR.      &
		vegmap(i,j) == 38 .OR. vegmap(i,j) == 39 .OR.      &
		vegmap(i,j) == 40 .OR. vegmap(i,j) == 41) then
		
         vegmap(i,j) = 10
	 
      else if ( vegmap(i,j) == 2 .OR. vegmap(i,j) == 3 ) then
   
         vegmap(i,j) = 11
	 
      else if ( vegmap(i,j) ==  23) then
   
         vegmap(i,j) = 12 
      
      else if ( vegmap(i,j) == 6 ) then
   
         vegmap(i,j) = 13
	 
      else if ( vegmap(i,j) == 1 .OR. vegmap(i,j) == 5 ) then
   
         vegmap(i,j) = 14
	 
      endif
   enddo
   enddo
   
   endif  ! 

!  O NDVI do arquivo está normalizado: N_NDVI=100(NDVI+1)
   do j = 1, NY
   do i = 1, NX
      vegfrac(i,j) = vegf(vegmap(i,j))
      ndvimap(i,j) = (ndvimap(i,j)/100) - 1
      laimap(i,j)= CalcLaiNDVI(vegmap(i,j), ndvimap(i,j))
   enddo
   enddo
   
   write(6,'(a)') ' OK, Mapas de vegetacao lidos !'
   close(1)
   close(2)

   return

end subroutine InitVegMaps

!
! --------------------------------------------------------------------
!  InitTopoMap()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...).
! --------------------------------------------------------------------
subroutine InitTopoMap (nx,ny,maskmap,topomap)
   implicit none
   include 'grid.inc'
   include 'globcst.inc'
    
   integer          ::   nx,ny
   integer          ::   maskmap(nx,ny)
   real             ::   topomap(nx,ny)
   
   integer          ::   i,j,k             ! Counters
   
!  debug auxiliary variables
   integer   :: STAT_ALLOC
   real      :: min,max
   integer   :: minx,miny,maxx,maxy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executavel
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(a)') ' Lendo mapa de topografia...'
   write(6,'(a)') topofile

!  Abre arquivos de MDT - Open DEM datasets
! --------------------------------------------------------------------
   open(1,file=topofile,status='unknown')
   
!  Le arquivos de topografia - Read topography Map
! --------------------------------------------------------------------
   do j = 1, ny
      read(1,*)(topomap(i,j) , i=1,nx)
   enddo
   
   write(6,'(a)') ' OK, Mapa de topografia lido!'
   close(1)

   return

end subroutine InitTopoMap

! --------------------------------------------------------------------
!  InitWtMaps()
!  Inicializa a profundidade do lençol e o nível de água no solo, como 
!  uma fração da profundidade total do solo.
! --------------------------------------------------------------------
subroutine InitWtMaps (nx,ny,nzsoil,sdepthmap,topomap,vegmap,soilmap, &
                       wetsfc,wetdp2,wetdp3,wtlevelmap,wtdepthmap)
                      
   implicit none
   include 'grid.inc'
   include 'soilcst.inc'
   include 'globcst.inc'
    
   integer            ::   nx,ny,nzsoil
   real               ::   sdepthmap(nx,ny)
   real               ::   topomap(nx,ny)
   integer            ::   vegmap(nx,ny)
   integer            ::   soilmap(nx,ny)
   real               ::   wetsfc(nx,ny)
   real               ::   wetdp2(nx,ny)
   real               ::   wetdp3(nx,ny)
   real               ::   wtlevelmap(nx,ny)
   real               ::   wtdepthmap(nx,ny)
   
   integer            ::   i,j,k,smp
!-----------------------------------------------------------------------
! Variáveis para a atualização da linha d'água
!-----------------------------------------------------------------------
   real                  ::   NRootLayers
   real,   allocatable   ::   LPorosity(:)
   real,   allocatable   ::   LFCap(:)
   real,   allocatable   ::   RootDepth(:)
   real,   allocatable   ::   Moisture(:)
   real(4)               ::   CalcWTableDepth
!   real(4)               ::   CalcRootDepth
   integer               ::   istatus   
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!---------------------------------------------------------------------    
! Aloca memória se necessário - Alocate memory if necessary
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
! --------------------------------------------------------------------     
   write(6,'(a)') ' Lendo mapa de profundidade da rocha...'
   write(6,'(a)') sdepthfile    
! --------------------------------------------------------------------  
!  Abre arquivo de profundidade média da rocha - Open rock depth file
! --------------------------------------------------------------------
   open(1,file=sdepthfile,status='unknown')

! --------------------------------------------------------------------
!  Le arquivo de profundidade média da rocha - Read rock depth map
! --------------------------------------------------------------------   
   do j = 1, ny
      read(1,*)(sdepthmap(i,j) , i=1,nx)
   enddo   

   write(6,'(a)') ' OK, Mapa de profundidade da rocha lido!'
   close(1)

!---------------------------------------------------------------------
!  Define o estado inicial da linha d'água em função das umidades 
!  iniciais atribuídas a Ws,W2 e W3. Futuramente, essa inicialização 
!  deve ser feita por algoritmo de assimilação de dados de umidade 
!  do solo mais apropriados.
!---------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx
   
      do k =1, nzsoil
         LPorosity(k) = porosity(soilmap(i,j))      ! constante ao longo de z
         LFCap(k) = wfc(soilmap(i,j))               ! constante ao longo de z
      enddo

      call CalcRootDepth(nzsoil,vegmap(i,j),NRootlayers,RootDepth,    &
                         sdepthmap(i,j),i,j)
                         
      Moisture(1) = wetsfc(i,j)
      Moisture(2) = wetdp2(i,j)
      Moisture(3) = wetdp3(i,j)

      wtdepthmap(i,j) = CalcWTableDepth(NRootLayers,sdepthmap(i,j),   &
                                        RootDepth,LFCap,LPorosity,    &
                                        Moisture,i,j)
           
! --------------------------------------------------------------------   
!  Define o nível da linha de água
! --------------------------------------------------------------------   
         wtlevelmap(i,j) = topomap(i,j) - wtdepthmap(i,j)
	 
   enddo
   enddo
   
   return

end subroutine InitWtMaps

!
! --------------------------------------------------------------------
!  InitSlopeAspectMaps()
!      esta rotina deverá ser utilizada futuramente para a
!     implementacao de leituras de diferentes tipos de 
!     arquivos de dados (ex.: BIN, ASC, NetCDF,...).
! --------------------------------------------------------------------
subroutine InitSlopeAspectMaps(nx,ny,maskmap,wtlevelmap,wtdepthmap,   &
                               topomap,aspectmap,slopemap,fgradmap,   &
                               dirmap,totaldirmap)
   implicit none
   include 'grid.inc'
   include 'globcst.inc'
    
   integer          ::   nx,ny
   integer          ::   maskmap(nx,ny)
   real             ::   wtlevelmap(nx,ny)
   real             ::   wtdepthmap(nx,ny)
   real             ::   topomap(nx,ny)
   real             ::   aspectmap(nx,ny)
   real             ::   slopemap(nx,ny)
   real             ::   fgradmap(nx,ny)
   integer          ::   dirmap(nx,ny,NDIRS)
   integer          ::   totaldirmap(nx,ny)

   integer          ::   xn, yn, n
   real             ::   neighbor_elev(NDIRS)
   integer          ::   XIndex, YIndex, valid_cell
   
   integer          ::   i,j,k             ! Counters
   
!  debug auxiliary variables
   integer   :: STAT_ALLOC
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do codigo executavel
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(a)') '  '
   write(6,'(a)') ' Inicializando slope and aspect maps! '
   write(6,'(1x,3a,i4)') slopefile,aspectfile,' slopeopt=',slopeopt
! --------------------------------------------------------------------
! Calculate slope, aspect, magnitude of subsurface flow gradient, and 
! fraction of flow flowing in each direction based on the land surface 
! slope.
! --------------------------------------------------------------------
! as opções devem ser futuramente movidas para NAMELISTS

   if (slopeopt == 0) then
!  Lê arquivos de slope e aspect pre-processados - Read Slope and 
!  Aspect preprocessed files
! --------------------------------------------------------------------  
      open(3,file=slopefile, status='unknown')
      open(4,file=aspectfile,status='unknown')
   
      do j = 1, NY
         read(3,*)(slopemap(i,j), i=1,NX)
         read(4,*)(aspectmap(i,j), i=1,NX)
      enddo

      write(6,'(a)') ' OK, mapas de declividade e aspecto lidos!'      
      close(3)
      close(4)
      
      do j = 1, ny
      do i = 1, nx
         if (maskmap(i,j) .eq. 1)  then
            do n = 1, NDIRS
               xn = i + XIndex(n)
               yn = j + YIndex(n)

               if (valid_cell(nx,ny,xn,yn) .eq. 1) then
                  if (maskmap(xn,yn) .eq. 1) then
                     neighbor_elev(n) = topomap(xn,yn)
!                    write(6,'(1x,3i4,f8.2)') xn,yn,n,neighbor_elev(n)
                  else
                     neighbor_elev(n) = 0   ! OUTSIDEBASIN
!                    write(6,'(1x,3i4,f8.2)') xn,yn,n,neighbor_elev(n)
                  endif
               else
                  neighbor_elev(n) = 0   ! OUTSIDEBASIN
               endif
            enddo

            call flow_fractions(slopemap(i,j),aspectmap(i,j),         &
                                fgradmap(i,j),wtdepthmap(i,j),        &
				neighbor_elev,dirmap(i,j,:),          &
				totaldirmap(i,j),i,j)
         endif
      enddo
      enddo            

   else if (slopeopt == 1) then
!  Calcula o "slope", "aspect" e "flow gradiets" em função da 
!  profundidade da linha d'água.
! --------------------------------------------------------------------  
      call HeadSlopeAspect(nx,ny,maskmap,topomap,wtlevelmap,slopemap, &
                           aspectmap,fgradmap,dirmap,totaldirmap)

      write(6,'(a)') ' OK, mapas de "declividades", "aspecto" e "gradientes" calculados - HeadSlopeAspect!'

   else if (slopeopt == 2) then
!  Calcula o "slope", "aspect" e "flow gradiets" em função da 
!  topografia do terreno.
! --------------------------------------------------------------------      
      call ElevationSlopeAspect(nx,ny,maskmap,topomap,wtdepthmap,     &
                                slopemap,aspectmap,fgradmap,dirmap,   &
                                totaldirmap)
      write(6,'(a)') ' OK, mapas de "declividades", "aspecto" e "gradientes" calculados - ElevationSlopeAspect!'
                                
   endif

   return

end subroutine InitSlopeAspectMaps

!
! --------------------------------------------------------------------
!  InitLnaTbMap() - Inicializa o mapa do indice topografico
!                   IT = Ln(a/TgB)
! --------------------------------------------------------------------
subroutine InitLnaTbMap (nx,ny,maskmap,lnatbmap)			     

   implicit none
   include 'grid.inc'
   include 'globcst.inc'
    
   integer          ::   nx,ny
   integer          ::   maskmap(nx,ny)
   real             ::   lnatbmap(nx,ny)   
   
   integer          ::   i,j,k             ! Counters   

   open(1,file=lnatbfile, status='unknown')
   
   do j = 1, NY
      read(1,*)(lnatbmap(i,j), i=1,NX)
   enddo
   
   write(6,'(a)') '    '
   write(6,'(a)') ' OK, mapa do indice topografico lido!'      
   close(1)

   return
   
end subroutine InitLnaTbMap   

!
! --------------------------------------------------------------------
!  InitHidroMap() - Inicializa o mapa de hidrografia da seguinte forma:
!                   células contendo hidrografia = 255
!                   células sem hidrografia      =   1  
! --------------------------------------------------------------------
subroutine InitHidroMap (nx,ny,maskmap,hidromap)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'
    
   integer          ::   nx,ny
   integer          ::   maskmap(nx,ny)
   integer          ::   hidromap(nx,ny)   
   
   integer          ::   i,j,k             ! Counters   

   open(1,file=hidrofile, status='unknown')
   
   do j = 1, NY
      read(1,*)(hidromap(i,j), i=1,NX)
   enddo
   
   write(6,'(a)') ' OK, mapa de hidrografia lido!'      
   close(1)   
   
   return
   
end subroutine InitHidroMap


!
! -------------------------------------------------------------
!  InitSvatpar()
!
! esta rotina inicializa os parâmetros utilizados nas diferentes
! rotinas de svats
!
! -------------------------------------------------------------
subroutine InitSvatpar (nx,ny,cdha,cdqa)

   implicit none
   
   integer   ::   nx,ny
   real      ::   cdha(nx,ny)
   real      ::   cdqa(nx,ny)
   
   integer   ::   i,j
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   include 'globcst.inc'
   include 'soilcst.inc'
!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   namelist /ephstc_svatpar/  Wl,Cdhlnd,Cdhwtr,Cdqlnd,Cdqwtr,C4lref,   &
                             C4b,Bvic
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(1x,a)') ' Inicializando parâmetros do Svat'
   REWIND(UNIT=5)
   READ (5, ephstc_svatpar)
!
!#######################################################################
!
!  Inicializando coeficientes de arrasto:
!
!#######################################################################
!
   do j=1, ny
   do i=1, nx
      cdha(i,j) = Cdhlnd    ! No momento não vamos diferenciar terra e água
      cdqa(i,j) = Cdqlnd
   enddo
   enddo

end subroutine InitSvatpar

!
! -------------------------------------------------------------
!  InitTopmData()
!
! esta rotina inicializa os dados e parâmetros necessarios para
! o TOPMODEL
!
! -------------------------------------------------------------
subroutine InitTopmData (AREA,NAC,IA,AC,ST,SUMAC,TL,NCH,ACH,D,ND,NR,  &
                         AR)

   implicit none   
   include 'globcst.inc'
   include 'topmodel.inc'

   
!  variáveis principais de conecao com o EPHSTC
!---------------------------------------------------------------------   
   real               ::   AREA
   integer            ::   NAC,IA,NSC
   real               ::   AC(NACMX)
   real               ::   ST(NACMX)
   real               ::   SUMAC
   real               ::   TL
   integer            ::   NCH
   real               ::   ACH(NACMX)
   real               ::   D(NACMX)
   real               ::   AR(NCHMX)
   integer            ::   ND,NR
   integer            ::   IMAP,IOUT
   character(128)     ::   SUBCAT
   
   real               ::   Q(nsteps)
   real               ::   QOF
   real               ::   QOUT
   real               ::   SUMQ
   
   real               ::   tarea

!  variáveis locais
!---------------------------------------------------------------------   
   integer            ::   i,j,k        ! counters
   
!
!#######################################################################
!  Listas de nomes - Namelists:
!#######################################################################
!
   namelist /ephstc_topmodel/  CHVDT,RVDT,Q0,itfile
!
!#######################################################################
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!#######################################################################
!
   write(6,'(a)')    '/-----------------------------------------------/   '
   write(6,'(1x,a)') '/---- Entrando na rotina InitTopmData-----------/'
!
!#######################################################################
!  Leitura das listas de nomes - Namelists reading:
!#######################################################################
!
!  Lê NAMELISTS 
   rewind(UNIT=5)
   read (5, ephstc_topmodel)
   
   open(110,file=itfile,form='formatted',status='unknown') 
      
!  READ IN SUBCATCHMENT TOPOGRAPHIC DATA
   read(110,*) NSC,IMAP,IOUT
      
   READ(110,"(A)")subcat
   READ(110,*) NAC,AREA

!  NAC IS NUMBER OF A/TANB ORDINATES
!  AREA IS SUBCATCHMENT AREA AS PROPORTION OF TOTAL CATCHMENT 
   READ(110,*) (AC(J),ST(J),J=1,NAC)

!  AC IS DISTRIBUTION OF AREA WITH LN(A/TANB)
!  ST IS LN(A/TANB) VALUE
   tarea = ac(1)
       
   do j=2,nac
      tarea = tarea + ac(j)
   enddo

!  CALCULATE AREAL INTEGRAL OF LN(A/TANB)
!  NB.  a/tanB values should be ordered from high to low with ST(1)
!  as an upper limit such that AC(1) should be zero, with AC(2) representing
!  the area between ST(1) and ST(2)
   TL=0.
   AC(1)=AC(1)/tarea
   SUMAC=AC(1)
      
   DO J=2,NAC
      AC(J)=AC(J)/tarea
      SUMAC=SUMAC+AC(J)
      TL=TL+AC(J)*(ST(J)+ST(J-1))/2
   enddo
   AC(NAC+1)=0.       
   
!  READ CHANNEL NETWORK DATA
   READ(110,*)NCH
   READ(110,*)(ACH(J),D(J),J=1,NCH)
!   do J=1,NCH
!      write(6,'(a,2f18.6)') 'ach= ',ACH(J),D(J)
!   enddo
!  ACH IS CUMULATIVE DISTRIBUTION OF AREA WITH DISTANCE D
!  FROM OUTLET.  D(1) is distance from subcatchment outlet
!  ACH(1) = 0.
   ND=0
   NR=0
!   write(6,'(a,3i8,f18.6)') 'NCH= ', NCH,ND,NR,AREA
   call DA2TLH (NCH,ACH,D,AREA,ND,NR,AR)
   write(6,'(a,3i8,f18.6)') 'NCH1= ', NCH,ND,NR,AREA

   write(6,'(a)') '    '
   write(6,'(a)') ' OK, Dados para o TOPMODEL iniciados!'      
   
   close(110)   
   
   return

end subroutine InitTopmData   
   
      
!
!
! -------------------------------------------------------------
!  InitObservData()
!
! esta rotina inicializa os dados observados para fins de 
! ajustes e calibração de parâmetros
!
! -------------------------------------------------------------
subroutine InitObsrvData (QOBS,RnflxOBS,ShflxObs,LhflxOBS,GhflxObs)

   implicit none
   include 'globcst.inc'   

   real      ::   QOBS(nsteps)
   real      ::   RnflxObs(nsteps)
   real      ::   ShflxObs(nsteps)
   real      ::   LhflxObs(nsteps)
   real      ::   GhflxObs(nsteps)

   integer   ::   i,j
   integer   ::   anoL,mesL,diamL,diajL,horaL,mintL
   real      ::   buffer1(16)                       ! Arquivo de Rn possui 16 colunas
   real      ::   buffer2(15)                       ! Arquivo de fluxos possui 15 colunas

!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
   namelist /ephstc_obsfl/  QOBSfile,RnOBSfile,SFOBSfile
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   write(6,'(1x,a)') ' Inicializando dados observados'
   REWIND(UNIT=5)
   READ (5, ephstc_obsfl)

   open(19,file=QOBSfile,form='formatted',status='unknown')
   open(20,file=RnOBSfile,form='formatted',status='unknown')
   open(21,file=SFOBSfile,form='formatted',status='unknown')
  
   do i = 1, nsteps
      read(19,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2,f16.5)') anoL, &
            mesL,diamL,diajL,horaL,mintL,QOBS(i)
!      write(6,'(i4,a,f12.6)') i,'QOBS=',QOBS(i)

      read(20,*) (buffer1(j), j=1,16)
      RnflxObs(i) = buffer1(10)

      read(21,*) (buffer2(j), j=1,15)
      ShflxObs(i) = buffer2(11)
      LhflxObs(i) = buffer2(10)
      GhflxObs(i) = buffer2(12)	    
	    
   enddo
   
   close (19)
   close (20)
   close (21)
  
end subroutine InitObsrvData



!subroutine InitObsrvData (QOBS)

!   implicit none
!   include 'globcst.inc'   

!   real             ::   QOBS(nsteps)
!   character(128)   ::   QOBSfile

!   real             ::   GhflxOBS(nsteps)
!   real             ::   ShflxOBS(nsteps)
!   real             ::   LhflxOBS(nsteps)
!   character(128)   ::   FlxOBSfile

!   integer   ::   i,j
!   integer   ::   anoL,mesL,diamL,diajL,horaL,mintL

!
!#######################################################################
!
!  Namelists:
!
!#######################################################################
!
!   namelist /ephst_/  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Inicio do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!   write(6,'(1x,a)') ' Inicializando dados observados'
!   REWIND(UNIT=5)
!   READ (5, ephst_)

!   QOBSfile='/home/afonso/workphd/dados/sgp97/vazao/sgp97_lw-vza550.dat'
!   FlxOBSfile='/home/afonso/workphd/dados/sgp97/fluxos/SfcFlux/sgp97_SfcFlux-LWF.dat'

!   open(19,file=QOBSfile,form='formatted',status='unknown')
  
!   do i = 1, nsteps
!      read(19,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2,f16.5)') anoL, &
!            mesL,diamL,diajL,horaL,mintL,QOBS(i)
!      write(6,'(i4,a,f12.6)') i,'QOBS=',QOBS(i)
!   enddo
   
!   close (19)
  
!end subroutine InitObsrvData

!end module init
