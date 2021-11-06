!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                    MODULE ROUTESFC                   ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!module routesfc
   

!declaração de variáveis 

!interfaces
!interface IRoute
!   module procedure CalcTransmissivity,CalcAvailableWater
!   module procedure XIndex
!end interface   

!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                SUBROUTINE ROUTESFCTOPM               ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine RoutSfcTOPM (nx,ny,IT,DT,maskmap,lnatbmap,RunOff,TArea,NAC,&
                        IA,AC,ST,SUMAC,TL,NCH,ACH,D,ND,NR,AR,Tlhflx, &
			TPrecip,TSatFlow,TRunOff,Q)
			

   implicit none   
   include 'grid.inc'
   include 'globcst.inc'
   include 'topmodel.inc'
      
!  variáveis principais de conecao com o EPHST
!---------------------------------------------------------------------   
   integer            ::   nx,ny
   integer            ::   IT
   integer            ::   DT
   real               ::   lnatbmap(nx,ny)
   integer            ::   maskmap(nx,ny)
   real               ::   RunOff(nx,ny)     ! Runoff distribuido   
   
   real               ::   Tlhflx            ! Evaporacao media na bacia
   real               ::   TPrecip           ! Precipitacao media na bacia
   real               ::   TSatFlow          ! Escoamento de base medio na bacia
   real               ::   TRunOff
   
   real               ::   Q(nsteps)
   real               ::   QOF
   real               ::   QOUT
   real               ::   SUMQ


!  variáveis locais
!---------------------------------------------------------------------   
   integer            ::   i,j,k        ! counters
   
!  variáveis do TOPMODEL
!---------------------------------------------------------------------      
   integer            ::   NSC,IMAP,IOUT,NAC
   integer            ::   IA,IB,IR
   character(80)      ::   SUBCAT
   real               ::   ACF,ACM,ACMAX
   real               ::   TArea,SUM,SUMAC
   integer            ::   TL,NCH
   real               ::   OF
   real               ::   EX(NCHMX)   ! recebe TRunoff agregado pelo indice topografico..
   real               ::   AC(NACMX)
   real               ::   ST(NACMX)
   real               ::   ACH(NACMX)  ! colocar no numero de classes maximo de IT
   real               ::   D(NACMX)    ! nos arquivos .inc ou input?
   integer            ::   ND,NR,ihour(NACMX)
   real               ::   AR(NCHMX)
   real               ::   CA(nsteps)
   real               ::   TRun

!  External functions
!---------------------------------------------------------------------

!  Diversos
!---------------------------------------------------------------------
   integer   ::   coef
   LOGICAL   ::   firstcall        ! First call flag of this subroutine
   SAVE firstcall
   DATA firstcall/.true./ 

!  Variáveis de depuração - Debug variables
!---------------------------------------------------------------------
   integer   ::   istatus
!---------------------------------------------------------------------
!  Includes
!---------------------------------------------------------------------
   include 'soilcst.inc'
   include 'alloc.inc'      ! Memory allocation declaration and interfaces
!
!#######################################################################
!  Listas de nomes - Namelists:
!#######################################################################
!
   namelist /ephst_topmodel/  CHVDT,RVDT,Q0,itfile
!
!#######################################################################
! 
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!
!#######################################################################
!
   write(6,'(a)')    '/-----------------------------------------------/   '
   write(6,'(1x,a)') '/---- Entrando na rotina RoutSfcTOPM------------/'

   if (firstcall) then

! Reinitialise discharge array
!---------------------------------------------------------------------
      SUM=0.
      Q=0.
      write(6,'(a,4f18.6)') '1 ',Q(1),Q0,CHVDT,TArea
      DO  I=1, ND
         Q(I) = Q(I) + Q0*TArea
  	 write(6,'(a,i8,2f18.6)') 'Q= ', I,Q(I),Q0
      enddo
      
      write(6,'(a,3i8)') 'IN-00= ',IN,ND,NR
      DO  I=1, NR
         SUM=SUM+AR(I)
         IN = ND + I 
         Q(IN)=Q(IN)+Q0*(TArea-SUM)
         write(6,'(i8,a,4f18.6,3i4)') IN,'  Q(in)= ',Q(IN),Q0,SUM,AR(I),NR,ND,I
      enddo
      
      IN = 0.
      
      firstcall = .false.
   endif  ! endif of firstcall

      
! AGREGA RUNOFF EM FUNCAO DO INDICE TOPOGRAFICO
!---------------------------------------------------------------------
!***Afonso Rio 26/01/06

      TRun=0.
      do IA=1,NAC
         EX(IA) = 0.
      enddo

      if (TRunOff > 0.0000001) then
      do j=1,ny
      do i=1,nx
         if (maskmap(i,j) == 1) then
         do IA = 1,NAC
	    if (IA > 1 .and. lnatbmap(i,j) < IA &
	               .and. lnatbmap(i,j) > (IA-1)) then
		       
               EX(IA) = EX(IA) + RunOff(i,j)*AC(IA)
!               write(6,'(2i4,a,4f18.8,3i4)') it,IA,' EX=',EX(IA),TRunOff,RunOff(i,j),lnatbmap(i,j),maskmap(i,j),i,j
	    endif

	 enddo
	 endif
      enddo
      enddo

      do IA=1,NAC
         TRun = TRun + EX(IA)
      enddo
      write(6,'(a,f18.8)') ' TRun=',TRun
      endif
      
!  Initialise contributing area counts
      do IA = 1, NAC
         ihour(ia)=0
      enddo            
      

!  START LOOP ON TIME STEPS

!      DO IT=1,NSTEP  ! Esse time step deve ser coordenado com o do EPHST
         
	 QOF=0.
         ACM=0.
!  START LOOP ON A/TANB INCREMENTS
         DO IA=1,NAC

            ACF = 0.5*(AC(IA)+AC(IA+1))

!  CALCULATION OF FLOW FROM FULLY SATURATED TArea
!  This section assumes that a/tanB values are ordered from high to low
!
            OF = 0.
!	    write(6,'(a,i8,4f18.6)') 'ACF= ',IA,ACF,AC(IA),ACM,EX(IA)	    
            IF(IA.GT.1)THEN
               IB=IA-1
	       IF(EX(IA) > 0.00000001)THEN
!  Both limits are saturated
		  OF = AC(IA)*(EX(IB)+EX(IA))/2
		  ACM=ACM+ACF
		  ihour(ib) = ihour(ib) + 1
                  write(6,'(a,6f18.6,3i8)') 'OF-1= ',OF,AC(IA),ACM,ACF,EX(IA),EX(IB),ihour(ib),IA,IB
	       ELSE
!  Check if lower limit saturated (higher a/tanB value)
	          IF(EX(IB).GT.0.00000001)THEN
                     ACF = ACF * EX(IB) / (EX(IB) - EX(IA))
		     OF=ACF * EX(IB)/2
		     ACM = ACM + ACF
		     ihour(ib) = ihour(ib) + 1
                     write(6,'(a,6f18.6,3i8)') 'OF-2=',OF,AC(IA),ACM,ACF,EX(IA),EX(IB),ihour(ib),IA,IB		     
	          ENDIF
	       ENDIF
            ENDIF
            QOF=QOF+OF   
!            write(6,'(a,2f18.6)') 'QOF= ',QOF,OF	    

!  Set contributing area plotting array
            CA(IT) = ACM
            IF(ACM.GT.ACMAX)ACMAX=ACM

         ENDDO  !  END OF A/TANB LOOP
	 

!  CALCULATE SATURATED ZONE DRAINAGE
!         QB=SZQ*EXP(-SBAR/SZM)   ! Essas variaveis do TOPMODEL sao substituidas pelas
!         SBAR=SBAR-QUZ+QB        ! calculadas pelo DHSVM
         QOUT = TSatFlow + QOF
         SUMQ = SUMQ + QOUT
         coef = 625270016/3600

         write(6,'(a)') '     '      
         write(6,'(a,6i8)') 'IN-0= ',IN,ND,IR,IT,nsteps,NR      	 
!  CHANNEL ROUTING CALCULATIONS
!  allow for time delay to catchment outlet ND as well as 
!  internal routing array

         if (IN.lt.nsteps) then
            DO IR = 1, NR
               IN=IT+ND+IR-1
               Q(IN)=Q(IN)+QOUT*AR(IR)
               write(6,'(i8,f18.8,a,5f18.8,2i8)') it,Q(it),'  Q(IN)= ',Q(IN),&
                                  QOUT,TSatFlow,QOF,AR(IR),IN,IR      
            ENDDO      
	 endif
!         write(6,'(a,5i8)') 'IN-1= ',IN,ND,IR,IT,nsteps      	 	 
!         IF(IN.GT.nsteps) WRITE(6,'(a)') '**Aviso, IN > NSTEPS **'    !GO TO 10  !???	 

          write(17,'(1x,i4,10f20.12)')it,Tlhflx,Q(it),QOUT,TSatFlow,QOF

!      ENDDO  !  END OF TIME STEP LOOP -- deve ser coordenado com o EPHST


end subroutine RoutSfcTOPM
   


!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                    SUBROUTINE DA2TLH                 ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine DA2TLH (NCH,ACH,D,TArea,ND,NR,AR)

   implicit none   
   include 'topmodel.inc' 
   
   integer   ::   NCH
   real      ::   ACH(NCHMX)   
   real      ::   D(NCHMX)
   real      ::   TArea
   real      ::   AR(NCHMX)   
   real      ::   TCH(NCHMX)
   integer   ::   NR,ND,TIME

   real      ::   A1,A2,SUMAR

!  variáveis locais
!---------------------------------------------------------------------   
   integer            ::   I,J,IR        ! counters
!
!#######################################################################
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!#######################################################################
   write(6,'(a)')    '/-----------------------------------------------/   '
   write(6,'(1x,a)') '/---- Entrando na rotina DA2TCH-------------/'
!  CONVERT DISTANCE/TArea FORM TO TIME DELAY HISTOGRAM ORDINATES
!  
      TCH(1) = D(1)/CHVDT
      write(6,'(a,3f18.6)') 'D1= ', D(1),CHVDT,TCH(1)
            
      DO J = 2,NCH
         TCH(J) = TCH(1) + (D(J) - D(1))/RVDT
	 write(6,'(a,3f18.6)') 'D= ', D(J),CHVDT,TCH(J)
      enddo
      
      NR = INT(TCH(NCH))
      IF (FLOAT(NR).LT.TCH(NCH))NR=NR+1
      ND = INT(TCH(1))
      NR = NR - ND
      write(6,'(a,2i4,f18.6)') 'NR= ', NR,ND,TCH(NCH)
      
      
      DO 20 IR=1,NR
         TIME = ND+IR
         IF(TIME.GT.TCH(NCH))THEN
	    AR(IR)=1.0
         ELSE
	 
	    DO J=2,NCH
	       IF(TIME.LE.TCH(J))THEN
	          AR(IR)=ACH(J-1)+(ACH(J)-ACH(J-1))*(TIME-TCH(J-1))/  &
                     (TCH(J)-TCH(J-1))
	          GOTO 20
	       ENDIF
            enddo
         ENDIF
   20 CONTINUE   ! para modernizar completamente essa rotina e preciso 
                 ! eliminar esse goto.....Por hora vamos ver se isso 
		 ! funciona....

      A1= AR(1)
      SUMAR=AR(1)
      AR(1)=AR(1)*TArea
      
      IF(NR.GT.1)THEN
         DO IR=2,NR
	    A2=AR(IR)
	    AR(IR)=A2-A1
	    A1=A2
	    SUMAR=SUMAR+AR(IR)
	    AR(IR)=AR(IR)*TArea
         enddo
      ENDIF

end subroutine DA2TLH


!conteúdos
!contains

!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                  SUBROUTINE ROUTESFC                 ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine RoutSfcDHSVM1 (dt,nx,ny,maskmap,dirmap,totaldirmap,RunOff)

   implicit none   
   include 'grid.inc'
   include 'globcst.inc'
   
   integer            ::   dt
   integer            ::   nx,ny
   integer            ::   maskmap(nx,ny)
   integer            ::   dirmap(nx,ny,NDIRS)
   integer            ::   totaldirmap(nx,ny)   
   real               ::   RunOff(nx,ny)

   real,allocatable   ::   surface(:,:)

!  variáveis locais
!---------------------------------------------------------------------   
   integer            ::   i,j,k        ! counters
	integer            ::   nxx,nyy

!  External functions
!---------------------------------------------------------------------
   integer            ::   XIndex,YIndex,valid_cell

!  Variáveis de depuração - Debug variables
!---------------------------------------------------------------------
   integer   ::   istatus
   integer   ::   teste

!---------------------------------------------------------------------
!  Includes
!---------------------------------------------------------------------
   include 'soilcst.inc'
   include 'alloc.inc'      ! Memory allocation declaration and interfaces
!
!#######################################################################
! 
!  Incício do código executável / BEGINING OF THE EXECUTAVEL CODE
!
!#######################################################################
!
   write(6,'(a)')    '/-----------------------------------------------/   '
   write(6,'(1x,a)') '/---- Entrando na rotina RoutSfc----------------/'
!   
!#######################################################################
!
!  Aloca memória para variáveis - Allocate memory for variables
!
!#######################################################################
   if (.not. allocated(surface)) ALLOCATE(surface  (nx,ny),STAT=istatus)
   CALL check_alloc_status(istatus, "ephst:surface")
   surface = 0

   do j = 1, ny
   do i = 1, nx
   
      if (maskmap(i,j) == 1) then
         surface(i,j) = RunOff(i,j)
         RunOff(i,j) = 0.0
      endif
      
   enddo
   enddo
   
   do j = 1, ny
   do i = 1, nx   
  
      if (maskmap(i,j) .eq. 1) then
         
         do k = 1, NDIRS
            nxx = XIndex(k) + i
            nyy = YIndex(k) + j
            if (valid_cell(nx,ny,nxx,nyy) .eq. 1) then
               if (i==sx .and. j==sy) then
                  write(6,'(a,i4,f18.12)') 'ka=', k,RunOff(nxx,nyy)
               endif               

               RunOff(nxx,nyy) = RunOff(nxx,nyy) + surface(i,j) *     &
               dirmap(i,j,k) / totaldirmap(i,j)
               
               if (i==sx .and. j==sy) then
                  write(6,'(a,3i4,2f18.12)') 'kd=', k,dirmap(i,j,k),&
                  totaldirmap(i,j),surface(i,j),RunOff(nxx,nyy)
               endif               
            endif
         enddo
         
      endif
   enddo
   enddo

   write(6,'(a)')    '/-----------------------------------------------/   '
   write(6,'(1x,a)') '/---- Saindo da rotina RoutSfc------------------/'

return

end subroutine RoutSfcDHSVM1


!end module route
