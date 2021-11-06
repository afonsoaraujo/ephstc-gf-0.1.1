!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                      MODULE CALC                     ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module calc
   
!use slope
!
!#######################################################################
!
!  Declaração de variáveis
!
!#######################################################################
!
!   implicit none
!   include 'globcst.inc'

!
!#######################################################################
!
!  Interfaces
!
!#######################################################################
!

!
!#######################################################################
!
!  Conteúdos - Contains
!
!#######################################################################
!
!contains

subroutine CalcRootDepth(nzsoil,vegtyp,NRootlayers,RootDepth,         &
                         TotalDepth,i,j)

   implicit none
   include 'soilcst.inc'
   include 'globcst.inc'
   
   integer           ::   nzsoil
   integer           ::   vegtyp
   integer           ::   NRootlayers
   real              ::   RootDepth(nzsoil)
   real              ::   TotalDepth
   
   real              ::   RootZDepth
   integer           ::   i,j
     
     
!   if (i == sx .and. j == sy) then
!      write(6,'(a)') '   '
!      write(6,'(1x,a)') ' Invocando a função CalcRootDepth...'
!   endif

   if ( rootzone(vegtyp) > TotalDepth ) then
      RootZDepth = TotalDepth
   else 
      RootZDepth = rootzone(vegtyp)
   endif
   
   RootZDepth = RootZDepth * (TotalDepth/d3)
       
   
!   if (i == sx .and. j == sy) then
!      write(6,'(a,3f16.6,i4)') 'rootzone= ',rootzone(vegtyp),d2,d3,nzsoil
!   endif
         
   if (RootZDepth .ge. d3) then
      RootDepth(1) = d1
      RootDepth(2) = d2-d1
      RootDepth(3) = d3-d2
      NRootLayers = 3
      
!   if (i == sx .and. j == sy) then
!      write(6,'(a,3f16.8,i4)') '1-',RootDepth(1),RootDepth(2),RootDepth(3),&
!                             NRootLayers
!   endif

   else if (RootZDepth > d2 .and. rootzone(vegtyp) <=  &
               d3) then

      RootDepth(1) = d1
      RootDepth(2) = d2-d1
      RootDepth(3) = RootZDepth - d2
      NRootLayers = 3

!      if (i == sx .and. j == sy) then      
!         write(6,'(a,3f16.8,i4)') '2-',RootDepth(1),RootDepth(2),RootDepth(3),&
!                               NRootLayers
!      endif

   else if (RootZDepth > d1 .and. rootzone(vegtyp) <=  &
               d2) then
      RootDepth(1) = d1
      RootDepth(2) = RootZDepth - d1
      RootDepth(3) = 0.0
      NRootLayers = 2

!      if (i == sx .and. j == sy) then      
!         write(6,'(a,4f18.12,i4)') '3-',RootDepth(1),RootDepth(2),RootDepth(3),&
!                                   rootzone(vegtyp),NRootLayers
!      endif
      
   else 
      RootDepth(1) = RootZDepth
      RootDepth(2) = 0.0
      RootDepth(3) = 0.0
      NRootLayers = 1
      
!      if (i == sx .and. j == sy) then      
!         write(6,'(a,4f18.12,i4)') '4-',RootDepth(1),RootDepth(2),RootDepth(3),&
!                                   rootzone(vegtyp),NRootLayers
!      endif      
      
   end if
!   write(6,'(a)') ' Fim da rotina CalcRootDepth!'
   return
   
end subroutine CalcRootDepth

!---------------------------------------------------------------------
! CalcTransmissivity()
!
!  Purpose : Calculates the transmissivity through the saturated part of
!            the soil profile
!
!  Required : 
!    float SoilDepth  - Total soil depth in [m]
!    float WaterTable - Depth of the water table below the soil surface in [m]
!    float LateralKs  - Lateral hydraulic conductivity in [m/s]
!    float KsExponent - Exponent that describes exponential decay of LateralKs
!                       with depth below the soil surface [adm]
!
!  Returns  : Transmissivity in [m2/s]
!---------------------------------------------------------------------
real(4) function CalcTransmissivity (SoilDepth, WaterTable,   &
                    LateralKs, KsExponent,i,j)  result (Transmissivity)

   real(4)   ::   SoilDepth
   real(4)   ::   WaterTable
   real(4)   ::   LateralKs
   real(4)   ::   KsExponent
   
   include 'globcst.inc'
   integer   ::  i,j

!   if (i==sx .and. j==sy) then
!      write(6,'(a)') '  '
!      write(6,'(1x,a)') ' Entrando na função CalcTransmissivity...'
!     write(6,'(2f12.8)') SoilDepth, Watertable
!   endif

   if (KsExponent .eq. 0.0) then
      Transmissivity = LateralKs * (SoilDepth - WaterTable)
!      if (i==sx .and. j==sy) then
!         write(6,'(a,5f12.8)') 'T1=',KsExponent, LateralKs, SoilDepth, Watertable, Transmissivity
!      endif
   else

! Essa formulação do código do DHSVM não está legal...
!      Transmissivity = (LateralKs / KsExponent) *                     &
!      (exp(-KsExponent * WaterTable) - exp(-KsExponent * SoilDepth))
      Transmissivity = (LateralKs * SoilDepth / KsExponent) * (1 -    &
                       (Watertable/SoilDepth)**KsExponent)
!      if (i==sx .and. j==sy) then
!         write(6,'(a,5f16.8)') 'T2=',KsExponent, LateralKs, SoilDepth, Watertable, Transmissivity
!      endif
   endif
   
   Transmissivity = Transmissivity * 100  ! obtendo valores muito baixos de trasmissividade....

end function CalcTransmissivity


!---------------------------------------------------------------------
! CalcAvailableWater()
!
!  Purpose : This routine calculates the amount of soil moisture
!            available for saturated flow below the groundwater table.
!            No flow is allowed to occur if the moisture content falls
!            below the field capacity.
!
!---------------------------------------------------------------------
real(4) function CalcAvailableWater(nzsoil,soiltyp,vegtyp,TotalDepth, &
                                    TableDepth,i,j) result (AvailableWater)
   implicit none
   integer                ::   nzsoil
   integer                ::   soiltyp
   integer                ::   vegtyp
   real(4)                ::   TotalDepth
   real(4)                ::   TableDepth

   real(4),allocatable    ::   LPorosit(:)
   real(4),allocatable    ::   LFCap(:)

   integer                ::   NRootLayers
   real(4)                ::   RootDepth(nzsoil)
   real(4)                ::   DeepFCap            ! field capacity of the layer below the 
                                                   !   deepest root layer
   real(4)                ::   DeepLayerDepth      ! depth of layer below deepest root zone layer
   real(4)                ::   DeepPorosit         ! Porosit of the layer below the deepest root layer
   real(4)                ::   Depth               ! depth below the ground surface (m)
   integer                ::   i,j,k               ! counter
   integer                ::   sx,sy

   include 'soilcst.inc'

   sx = sx
   sy =  sy

!   if (i == sx .and. j == sy) then
!      write(6,'(a)')  '   '
!      write(6,'(1x,a)') ' Invocando a função CalcAvailableWater...'
!   endif

   if(.not. allocated(LPorosit))   allocate(LPorosit(nzsoil))
   if(.not. allocated(LFCap))      allocate(LFCap(nzsoil))
   
   !Atribuição de valores
   AvailableWater = 0
   Depth = 0.0

!   if (i == sx .and. j == sy) then
!      write(6,'(a,4f12.6,2i4)') 'AWater0=',AvailableWater,Depth,TotalDepth,        &
!                              TableDepth,soiltyp,vegtyp
!   endif
                        
   do k =1, nzsoil
      LPorosit(k) = Porosity(soiltyp)           ! constante ao longo de z
      LFCap(k) = wfc(soiltyp)                   ! constante ao longo de z
   enddo

   call CalcRootDepth(nzsoil,vegtyp,NRootlayers,RootDepth,TotalDepth, &
                      i,j) 
!   if (i == sx .and. j == sy) then
!      write(6,'(a,i4,f16.8)') 'NRoots=', NRootlayers, TotalDepth
!   endif

   k = 1
   do while ((k .le. NRootLayers) .and. (Depth .lt. TotalDepth))
!      if (i == sx .and. j == sy) then
!         write(6,'(a,i4,3f16.6)') 'k=',k,RootDepth(k),TotalDepth,Depth
!      endif
      
      if (RootDepth(k) .lt. (TotalDepth - Depth)) then
         Depth = Depth + RootDepth(k)
      else
         Depth = TotalDepth
      endif
      
!      if (i == sx .and. j == sy) then
!         write(6,'(a,5f12.6)') 'AWater1=',AvailableWater,Depth,TotalDepth,        &
!                                TableDepth, RootDepth(k)
!      endif      
      
      if (Depth .gt. TableDepth) then
         if ((Depth - TableDepth) .gt. RootDepth(k)) then
            AvailableWater = AvailableWater + (LPorosit(k) -         &
                             LFCap(k)) * RootDepth(k)
         else
            AvailableWater = AvailableWater + (LPorosit(k) -         &
                             LFCap(k)) * (Depth-TableDepth)
         endif

!         if (i == sx .and. j == sy) then
!            write(6,'(a,i4,5f16.6)') 'k=',k,RootDepth(k),TotalDepth,&
!                                      Depth,TableDepth,AvailableWater
!         endif            

      endif
      k = k + 1
   enddo
   
   if (Depth .lt. TotalDepth) then
      DeepPorosit = LPorosit(NRootLayers)
      DeepFCap = LFCap(NRootLayers)
      DeepLayerDepth = TotalDepth - Depth
      Depth = TotalDepth

!      if (i == sx .and. j == sy) then
!         write(6,'(a,4f18.12)') 'DP=',DeepPorosit,DeepFCap,DeepLayerDepth,&
!                                   Depth
!      endif                  
      
      if ((Depth - TableDepth) .gt. DeepLayerDepth) then
         AvailableWater = AvailableWater + (DeepPorosit - DeepFCap)  &
            * DeepLayerDepth 
      else
         AvailableWater = AvailableWater + (DeepPorosit - DeepFCap)  &
            * (Depth - TableDepth)          
      endif

!      if (i == sx .and. j == sy) then
!         write(6,'(a,f18.12)') 'AWater3=',AvailableWater
!         write(6,'(a)') '  '
!      endif                  
      
   endif
     
!  ??? assert(AvailableWater >= 0.0);

   if(allocated(LPorosit))   deallocate(LPorosit)
   if(allocated(LFCap))      deallocate(LFCap)

end function CalcAvailableWater 

!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                    MODULE WTDEPTH                    ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module WTDEPTH

!#include <assert.h>
!#include <stdio.h>
!#include <stdlib.h>
!#include "settings.h"
!#include "soilmoisture.h"

!-----------------------------------------------------------------------
real(4)  function CalcWTableDepth(NRootLayers,TotalDepth,RootDepth,   &
                                  LFCap,LPorosity,Moisture,i,j)          &
                                  result (WTableDepth)

   implicit none
   include 'soilcst.inc'
   include 'globcst.inc'

   integer   ::   NRootLayers
   real      ::   TotalDepth
   real      ::   LFCap(NRootLayers)
   real      ::   RootDepth(NRootLayers)
   real      ::   LPorosity(NRootLayers)
   real      ::   Moisture(NRootLayers)
   real(4)   ::   Deepwfc		      ! field capacity of the layer below the
	                                 ! deepest rootlayer
   real(4)   ::   DeepLayerDepth		! depth of layer below deepest root zone
				                        !    layer
   real(4)   ::   DeepPorosit		! Porosit of the layer below the deepest
				                        !   root layer
   real(4)   ::   MoistureTransfer	! amount of soil moisture transferred from
                                    !    the current layer to the layer above (m)
   real(4)   ::   TotalStorage
   real(4)   ::   Excesswfc
   real(4)   ::   TotalExcesswfc
   
   real(4)   ::   CalcMoistureTransfer

!-----------------------------------------------------------------------
   integer   ::   i,j,k
!-----------------------------------------------------------------------	
!  Includes
!-----------------------------------------------------------------------	
!   include 'globcst.inc'
!   include 'phycst.inc'
!   include 'soilcst.inc'
!-----------------------------------------------------------------------
!  Purpose    :   This function calculates the depth of the water table
!                 below the ground surface, based on the amount of soil
!                 moisture in the root zone layers. 
!
!  Required     :
!    int NRootLayers  - Number of soil lyers
!    float TotalDepth - Total depth of the soil profile (m)
!    float *RootDepth - Depth of each of the soil layers (m)
!    float *Porosit  - Porosit of each soil layer
!    float *wfc      - Field capacity of each soil layer
!    float *Adjust    - Correction for each layer for loss of soil storage
!                       due to channel/road-cut.  Multiplied with RootDepth
!                       to give the layer thickness for use in calculating
!                       soil moisture 
!
!  Returns      :
!    float WTableDepth - Depth of the water table below the soil surface
!
!  Modifies     :
!    float *Moisture - Amount of soil moisture in each layer 
!
!  Comments     :   
!    Since no unsaturated flow is allowed to occur when the moisture content
!    is below the field capacity, and because lateral saturated flow is the
!    only mechanism by which water can disappear from the soil below the
!    deepest root layer, the soil moisture content in the soil below the
!    deepest root layer can never fall below field capacity.  The water
!    immediately above the water table is assumed to be at field capacity.
!
!    Changes have been made to account for the potential loss of soil storage
!    in a grid cell due to a road-cut or channel.
!-----------------------------------------------------------------------
!   if (i == sx .and. j == sy) then
!      write(6,'(a)') '   '
!      write(6,'(1x,a)') ' Invocando a função CalcWtDepth...'
!   endif
   
   TotalExcesswfc = 0.0
   TotalStorage   = 0.0
   
   MoistureTransfer = 0.0
   DeepLayerDepth = TotalDepth
   DeepPorosit = LPorosity(NRootLayers)
   Deepwfc = LFCap(NRootLayers)  
   
!   if (i == sx .and. j == sy) then
!       write(6,'(a,f18.12,i4)') 'DLayerDepth= ', DeepLayerDepth,NRootLayers
!   endif

   do k = 1, NRootLayers
!      if (i == sx .and. j == sy) then
!         write(6,'(a,2f18.12,i4)') 'DLayerDepth= ', DeepLayerDepth,RootDepth(k),k
!      endif   
      DeepLayerDepth = DeepLayerDepth - RootDepth(k)
!      if (i == sx .and. j == sy) then
!         write(6,'(a,2f18.12,i4)') 'DLayerDepth= ', DeepLayerDepth,RootDepth(k),k
!      endif         
	end do

   MoistureTransfer=  CalcMoistureTransfer(NRootLayers,TotalDepth,    &
                                           RootDepth,LFCap,LPorosity, &
                                           Moisture(1),Moisture(2),   &
                                           Moisture(3),i,j)
    
!   if (i == sx .and. j == sy) then
!       write(6,'(a,f18.12)') 'Mtransfer= ', MoistureTransfer
!       write(6,'(a)') '  '
!   endif   

   if (MoistureTransfer > 0.0) then
!      /* Surface ponding occurs */
!-----------------------------------------------------------------------
      WTableDepth = - MoistureTransfer
   else
!-----------------------------------------------------------------------
!    /* Warning added by Pascal Storck, 08/15/2000 */
!    /* Based on a single bad parameter in a DHSVM input file (a third layer
!       vertical hydraulic conductivity that was 10 times smaller than the layer
!       above it), it was noted that DHSVM can develop what are basically perched
!       water tables.  These perched water tables greatly complicate the calculation
!       of the pixel water table depth because the soil below the perched table
!       is not completely saturated above field capacity.
!       For example, if we have three soil layers and a deep layer, all 1 meter thick,
!       and we saturate the second layer from the surface, what is, or should be, the water
!       table depth. Should we allow subsurface flow to occur, should we include the saturation
!       of disconnected overlying layers in the calculation of the hydraulic gradient.
!       At this point, just be cautious.  Using any combination of soil parameters or 
!       intial water states which can cause the lower layers of the soil profile
!       to drain more quickly than water can flow down through the matrix will
!       result in mass balance problems.  I.e. water will be forced out of the cell
!       to the downslope, this water will be taken from the deepest soil layer, which
!       can cause the deep layer soil moisture to go negative. */
!-----------------------------------------------------------------------
      
      TotalStorage = TotalStorage + DeepLayerDepth * (DeepPorosit     &
                      - Deepwfc)

      Excesswfc = DeepLayerDepth * (Moisture(3) - Deepwfc)

!      if (i == sx .and. j == sy) then
!          write(6,'(a,6f18.12)') 'TStorage0=',TotalStorage,DeepLayerDepth,&
!		                           DeepPorosit,Deepwfc,Moisture(3),&
!                                 Excesswfc
!      endif   


      if (Excesswfc < 0.0) then
         Excesswfc = 0.0
		end if
    
	   TotalExcesswfc = TotalExcesswfc + Excesswfc

      do k = 1, NRootLayers
         TotalStorage = TotalStorage + RootDepth(k) * (LPorosity(k)   &
                        - LFCap(k))
				
         Excesswfc = RootDepth(k) * (Moisture(k) - LFCap(k))
         
			if (Excesswfc < 0.0) then
	         Excesswfc = 0.0
			end if
         TotalExcesswfc = TotalExcesswfc + Excesswfc
!         if (i == sx .and. j == sy) then
!            write(6,'(a,6f18.12)') 'TStorage1=',TotalStorage,RootDepth(k),&
!		                             LPorosity(k),LFCap(k),Moisture(k),Excesswfc
!         endif            
			
      end do
      
		WTableDepth = TotalDepth * (1 - TotalExcesswfc / TotalStorage)
!      if (i==sx .and. j==sy) then         
!         write(6,'(a,4f20.12)') 'WTableDepth1',WTableDepth,TotalDepth,TotalExcesswfc,TotalStorage
!      endif      
      
		if (WTableDepth < 0) then
         WTableDepth = -(TotalExcesswfc - TotalStorage)
!         if (i==sx .and. j==sy) then         
!            write(6,'(a,5f16.6)') 'WTableDepth11',WTableDepth,TotalDepth,TotalExcesswfc,TotalStorage
!         endif
		end if
   end if
   
!   assert(WTableDepth <= TotalDepth)
  
end function CalcWTableDepth


!-----------------------------------------------------------------------
real(4)  function CalcMoistureTransfer (NRootLayers,TotalDepth,       &
                                        RootDepth,LFCap,LPorosity,    &
                                        wetsfc,wetdp2,wetdp3,i,j)     &
                                        result (MoistureTransfer)
   implicit none
   include 'soilcst.inc'

   integer   ::   NRootLayers
   real      ::   TotalDepth
   real      ::   LFCap(NRootLayers)
   real      ::   RootDepth(NRootLayers)
   real      ::   LPorosity(NRootLayers)
   real      ::   wetsfc
   real      ::   wetdp2
   real      ::   wetdp3
   real(4)   ::   Deepwfc		      ! field capacity of the layer below the
	                                 ! deepest rootlayer
   real(4)   ::   DeepLayerDepth		! depth of layer below deepest root zone
				                        !    layer
   real(4)   ::   DeepPorosit		   ! Porosit of the layer below the deepest
				                        !   root layer

!-----------------------------------------------------------------------
   integer   ::   i,j,k
   integer   ::   sx,sy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   sx=sx
   sy=sy
!-----------------------------------------------------------------------
!  /* Redistribute soil moisture.  I.e. water from supersaturated layers 
!     is transferred to the layer immediately above */
!-----------------------------------------------------------------------
   
   if ( wetdp3 >= LPorosity(3) ) then
		   
      MoistureTransfer = (wetdp3 - LPorosity(3)) * d3 
			  
      wetdp3 = LPorosity(3)

      wetdp2 = wetdp2 + MoistureTransfer / d2 
      if ( wetdp2 >= LPorosity(2) ) then
         MoistureTransfer = (wetdp2 - LPorosity(2)) * d2
         wetdp2 = LPorosity(2)
      else
         MoistureTransfer = 0.0
      endif
      
      wetsfc = wetsfc + MoistureTransfer / d1
      if ( wetsfc >= LPorosity(1) ) then
         MoistureTransfer = (wetsfc - LPorosity(1)) * d1
         wetsfc = LPorosity(1)
      else
         MoistureTransfer = 0.0
      endif

   else if ( wetdp2 >= LPorosity(2) ) then
		   
      MoistureTransfer = (wetdp2 - LPorosity(2)) * d2 
			  
      wetdp2 = LPorosity(2)
      
      wetsfc = wetsfc + MoistureTransfer / d1
      if ( wetsfc >= LPorosity(1) ) then
         MoistureTransfer = (wetsfc - LPorosity(1)) * d1
         wetsfc = LPorosity(1)
      else
         MoistureTransfer = 0.0
      endif
   else
      MoistureTransfer = 0.0
   end if

end function CalcMoistureTransfer


!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION GETLAI                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

real(4) function CalcLaiNDVI(VegTyp,NDVI) result(LAI)

!SUBROUTINE CalcLaiNDVI(VegMap,VType)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the Leaf Area Index (lai) from NDVI data and interpolate
!  them into the model domain.
!
!  The calculation of LAI depends on vegetation type: grass or tree.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!
!  3/22/94
!
!  ADPTAÇÃO: Afonso Auguato M. de Araujo
!  5/11/2004
!  Esta rotina foi adaptada a partir do código do modelo ARPS.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!
!  vtyp     Vegetation type;
!  ndvi     NDVI data array
!
!  OUTPUT:
!
!  lai      LAI array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
   implicit none
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
   integer                        ::   VegTyp
   real(4)                        ::   NDVI
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   if ( VegTyp == 9 .OR. VegTyp == 14)  then  ! água e gelo
      LAI = 0.0
!      write(6, '(a,i4,2f10.2)') 'VType=',VegTyp, LAI, NDVI
   else if ( VegTyp == 1 .OR. VegTyp == 2 .OR. VegTyp == 3 .OR.       &
      VegTyp == 4 .OR. VegTyp == 5 .OR. VegTyp == 9 .OR.              &
      VegTyp == 10 .OR. VegTyp == 11 .OR. VegTyp == 12 .OR.           &
      VegTyp == 13 ) then

      if (NDVI .le. 0.0) then
         NDVI = 0.0
      endif
      
      LAI = - LOG( ( 1. - NDVI/.915 ) / .83 ) / 0.96
!      write(6, '(a,i4,2f10.2)') 'VType1=',VegTyp, LAI, NDVI            

   else
      
      if (NDVI .eq. 0.0) then
         NDVI = 0.00000001
      end if

      LAI = 1.625 * EXP( NDVI / 0.34 )
!      write(6, '(a,i4,2f10.2)') 'VType1=',VegTyp, LAI, NDVI            

   end if

   LAI = MAX( LAI, 0.00000001 )

   end function CalcLaiNDVI

!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                  SUBROUTINE CALCRNFLX                ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine CalcRnflx (nx,ny,vegmap,tair,tsfc,radsw,radlw,rnflx,qvair) 

!-----------------------------------------------------------------------
!
!  OBJETIVO - PURPOSE:
!
!  Calculo do fluxo de radiação líquida
!
!-----------------------------------------------------------------------
   implicit none
   
   integer   ::   nx,ny
   integer   ::   vegmap (nx,ny)
   real      ::   tair   (nx,ny)
   real      ::   tsfc   (nx,ny)
   real      ::   radsw  (nx,ny)
   real      ::   radlw  (nx,ny)   
   real      ::   rnflx  (nx,ny)
   real      ::   qvair  (nx,ny)
   
   real      ::   emissac        ! emissividade calculada segundo Brutsaert
   
   include 'globcst.inc'
   include 'phycst.inc'
   include 'soilcst.inc'
   
   integer   ::   i,j,k

   ! Afonso Araujo 11/01/2005 - runopt = 1 => Rodada desacoplada
   do j=1, ny
   do i=1, nx
   
      emissac = 0.643 * (qvair(i,j)/tair(i,j))**(1/7)		       ! Brutsaert

!      rnflx(i,j) = radsw(i,j)*(1-albdveg(vegmap(i,j)))+emissa*sbcst*(tair(i,j)**4)  &  ! ARPS
!                    - emissg*sbcst*(tsfc(i,j)**4)
		    
      rnflx(i,j) = radsw(i,j)*(1-albdveg(vegmap(i,j))) + emissg*1.18*radlw(i,j)  & 
                    - emissg*sbcst*(tsfc(i,j)**4)
		    
!      rnflx(i,j) = radsw(i,j)*(1-albdveg(vegmap(i,j))) + emissac*sbcst*(tair(i,j)**4)  & 
!                    - emissg*sbcst*(tsfc(i,j)**4)		    
                    
                   
!      if (i == sx .and. j == sy) then
!         write(6,'(a)') '   '
!	 write(6,'(2i4)') i,j
!         write(6,'(a,7f18.8,i4)') 'RnFlx=',rnflx(i,j),radsw(i,j),albdveg(vegmap(i,j)),&
!                             qvair(i,j),tair(i,j),tsfc(i,j),emissac,vegmap(i,j)
!      endif
                    
   
   enddo
   enddo

   return

end subroutine CalcRnflx


!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                  SUBROUTINE CALCqvair                ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine Calcqvair (nx,ny,pres,tair,relh,qvair)

!-----------------------------------------------------------------------
!
!  OBJETIVO - PURPOSE:
!
!  Calculo do fluxo de radiação líquida
!
!-----------------------------------------------------------------------
   implicit none
   
   integer   ::   nx,ny
   real      ::   pres   (nx,ny)
   real      ::   tair   (nx,ny)
   real      ::   relh   (nx,ny)
   real      ::   qvair  (nx,ny)
   real      ::   qvsata (nx,ny)

   integer   ::   i,j
	
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Início do código executável - Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
   CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, pres,tair,qvsata)
	 
   do j = 1, ny-1
   do i = 1, nx-1
      qvair(i,j) = relh(i,j) * qvsata(i,j)
   enddo
   enddo
   
   return

end subroutine Calcqvair


!end module calc


