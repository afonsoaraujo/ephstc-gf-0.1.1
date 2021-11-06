!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                     MODULE ROUTE                     ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!module route
   

!declaração de variáveis 


!interfaces
!interface IRoute
!   module procedure CalcTransmissivity,CalcAvailableWater
!   module procedure XIndex
!end interface   


!conteúdos
!contains

!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                SUBROUTINE ROUTESUBSFC                ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
subroutine RoutSubSfc (dt,nx,ny,nzsoil,maskmap,topomap,aspectmap,     &
                       slopemap,fgradmap,dirmap,soilmap,vegmap,       &
                       sdepthmap,wtdepthmap,wtlevelmap,OutFlow,       &
                       SatFlow)

   implicit none   
   include 'grid.inc'
   include 'globcst.inc'
   
   integer            ::   dt
   integer            ::   nx,ny,nzsoil
   integer            ::   maskmap(nx,ny)
   real               ::   topomap(nx,ny)
   real               ::   aspectmap(nx,ny)
   real               ::   slopemap(nx,ny)
   real               ::   fgradmap(nx,ny)
   integer            ::   dirmap(nx,ny,NDIRS)   
   integer            ::   soilmap(nx,ny)
   real               ::   sdepthmap(nx,ny)
   integer            ::   vegmap(nx,ny)
   real               ::   wtdepthmap(nx,ny)
   real               ::   wtlevelmap(nx,ny)   
   real               ::   OutFlow(nx,ny)
!   real               ::   InFlow(nx,ny)
   real               ::   Satflow(nx,ny)

!  variáveis locais
!---------------------------------------------------------------------   
   integer            ::   i,j,k        ! counters
	integer            ::   nxx,nyy
   integer            ::   Layer
   real(4)            ::   BankHeight
   real(4)            ::   Adjust
   real(4)            ::   fract_used
   real(4)            ::   Depth
   real(4)            ::   water_out_road;
   real(4)            ::   Transmissivity
   real(4)            ::   AvailableWater
!  External functions
!---------------------------------------------------------------------
   real(4)            ::   CalcTransmissivity
   real(4)            ::   CalcAvailableWater
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
!   write(6,'(a)')    '/-----------------------------------------------/   '
!   write(6,'(1x,a)') '/---- Entrando na rotina RoutSubSfc-------------/'
!   
!#######################################################################
!
!  Aloca memória para variáveis - Allocate memory for variables
!
!#######################################################################
!-----------------------------------------------------------------------
!  /* reset the saturated subsurface flow to zero */
!-----------------------------------------------------------------------
   do j = 1, ny 
   do i = 1, nx
      if (maskmap(i,j) .eq. 1) then
         OutFlow(i,j) = 0.0
         SatFlow(i,j) = 0.0
      endif
   enddo
   enddo
!-----------------------------------------------------------------------
!  /* next sweep through all the grid cells, calculate the amount of
!     flow in each direction, and divide the flow over the surrounding
!     pixels */
!-----------------------------------------------------------------------
   do j = 1, ny
   do i = 1, nx

      if (maskmap(i,j) .eq. 1) then
		
         Adjust = 1.0                 ! Network(y,x)%Adjust
!-----------------------------------------------------------------------
! Determina as frações de escoamento em cada direção
!-----------------------------------------------------------------------
         fract_used = 0.0
         do k = 1, NDIRS
            fract_used = fract_used + dirmap(i,j,k)
!            if (i==sx .and. j==sy) then
!               write(6,'(a,3i4,f16.8,i4)') 'RSub-',i,j,k,fract_used,dirmap(i,j,k)
!            endif
         enddo
         fract_used = fract_used / 255.0 
!        write(6,'(f12.6)') fract_used
!-----------------------------------------------------------------------
!  only bother calculating subsurface flow if water table is above bedrock
!-----------------------------------------------------------------------
!         if (i==sx .and. j==sy) then
!            write(6,'(a,3f16.6)') 'R1-wtdep=', wtdepthmap(i,j),sdepthmap(i,j),wtlevelmap(i,j)
!         endif
         
         if (wtdepthmap(i,j) .lt. sdepthmap(i,j)) then
            if (wtdepthmap(i,j) .gt. bankheight) then
               Depth = wtdepthmap(i,j)
            else
               Depth = bankheight
            endif
!            write(6,'(1x, f8.4)') Depth
            
            Transmissivity= CalcTransmissivity(sdepthmap(i,j),Depth,  &
                            KsLat(soilmap(i,j)),                      &
                            KsLatExp(soilmap(i,j)),i,j)

            Outflow(i,j) = (Transmissivity*fract_used*fgradmap(i,j) * &
                            dt) / (dx * dy)

!            if (i==sx .and. j==sy) then
!               write(6,'(a,4f16.12,a,i6)') 'R2-', Transmissivity,fgradmap(i,j), &
!                     fract_used,Outflow(i,j),'  soilmap=', soilmap(i,j)
!            endif
!---------------------------------------------------------------------
!  check whether enough water is available for redistribution */
!---------------------------------------------------------------------
            AvailableWater = CalcAvailableWater(nzsoil,soilmap(i,j),  &
                                vegmap(i,j),sdepthmap(i,j),           &
                                wtdepthmap(i,j),i,j)

!            if (i==sx .and. j==sy) then
!               write(6,'(a,4f16.12)') 'R3-',OutFlow(i,j),AvailableWater,&
!                                            wtdepthmap(i,j),sdepthmap(i,j)
!            endif

            if ( OutFlow(i,j) .gt. AvailableWater) then
               OutFlow(i,j) = AvailableWater
            else
               Outflow(i,j) = Outflow(i,j)
            endif
         else
            Depth = sdepthmap(i,j)
            OutFlow(i,j) = 0.0
         endif

!---------------------------------------------------------------------
!  Subsurface Component - Decrease water change by outwater
!---------------------------------------------------------------------
!         if (i==sx .and. j==sy) then
!            write(6,'(a,3f16.8,2i4)') 'SatFlow0=',SatFlow(i,j),Outflow(i,j),AvailableWater,i,j
!         endif

         SatFlow(i,j) = SatFlow(i,j) - OutFlow(i,j)     ! aqui, SatFlow deve ser sempre negativo!
         
!         if (i==sx .and. j==sy) then
!            write(6,'(a,2f16.8,2i4)') 'SatFlow0=', SatFlow(i,j),Outflow(i,j),i,j
!         endif

!---------------------------------------------------------------------
!  Assign the water to appropriate surrounding pixels
!---------------------------------------------------------------------
         OutFlow(i,j) = OutFlow(i,j)/255.0

!         if (i==sx .and. j==sy) then
!            write(6,'(a)') '  '
!            write(6,'(2i4,a,6f18.12)') i,j,'  SatFlow=',SatFlow(i,j),Outflow(i,j),wtdepthmap(i,j),sdepthmap(i,j),wtlevelmap(i,j),AvailableWater
!         endif         

         do k = 1, NDIRS
            nxx = XIndex(k) + i
            nyy = YIndex(k) + j
            if (valid_cell(nx,ny,nxx,nyy) .eq. 1) then

               SatFlow(nxx,nyy) = SatFlow(nxx,nyy) + OutFlow(i,j) *        &
                                  dirmap(i,j,k)

            endif
!            if (i==sx .and. j==sy) then
!                write(6,'(a)') '  '	    
!               write(6,'(3i4,3f16.8)') k,valid_cell(nx,ny,nxx,nyy),dirmap(i,j,k),SatFlow(nxx,nyy),Outflow
!                write(6,'(3i4,2f16.8,2i6)') k,i,j,wtlevelmap(i,j),wtlevelmap(nxx,nyy),dirmap(i,j,k),dirmap(nxx,nyy,k)
!            endif
         enddo
!         if (i==sx .and. j==sy) then
!            write(6,'(a,4f18.12,2i4)') 'SatFlow3=', SatFlow(i,j),OutFlow(i,j),OutFlow(sx,sy-1),OutFlow(sx-1,sy),i,j
!         endif         
!         if (i==sx-1 .and. j==sy-1) then
!            write(6,'(a,2f18.12,2i4)') 'SatFlow30=', SatFlow(sx,sy),OutFlow(i,j),i,j
!         endif
!         if (i==sx .and. j==sy-1) then
!            write(6,'(a,2f18.12,3i4)') 'SatFlow31=', SatFlow(sx,sy),OutFlow(i,j),dirmap(i,j,3),i,j
!         endif
!         if (i==sx-1 .and. j==sy) then
!            write(6,'(a,2f18.12,3i4)') 'SatFlow32=', SatFlow(sx,sy),OutFlow(i,j),dirmap(i,j,2),i,j
!         endif
!         if (i==sx+1 .and. j==sy) then
!            write(6,'(a,2f18.12,3i4)') 'SatFlow33=', SatFlow(sx,sy),OutFlow(i,j),dirmap(i,j,4),i,j
!         endif                                    
!         if (i==sx .and. j==sy+1) then
!            write(6,'(a,2f18.12,3i4)') 'SatFlow34=', SatFlow(sx,sy),OutFlow(i,j),dirmap(i,j,1),i,j
!         endif
      endif

   enddo
   enddo
   
!   write(6,'(1x,a)') ' Saindo da rotina RoutSubSfc...'
!   write(6,'(a)') '   '

return

end subroutine RoutSubSfc


!end module route
