!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                     MODULE SLOPE                     ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!module slope

!interfaces
! --------------------------------------------------------------------
   
!conteúdos
! --------------------------------------------------------------------
!contains

! --------------------------------------------------------------------
!  HeadSlopeAspect
!  This computes slope and aspect using the water table elevation. 
! --------------------------------------------------------------------
subroutine HeadSlopeAspect(nx,ny,maskmap,topomap,wtlevelmap,slopemap, &
                           aspectmap,fgradmap,dirmap,totaldirmap)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'
   integer   ::   nx,ny
   integer   ::   maskmap(nx,ny)
   real(4)   ::   topomap(nx,ny)
   real(4)   ::   wtlevelmap(nx,ny)
   real(4)   ::   slopemap(nx,ny)
   real(4)   ::   aspectmap(nx,ny)
   real(4)   ::   fgradmap(nx,ny)
   integer   ::   dirmap(nx,ny,NDIRS)
   integer   ::   totaldirmap(nx,ny)
   
   integer   ::   x, xn
   integer   ::   y, yn
   integer   ::   n
   real      ::   neighbor_elev(NDIRS)
   
   integer   ::   XIndex, YIndex, valid_cell
   
! Início do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if (x==sx .and. y==sy) then
         write(6,'(a)') '   '
         write(6,'(a)') ' Entrando na rotina HeadSlopeAspect!'
   endif
   
   do y = 1, ny
   do x = 1, nx 
      if (maskmap(x,y) .eq. 1) then
         do n = 1, NDIRS
            xn = x + XIndex(n)
            yn = y + YIndex(n)
            if (valid_cell(nx,ny,xn,yn) .eq. 1) then
               if (maskmap(xn,yn) .eq. 1) then
                  neighbor_elev(n) = wtlevelmap(xn,yn)   
               else
                  neighbor_elev(n) = 0  ! OUTSIDEBASIN
               endif
            else
               neighbor_elev(n) = 0  ! OUTSIDEBASIN
            endif
         enddo
         
         call slope_aspect(wtlevelmap(x,y),neighbor_elev,slopemap(x,y),&
                           aspectmap(x,y),x,y)
         
         call flow_fractions(slopemap(x,y),aspectmap(x,y),            &
                             fgradmap(x,y),wtlevelmap(x,y),           &
			     neighbor_elev,dirmap(x,y,:),             &
			     totaldirmap(x,y),x,y)
      endif
   enddo
   enddo
   
   return

end subroutine HeadSlopeAspect


! -------------------------------------------------------------
!  ElevationSlopeAspect
!  This computes slope and aspect using the topgraph. 
! -------------------------------------------------------------
subroutine ElevationSlopeAspect(nx,ny,maskmap,topomap,wdepthmap,      &
                                slopemap,aspectmap,fgradmap,dirmap,   &
                                totaldirmap)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'
   integer   ::   nx,ny
   integer   ::   maskmap(nx,ny)
   real(4)   ::   topomap(nx,ny)
   real(4)   ::   wdepthmap(nx,ny)
   real(4)   ::   slopemap(nx,ny)
   real(4)   ::   aspectmap(nx,ny)
   real(4)   ::   fgradmap(nx,ny)
   integer   ::   dirmap(nx,ny,NDIRS)
   integer   ::   totaldirmap(nx,ny) 

   integer   ::   x, xn
   integer   ::   y, yn
   integer   ::   n, j, i
   real      ::   neighbor_elev(NDIRS)
   
   integer   ::   XIndex, YIndex, valid_cell
   real(4)   ::   FlowGrad
   
   integer   ::   Dir(nx,ny,NDIRS)
   integer   ::   TotalDir(nx,ny)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Início do código executável
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if (x==sx .and. y==sy) then
      write(6,'(a)') ' Entrando na rotina ElevationSlopeAspect!'
   endif

   ! fill neighbor array
   do y = 1, ny
   do x = 1, nx
     
      if (maskmap(x,y) .eq. 1)  then
         do n = 1, NDIRS
            xn = x + XIndex(n)
            yn = y + YIndex(n)
!	    write(102,*) x,y,xn,yn,n,XIndex(n),YIndex(n),valid_cell(nx,ny,xn,yn)

            if (valid_cell(nx,ny,xn,yn) .eq. 1) then
               if (maskmap(xn,yn) .eq. 1) then
                  neighbor_elev(n) = topomap(xn,yn)
!                   write(6,'(a,1x,4i4,2f8.2)') ' 1',xn,yn,n,maskmap(xn,yn),topomap(xn,yn),neighbor_elev(n)
               else
                  neighbor_elev(n) = 0   ! OUTSIDEBASIN
!                   write(6,'(a,1x,4i4,2f8.2)') ' 2',xn,yn,n,maskmap(xn,yn),topomap(xn,yn),neighbor_elev(n)
               endif
            else
               neighbor_elev(n) = 0   ! OUTSIDDOMAIN
            endif
         enddo

         call slope_aspect(topomap(x,y),neighbor_elev,slopemap(x,y),  &
                           aspectmap(x,y),x,y)

         call flow_fractions(slopemap(x,y),aspectmap(x,y),            &
                             fgradmap(x,y),topomap(x,y),neighbor_elev,&
                             dirmap(x,y,:),totaldirmap(x,y),x,y)
      endif
   enddo
   enddo

   return

end subroutine ElevationSlopeAspect


! --------------------------------------------------------------------
!  slope_aspect
!  Calculation of slope and aspect given elevations of cell and 
!  neighbors
! --------------------------------------------------------------------
subroutine slope_aspect(celev, nelev, slope, aspect,i,j)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'

   real      ::   celev
   real      ::   nelev(NDIRS)
   real      ::   slope
   real      ::   aspect

   integer   ::   n,i,j
   real      ::   dzdx, dzdy

!   write(6,'(a)') ' Entrando na rotina slope_aspect...!'
!   write(6,'(3f8.2)') dx,dy,celev
!   do n=1, NDIRS
!      write(6,'(1x,i4,f8.2)') n,nelev(n)
!   enddo


   select_algoritm: select case (NDIRS)

   case (8)
   ! for eight neighbors, this is exactly the 
   ! same algorithm that Arc/Info uses

      do n=1, NDIRS 
!         write(6,'(1x,i4,f8.2)') n,nelev(n)
         if (nelev(n) .eq. 0) then ! OUTSIDEBASIN
            nelev(n) = celev
         endif
      enddo

!      write(6,'(7f8.2,i4)') dzdx,nelev(8),nelev(7),nelev(6),nelev(2),   &
!                              nelev(3),nelev(4),DX
!      write(6,'(7f8.2,i4)') dzdx,nelev(8),nelev(1),nelev(2),nelev(6),   &
!                              nelev(5),nelev(4),DY


      dzdx = ((nelev(8) + 2 * nelev(7) + nelev(6)) -                  &
        (nelev(2) + 2 * nelev(3) + nelev(4))) / (8 * dx)
      dzdy = ((nelev(8) + 2 * nelev(1) + nelev(2)) -                  &
        (nelev(6) + 2 * nelev(5) + nelev(4))) / (8 * dy)


   case (4)

      if (nelev(2) .eq. 0 .and. nelev(4) .eq. 0) then
         dzdx = 0.0
      else if (nelev(2) .eq. 0) then
         dzdx = (nelev(4) - celev) / dx
      else if (nelev(4) .eq. 0) then
         dzdx = (celev - nelev(2)) / dx
      else 
         dzdx = (nelev(4) - nelev(2)) / (2 * dx)
      endif

      if (nelev(1) .eq. 0 .and. nelev(3) .eq. 0) then
         dzdy = 0.0
      else if (nelev(3) .eq. 0) then
         dzdy = (celev - nelev(1)) / dy
      else if (nelev(1) .eq. 0) then
         dzdy = (nelev(3) - celev) / dy
      else
         dzdy = (nelev(3) - nelev(1)) / (2 * dy)
      endif
    
!  default:
!     assert(0);          /* nothing else works */

  end select select_algoritm

!   if (i == sx .and. j == sy) then
!      write(6,'(a,2f18.12)') 'dzdx=', dzdx, dzdy
!      write(6,'(a,5f16.8)'), 'nelev=',nelev(1),nelev(2),nelev(3),nelev(4),celev
!   endif
    
   slope = sqrt(dzdx * dzdx + dzdy * dzdy)
   if ((dzdx .eq. 0.0) .and. (dzdy .eq. 0.0)) then
      aspect = 0.0
   else
      aspect = atan2(dzdx,dzdy)
   endif
   
!   if (i==sx .and. j==sy) then
!     write(6,'(2i4,a,2f16.8)'), i,j,'slope=',slope,aspect
!   endif   
!   write(6,'(a)') ' Saindo de slope_aspect'

return

end subroutine slope_aspect


! --------------------------------------------------------------------
!  FlowGrad()
! --------------------------------------------------------------------
real(4) function FlowGrad(slope,aspect,nelev,i,j) result(FGrad)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'

   real      ::   slope
   real      ::   aspect
   real      ::   tanB
   real      ::   nelev(NDIRS)

   real      ::   cosine
   real      ::   sine
   real      ::   total_width 
   real      ::   effective_width
   integer   ::   n,i,j

! --------------------------------------------------------------------
!   write(6,'(a)') ' Entrando na função FlowGrad...!'
! --------------------------------------------------------------------

!  inicializa constantes
! --------------------------------------------------------------------
   cosine = cos(aspect)
   sine   = sin(aspect) 

   select_algoritm: select case (NDIRS)
   
      case (4)
      
      if ((cosine .gt. 0 .and. nelev(1) .eq. 0) .or.                  &
          (cosine .lt. 0 .and. nelev(3) .eq. 0)) then
   
         cosine = -cosine
      endif
    
      if ((sine .gt. 0 .and. nelev(2) .eq. 0) .or.                    &
          (sine .lt. 0 .and. nelev(4) .eq. 0)) then

          sine = -sine
      endif
!      if (i==sx .and. j==sy) then
!         write(6,'(a)') '  '
!         write(6,'(a,6f16.8)'), 'fgrad1=',cosine,sine,nelev(1),nelev(2),nelev(3),nelev(4)
!      endif   

! compute flow widths
! --------------------------------------------------------------------
      total_width = abs(cosine) * dx + abs(sine) * dy
      tanB = sine
      FGrad = slope * total_width
      
!      if (i==sx .and. j==sy) then
!         write(6,'(a,4f16.8)') 'fgrad2=',total_width,slope,aspect,FGrad
!      endif
      
!      case (8)
!      assert(0) /*can't handle this yet*/
      
   end select select_algoritm
   
   
end function FlowGrad


! -------------------------------------------------------------
!  flow_fractions
!  Computes subsurface flow fractions given the slope and aspect 
! -------------------------------------------------------------
subroutine flow_fractions(slope,aspect,fgrad,celev,nelev,dir,total_dir,i,j)

   implicit none
   include 'grid.inc'
   include 'globcst.inc'

   real      ::   slope
   real      ::   aspect
   real      ::   fgrad
   real      ::   celev
   real      ::   nelev(NDIRS)
   integer   ::   dir(NDIRS)

   integer   ::   total_dir ! foi suprimido, da passagem de parâmetros
                                    ! pois não achei utilidade 
                                    ! em outras rotinas, por enquanto...
   real      ::   cosine
   real      ::   sine
   real      ::   total_width 
   real      ::   effective_width
   integer   ::   n,i,j
   real(4)   ::   FlowGrad


!   if (i==sx .and. j==sy) then
!      write(6,'(a)') ' Entrando na rotina flow_fractions...!'
!      write(6,'(4f8.2)') DX,DY,slope,aspect
!   endif
   

   !inicializa constantes
   cosine = cos(aspect)
   sine   = sin(aspect)
!   write(6,'(4f12.5)') slope,aspect,cosine,sine
   
   select_algoritm: select case (NDIRS)

      case (4)

      ! fudge any cells which flow outside the basin by just pointing 
      ! the aspect in the opposite direction
!      if (i == sx .and. j == sy) then 
!         write(6,'(a)')  '   '
!         write(6,'(f12.5)') celev
!         write(6,'(3f12.5)') cosine, nelev(1), nelev(3)
!         write(6,'(3f12.5)') sine, nelev(2), nelev(4)
!      endif	 
	 
      if ((cosine .gt. 0 .and. nelev(1) .eq. 0) .or.                  &
          (cosine .lt. 0 .and. nelev(3) .eq. 0)) then
      
         cosine = -cosine
      endif
    
      if ((sine .gt. 0 .and. nelev(2) .eq. 0) .or.                    &
          (sine .lt. 0 .and. nelev(4) .eq. 0)) then

         sine = -sine
      endif
!      write(6,'(2f12.5)') cosine, sine

      ! compute flow widths

      total_width = abs(cosine) * dx + abs(sine) * dy
 
!      if (i==sx .and. j==sy) then
!         write(6,'(a,3f16.8)') 'total_w=', total_width,cosine, sine
!      endif
 
      if (slope == 0.0 .and. aspect == 0.0) then
         fgrad = 0.0
      else
         fgrad = FlowGrad(slope,aspect,nelev,i,j)
      endif
      total_dir = 0

      do n=1, NDIRS
         select_width: select case (n)
            case (1)
               if (cosine .gt. 0 .and. celev > nelev(1)) then
                  effective_width = cosine * dx
               else
                  effective_width = 0.0
               endif
         
            case (3)
               if (cosine .lt. 0 .and. celev > nelev(3)) then
                  effective_width = -cosine * dx
               else
                  effective_width = 0.0
               endif

            case (2)
               if (sine .gt. 0 .and. celev > nelev(2)) then
                  effective_width = sine * dx 
               else
                  effective_width = 0.0
               endif

            case (4)
               if (sine .lt. 0 .and. celev > nelev(4)) then 
                  effective_width = -sine * dx
               else if (sine .lt. 0 .and. nelev(2) .eq. 0) then
                  effective_width = 0.0
               else
                  effective_width = 0.0		  
               endif

            !default:
            !   assert(0);      /* How can this happen? */
         end select select_width


         dir(n) = ((effective_width / total_width) * 255.0 + 0.5)
         total_dir = total_dir + dir(n)


!         if (i==sx .and. j==sy) then
!            write(6,'(a)') '  '
!	    write(6,'(2i4,4f10.0)') i,j,nelev(1),nelev(2),nelev(3),nelev(4)
!            write(6,'(1x,a,5f16.8,3i6)') 'ewidth=',effective_width, total_width,&
!                                          cosine,sine,fgrad,n,dir(n),total_dir
!         endif
      enddo

  
!      case (8)
!          assert(0)          /* can't handle this */
!
!      default:
!          assert(0)          /* other cases don't work either */

   end select select_algoritm
!   write(6,'(a)') ' Saindo de flow_fractions!'

   return

end subroutine flow_fractions



integer function XIndex (Kdir) result (xneighbor)

   implicit none
   integer               ::   Kdir
   integer,allocatable   ::   Pxneighbor(:)

   include 'grid.inc'

   if(.not. allocated(Pxneighbor))   allocate(Pxneighbor(NDIRS))
   
   if (NDIRS .eq. 4) then
      Pxneighbor = (/ 0, 1, 0, -1 /)
   elseif (NDIRS .eq. 8) then
      Pxneighbor = (/ -1, 0, 1, 1, 1, 0, -1, -1 /)
   endif

   xneighbor = Pxneighbor(Kdir)
   
   deallocate(Pxneighbor)

end function XIndex


integer function YIndex (Kdir) result (yneighbor)

   implicit none
   integer                ::   Kdir
   integer, allocatable   ::   Pyneighbor(:)

   include 'grid.inc'

   if(.not. allocated(Pyneighbor))   allocate(Pyneighbor(NDIRS))
   
   if (NDIRS .eq. 4) then
      Pyneighbor = (/ -1, 0, 1, 0 /)
   elseif (NDIRS .eq. 8) then
      Pyneighbor = (/ 1, 1, 1, 0, -1, -1, -1, 0 /)
   endif

   yneighbor = Pyneighbor(Kdir)
   
   deallocate(Pyneighbor)

end function YIndex


! -------------------------------------------------------------
!  valid_cell
!  Checks to see if grid indices, x and y, are within the grid 
!  defined by the specified Map
! -------------------------------------------------------------
integer function valid_cell(nx,ny,xn,yn)

   implicit none
   integer   ::   nx,ny
   integer   ::   xn,yn
   
   include 'grid.inc'
   
!   if (x==170 .and. y==240) then
!      write(6,'(a,2i4)') 'Entrando em valid_cell', x,y
!   endif

   if (xn .ge. 0 .and. yn .ge. 0 .and. xn .lt. NX .and. yn .lt. NY) then
      valid_cell = 1       ! INSIDEDOMAIN
   else
      valid_cell = 0       ! OUTSIDEDOMAIN
   endif

!   if (x==sx .and. y==sy) then
!      write(6,'(a,5i4)') 'valid_cell', nx,ny,xn,yn,valid_cell
!   endif
   
end function valid_cell

!end module slope
