program geramaskbin

   implicit none
   integer           ::   nx,ny,i,j
   parameter (nx=572, ny=468)
   
   character(128)    ::   maskfile,maskout
   integer           ::   maskmap(nx,ny)
   
   maskfile="/home/afonsoaraujo/work/bdshared/grijalva/sfc/grijalva_maskUTM1km.dat"

   call InitMaskMap (nx,ny,maskfile,maskmap)

   do j = 1, ny
   do i = 1, nx
      write (6,'(2i4,i6)') i,j,maskmap(i,j)
   enddo
   enddo

    maskout = "/home/afonsoaraujo/work/bdshared/grijalva/sfc/grijalva_maskUTM1km.bin"
    open (unit=2,file=maskout,form='unformatted',access='direct',recl=(nx*ny*1))
    
    write(2,rec=1) maskmap



end program geramaskbin


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
