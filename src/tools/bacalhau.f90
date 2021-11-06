PROGRAM bacalhau

   implicit none
!   include 'grid.inc'

   integer             ::   lun,ncol,nlin,i,j,k
   parameter(ncol=572,nlin=468)
   
   real                ::   outarray(ncol,nlin)
   real                ::   inbuf(ncol)
   
   integer             ::   STAT_IOS
   integer             ::   STAT_ALLOC


   OPEN (UNIT=1,FILE="grijalva_maskUTM1km.bin",ACCESS='DIRECT',FORM='unformatted',  &
          RECL=ncol*1, iostat=STAT_IOS)

   if (STAT_IOS .eq. 0) then
      write(6,'(a)') ' Arquivo aberto com sucesso! '
   else
      write(6,'(a)') ' Falha ao abrir arquivo! '
   endif

   do j=1, nlin
      read (lun, REC = j) inbuf
      do i = 1, ncol
         outarray(i,j) = inbuf(i)
!         if (forc == 1) then
!            write(6,'(2i4,f8.3)')  i,j, outarray(i,j)
!         endif
      enddo
   enddo

   close (1)


end program bacalhau
