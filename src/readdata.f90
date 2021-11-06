!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                   MODULE readdata                    ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module readdata

!declaração de variáveis 


!contains

! -------------------------------------------------------------
! readbinf - lê arquivo binário contendo valores reais
!
! -------------------------------------------------------------

subroutine readbinf (lun,forc,datafile,ncol,nlin,outarray)

   implicit none
   include 'grid.inc'

   character(128)      ::   datafile
   integer             ::   forc
   integer             ::   lun,ncol,nlin,i,j,k
   real, allocatable   ::   inbuf(:)
   real                ::   outarray(ncol,nlin)


   
   
!  variáveis de depuração - debug variables
   integer             ::   STAT_IOS
   integer             ::   STAT_ALLOC

   allocate(inbuf(NX_Forc), stat=STAT_ALLOC)
   if (STAT_ALLOC .eq. 0) then
!      write(6,'(a)') ' Variável alocada com sucesso2! '
   else
      write(6,'(a,i4)') ' Erro de alocação - ', STAT_ALLOC      
   endif
   
!   write(6,'(a)') datafile
		
   OPEN (UNIT=lun, FILE=datafile,ACCESS='DIRECT',FORM='unformatted',  &
          RECL=ncol*1, iostat=STAT_IOS)

   if (STAT_IOS .eq. 0) then
!      write(6,'(a)') ' Arquivo aberto com sucesso! '
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

   deallocate(inbuf)
   close (lun)

   return

end subroutine readbinf

!end module readdata
