!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                    MODULE NESTSTEP                   ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!module nextstep


! -------------------------------------------------------------
!  NextStep()
!
! -------------------------------------------------------------
subroutine NextStep(ano,mes,diam,diaj,hora,mint,tairfile,presfile,    &
                    relhfile,wspdfile,rsifile,rsrfile,precipfile,     &
                    nx,ny,maskmap,topomap,wtdepthmap,wtlevelmap,      &
                    slopemap,aspectmap,fgradmap)

   implicit none
   
   integer          ::   it
   integer          ::   ano
   integer          ::   mes
   integer          ::   diam
   integer          ::   diaj
   integer          ::   hora
   integer          ::   mint
   character(128)   ::   tairfile,presfile,relhfile,wspdfile,rsifile, &
                         rsrfile,precipfile
   integer          ::   nx,ny                         
   integer          ::   maskmap(nx,ny)
   real(4)          ::   topomap(nx,ny)
   real(4)          ::   wtdepthmap(nx,ny)
   real(4)          ::   wtlevelmap(nx,ny)
   real(4)          ::   slopemap(nx,ny)
   real(4)          ::   aspectmap(nx,ny)
   real(4)          ::   fgradmap(nx,ny)
   
   real(4)          ::   HeadSlopeAspect
!
!#######################################################################
!
!  Include:
!
!#######################################################################
!   
   include 'calendario.inc'
!
!#######################################################################
!
!  Variáveis locais - Local variables
!
!#######################################################################
!   
   integer   ::   i,j
   integer   ::   bis,aux1,aux2,ii,jj   
   
!-----------------------------------------------------------------------
! Prepara variáveis para o próximo passo de tempo

write(6,'(1x,a)') 'Im here 1 !'	       		    

      call bissexto(ano,bis)
      hora = hora + 1
      if (hora .gt. 23) then
         hora  = 0
         diam = diam + 1
         diaj= diaj+1         
      endif
!      write(6,'(6i6)') ano,mes,diam,diaj,hora,NNDias(bis,mes)
      if (diam .gt. NNDias(bis,mes)) then
         diam = 1
         mes   = mes + 1
      endif
      if (mes .gt. 12) then
         mes   = 1
         ano   = ano + 1
      endif
      
write(6,'(1x,a)') 'Im here 2 !'	       		    
!-----------------------------------------------------------------------
! Ajusta nomes dos arquivos dos forçantes
      aux1 = scan(tairfile,'.',.true.)
      aux2 = len_trim(tairfile)
      write(6,'(2i4)') aux1, aux2
      write(6,'(i4.4,3i2.2)') hora,diam,mes,(ano-1900)
      write(tairfile(aux1+1:aux2),'(i4.4,3i2.2)') hora,      &
         diam,mes,(ano-1900)
	 
write(6,'(1x,a)') 'Im here 3 !'	       		    	 
         
      aux1 = scan(presfile,'.',.true.)
      aux2 = len_trim(presfile)
      write(presfile(aux1+1:aux2),'(i4.4,3i2.2)') hora,      &
         diam,mes,(ano-1900)

      aux1 = scan(relhfile,'.',.true.)
      aux2 = len_trim(relhfile)
      write(relhfile(aux1+1:aux2),'(i4.4,3i2.2)') hora, &
         diam,mes,(ano-1900)
         
      aux1 = scan(wspdfile,'.',.true.)
      aux2 = len_trim(wspdfile)
      write(wspdfile(aux1+1:aux2),'(i4.4,3i2.2)') hora, &
         diam,mes,(ano-1900)

      aux1 = scan(rsifile,'.',.true.)
      aux2 = len_trim(rsifile)
      write(rsifile(aux1+1:aux2),'(i4.4,3i2.2)') hora,  &
         diam,mes,(ano-1900)         

      aux1 = scan(rsrfile,'.',.true.)
      aux2 = len_trim(rsrfile)
      write(rsrfile(aux1+1:aux2),'(i4.4,3i2.2)') hora,  &
         diam,mes,(ano-1900)

      aux1 = scan(precipfile,'.',.true.)
      aux2 = len_trim(precipfile)
      write(precipfile(aux1+1:aux2),'(i4.4,3i2.2)') hora,&
         diam,mes,(ano-1900)
	 
!      write(6,'(a,i4)') ' diaj=', diaj
!      write(6,'(a)') tairfile
!      write(6,'(a)') presfile
!      write(6,'(a)') relhfile
!      write(6,'(a)') wspdfile
!      write(6,'(a)') rsifile
!      write(6,'(a)') rsrfile
!      write(6,'(a)') precipfile

!-----------------------------------------------------------------------
! if the flow gradient is based on the water table, recalculate the 
! water table gradients and the flow directions.
  
!   do i=1, nx
!   do j=1, ny
!      wtlevelmap(i,j) = topomap(i,j) - wtdepthmap(i,j)
!   enddo
!   enddo
!   call HeadSlopeAspect(nx,ny,maskmap,topomap,wtlevelmap,slopemap, &
!                        aspectmap,fgradmap)

end subroutine NextStep
