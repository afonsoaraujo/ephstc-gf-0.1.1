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
subroutine NextStep(ano,mes,diam,diaj,hora,mint,dtforc,forcflopt,     &
                    tairfile,presfile,relhfile,wspdfile,rsifile,      &
		    rsrfile,precipfile,nx,ny,maskmap,topomap,         &
		    wtdepthmap,wtlevelmap,slopemap,aspectmap,fgradmap)

   implicit none
   
   integer          ::   it
   integer          ::   ano
   integer          ::   mes
   integer          ::   diam
   integer          ::   diaj
   integer          ::   hora
   integer          ::   mint
   integer          ::   dtforc
   integer          ::   forcflopt
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
!-----------------------------------------------------------------------
      call bissexto(ano,bis)
      hora = hora + (dtforc/(60*60))
      if (hora .gt. 23) then
         hora  = 0
         diam = diam + 1
         diaj= diaj+1         
      endif
      if (diam .gt. NNDias(bis,mes)) then
         diam = 1
         mes   = mes + 1
      endif
      if (mes .gt. 12) then
         mes   = 1
         ano   = ano + 1
      endif
!-----------------------------------------------------------------------
! Ajusta nomes dos arquivos dos forçantes
!-----------------------------------------------------------------------
      if (forcflopt == 1) then
         aux1 = scan(tairfile,'.',.true.)
         aux2 = len_trim(tairfile)
      
         write(6,'(a,3i4)') aux1,aux2,(ano-2000)
         write(tairfile(aux1+1:aux2),'(4i2.2)') hora,      &
               diam,mes,(ano-2000)
         
         aux1 = scan(presfile,'.',.true.)
         aux2 = len_trim(presfile)
         write(presfile(aux1+1:aux2),'(4i2.2)') hora,      &
               diam,mes,(ano-1900)

         aux1 = scan(relhfile,'.',.true.)
         aux2 = len_trim(relhfile)
         write(relhfile(aux1+1:aux2),'(4i2.2)') hora, &
               diam,mes,(ano-1900)
         
         aux1 = scan(wspdfile,'.',.true.)
         aux2 = len_trim(wspdfile)
         write(wspdfile(aux1+1:aux2),'(4i2.2)') hora, &
               diam,mes,(ano-1900)

         aux1 = scan(rsifile,'.',.true.)
         aux2 = len_trim(rsifile)
         write(rsifile(aux1+1:aux2),'(4i2.2)') hora,  &
               diam,mes,(ano-1900)         

         aux1 = scan(rsrfile,'.',.true.)
         aux2 = len_trim(rsrfile)
         write(rsrfile(aux1+1:aux2),'(4i2.2)') hora,  &
               diam,mes,(ano-1900)

         aux1 = scan(precipfile,'.',.true.)
         aux2 = len_trim(precipfile)
         write(precipfile(aux1+1:aux2),'(4i2.2)') hora,&
               diam,mes,(ano-1900)
      endif
      
      if (forcflopt == 2) then
         aux1 = scan(tairfile,'.',.true.)-8
         aux2 = len_trim(tairfile)
         if (ano > 2000) then
            write(tairfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(tairfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif
         
         aux1 = scan(presfile,'.',.true.)-8
         aux2 = len_trim(presfile)
         if (ano > 2000) then
            write(presfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(presfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif

         aux1 = scan(relhfile,'.',.true.)-8
         aux2 = len_trim(relhfile)
         if (ano > 2000) then
            write(relhfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
         
	    write(relhfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif
         
         aux1 = scan(wspdfile,'.',.true.)-8
         aux2 = len_trim(wspdfile)
         if (ano > 2000) then
            write(wspdfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(wspdfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif

         aux1 = scan(rsifile,'.',.true.)-8
         aux2 = len_trim(rsifile)
         if (ano > 2000) then
            write(rsifile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(rsifile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif

         aux1 = scan(rsrfile,'.',.true.)-8
         aux2 = len_trim(rsrfile)
         if (ano > 2000) then
            write(rsrfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(rsrfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif

         aux1 = scan(precipfile,'.',.true.)-8
         aux2 = len_trim(precipfile)
         if (ano > 2000) then
            write(precipfile(aux1:aux2),'(4i2.2,a)') (ano-2000), &
	          mes,diam,hora,'.r4'
         else
            write(precipfile(aux1:aux2),'(4i2.2,a)') (ano-1900), &
	          mes,diam,hora,'.r4'
         endif
      endif
      
!-----------------------------------------------------------------------
! if the flow gradient is based on the water table, recalculate the 
! water table gradients and the flow directions.
!-----------------------------------------------------------------------  
!   do i=1, nx
!   do j=1, ny
!      wtlevelmap(i,j) = topomap(i,j) - wtdepthmap(i,j)
!   enddo
!   enddo
!   call HeadSlopeAspect(nx,ny,maskmap,topomap,wtlevelmap,slopemap, &
!                        aspectmap,fgradmap)

end subroutine NextStep
