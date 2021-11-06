program extrai1d

   character(128)   ::   filein,fileout,fileoutaux,initime
   integer(2)       ::   lun,nx,ny
   parameter (nx=43,ny=30)
   integer          ::   nstations
   parameter (nstations=9)

   real             ::   array2d(nx,ny)
   real             ::   array1d(nstations,3,3)
   
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/calendario.inc'
   type (calendario)     ::   calend
   
   integer   ::   station(nstations)
   !data station / 999 /
   data station / 999,122,133,144,146,148,153,550,902 /
   
   integer   ::   pixel(nstations,2)
   integer   ::   kk, tstop
   parameter (tstop = 765)

!1km
!   data pixel /41, 5 /
   data pixel /41,33,17,36,26,17,10,38,30, 5,6,9,16,16,14,19,6,7/


! -------------------------------------------------------------
! Início do código
!
! -------------------------------------------------------------

   filein ='/media/DATA/home/afonso/workphd/dados/sgp97/forcantes/sgp97_lw1km-forcantes/sgp97_lw-pres.02190697'
   fileout   = '../data/forc/sgp97lw1d/sgp97_lw-pres1d-'

   initime   = '1997-06-19.02:00:00'

   call InitCalend (initime,calend)
   

   do kk = 1, tstop

      call readbinF (1,filein,nx,ny,array2d)

      call f_extrai1d (nx,ny,array2d,nstations,station,pixel,array1d)

      do k=1, nstations
	     fileoutaux = fileout
	     aux  = len_trim(fileoutaux)
         write(fileoutaux(aux+1:aux+19),'(i3,a,3i2.2,i4.4,a)') station(k),'.',calend%hora,  &
               calend%diam,calend%mes,calend%ano,'.bin'
	     write (6,'(a)') fileoutaux
		 call writebinF (1,fileoutaux,3,3,array1d)
      enddo
	  
      call nextstep (calend,nstations,station,filein,fileout)
   
   enddo
		
end program extrai1d


!---------------------------------------------------------------------------
! Inicializa calendario e estrutura de horarios e datas
!---------------------------------------------------------------------------
subroutine InitCalend(initime,calend)

   implicit none
   
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/calendario.inc'
   
   character(128)      ::   initime
   type (calendario)   ::   calend

   integer             ::   bis
   
!  Inicializa variáveis de calendario
   read( initime(1:4) , '(I4)') calend%ano
   read( initime(6:7) , '(I2)') calend%mes
   read( initime(9:10) , '(I2)') calend%diam

   call bissexto( calend%ano,bis )
   call datacor(bis, calend%mes, calend%diam, calend%diaj)

   write(6,'(4i4)') calend%ano, calend%mes,calend%diam,calend%diaj
      
   read( initime(12:13) , '(I2)') calend%hora
   read( initime(15:16) , '(I2)') calend%mint   

!   write(6,'(i4,1x,i2.2,1x,i2.2,1x,i3,1x,i2.2,1x,i2.2)')               &
!         (calend%ano-1900),calend%mes,calend%diam,calend%diaj,         &
!	       calend%hora,calend%mint   

end subroutine InitCalend


! -------------------------------------------------------------
! readbinf - lê arquivo binário contendo valores reais
!
! -------------------------------------------------------------
subroutine readbinF (lun,datafile,ncol,nlin,array)

   implicit none
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/grid.inc'

   character(128)      ::   datafile
   integer(2)          ::   lun,ncol,nlin,i,j,k
   real                ::   array(ncol,nlin)
   
!  variáveis de depuração - debug variables
   integer             ::   STAT_IOS
   integer             ::   STAT_ALLOC

   write(6,'(a)')   'Lendo arquivo binario com reais...: '
   write(6,'(a)')    datafile

   OPEN (UNIT=lun, FILE=datafile,ACCESS='DIRECT',FORM='unformatted',  &
          RECL=ncol*nlin*4, iostat=STAT_IOS)

   if (STAT_IOS .eq. 0) then
      write(6,'(a)') '  Arquivo binario aberto com sucesso! '
   else
      write(6,'(a)') '  Falha ao abrir arquivo! '
   endif

   read(lun,rec=1) array
   write(6,'(a)') '  Arquivo binario lido sucesso! '

   close (lun)

   return

end subroutine readbinF


! -------------------------------------------------------------
! f_extrai1d - extrai um pixel de dados de arquivos binarios
!
! -------------------------------------------------------------
subroutine f_extrai1d (ncol,nlin,array2d,nstations,station,pixel,array1d)

   implicit none

   integer             ::   ncol,nlin
   real                ::   array2d(ncol,nlin)
   integer             ::   nstations
   integer             ::   station(nstations)
   integer             ::   pixel(nstations,2)
   real                ::   array1d(nstations,3,3)

   integer             ::   i,j,k

   do k = 1, nstations
   do j = 1, 3
      do i = 1, 3
         array1d (k,i,j) = array2d( (pixel(k,1)+(-2+i)),(pixel(k,2)+(-2+j)) )
!		 write (6,'(6i4,2f12.4)') i,j,pixel(k,1),pixel(k,2),(pixel(k,1)+(-2+i)),  &
!		        (pixel(k,2)+(-2+j)),array1d(k,i,j),array2d( (pixel(k,1)+(-2+i)),(pixel(k,2)+(-2+j)) )
      enddo
   enddo
   enddo

   return

end subroutine f_extrai1d

! -------------------------------------------------------------
! writebinf - escreve arquivo binario contendo valores reais
!
! -------------------------------------------------------------
subroutine writebinF (lun,fileout,ncol,nlin,array)

   implicit none

   character(128)     ::   fileout
   integer(2)         ::   lun,ncol,nlin,i,j,k
	
   real               ::   array(ncol,nlin)
	   
!  variáveis de depuração - debug variables
   integer(2)          ::   STAT_IOS
   integer(2)          ::   STAT_ALLOC

   write(6,'(a)')   ' '
   write(6,'(a)')   'Escrevendo arquivo binF'
	
   open(unit=lun,file=fileout,access='direct',form='unformatted',  &
        recl=ncol*nlin,iostat=STAT_IOS)

   write(lun,rec=1) ((array(i,j),i=1,ncol),j=1,nlin)
!   ((wetdp2(i,j),i=1,nx),j=1,ny)
   
   close(lun)
   write(6,'(a)') '  Arquivo binario escrito sucesso! '

   return

end subroutine writebinF

! -------------------------------------------------------------
!  NextStep()
!
! -------------------------------------------------------------
subroutine nextstep(calend,nstations,station,filein,fileout)

   implicit none
   
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/calendario.inc'
   type (calendario)   ::   calend

   integer             ::   nstations
   integer             ::   station(nstations)
   character(128)      ::   filein
   character(128)      ::   fileout,fileoutaux
!
!#######################################################################
!
!  Variáveis locais - Local variables
!
!#######################################################################
!   
   integer   ::   i,j,k
   integer   ::   bis,aux,aux1,aux2
!-----------------------------------------------------------------------
! Prepara variáveis para o próximo passo de tempo

   call bissexto( calend%ano,bis )

   calend%hora = calend%hora + 1
   if (calend%hora .gt. 23) then
      calend%hora  = 0
      calend%diam = calend%diam + 1
      calend%diaj= calend%diaj + 1         
      endif
!      write(6,'(6i6)') ano,mes,diam,diaj,hora,NNDias(bis,mes)
      if (calend%diam .gt. NNDias(bis,calend%mes)) then
         calend%diam = 1
         calend%mes  = calend%mes + 1
      endif
      if (calend%mes .gt. 12) then
         calend%mes  = 1
         calend%ano  = calend%ano + 1
      endif

      aux1 = scan(filein,'.',.true.)
      aux2 = len_trim(filein)
      write(filein(aux1+1:aux2),'(4i2.2)') calend%hora,      &
         calend%diam,calend%mes,(calend%ano-1900)

      do k=1, nstations
	     fileoutaux = fileout
	     aux  = len_trim(fileoutaux)
         write(fileoutaux(aux+1:aux+19),'(i3,a,3i2.2,i4.4,a)') station(k),'.',calend%hora,  &
               calend%diam,calend%mes,calend%ano,'.bin'
	     write (6,'(a)') fileoutaux
      enddo
         
end subroutine nextstep


! -------------------------------------------------------------------------
!   --> Bissexto: Diz se um ano e' bissexto ou nao
!
!          Entrada:
!             Ano      -- Numero do ano de 0 em diante
!
!          Saida:
!             Bissexto -- Um valor logico ( TRUE ou FALSE )
! -------------------------------------------------------------------------
subroutine bissexto ( ano, bis )

   implicit none
   
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/calendario.inc'

   integer   ::   ano
   integer   ::   bis

   IF ( (MOD(ano,4)   .eq. 0) .and. (MOD(ano,100) .ne. 0) .or.        &
      (MOD(ano,400) .eq. 0) ) THEN
       bis = 2
   ELSE
       bis = 1
   ENDIF

   RETURN   

end subroutine bissexto


!-----------------------------------------------------------------------
!  --> DataCor: Dado o fato de o ano ser ou nao bissexto e um par mes e
!               dia, retorna um dia corrido entre 1 e 366
!
!         Entrada:
!            Bis    -- Se o ano e' ou nao bissexto
!            Mes    -- Mes do ano
!            Dia    -- Dia do mes
!
!         Saida:
!            DiAno  -- Dia do ano
!------------------------------------------------------------------------
subroutine datacor ( bis, mes, dia, diano)

   implicit none

   integer   ::   bis
   integer   ::   mes
   integer   ::   dia
   integer   ::   diano

   integer    ::  k
   
   include '/media/DATA/home/afonso/workposdoc/modelos/ephstc1.1.0/include/calendario.inc'

       
   diano = dia
   k = 1 ;
   do while  ( k .lt. mes )
      diano = diano + NNDias(bis, k)
      k = k + 1
   end do

   RETURN

END subroutine datacor
