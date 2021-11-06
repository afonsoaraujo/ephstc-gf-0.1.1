!
!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                      EPHST.INC                       ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'ephst.inc' for ARPS
!
!  This file contains most of the model run control parameters
!  and parameters used by various physics packages.
!
!  These parameters are allocated in named common blocks and are
!  accessible to subroutines that include this file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Afonso A. M. de Araujo
!  14/02/2005
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

   integer,allocatable   ::   maskmap(:,:)   
   real,   allocatable   ::   topomap(:,:)
   real,   allocatable   ::   aspectmap(:,:)
   real,   allocatable   ::   slopemap(:,:)
   integer,allocatable   ::   soilmap(:,:)
   real,   allocatable   ::   sdepthmap(:,:,:)
   integer,allocatable   ::   vegmap(:,:)
   real,   allocatable   ::   vegfrac(:,:)
   real,   allocatable   ::   ndvimap(:,:)
   real,   allocatable   ::   laimap(:,:)
   real,   allocatable   ::   wtdepthmap(:,:)
   real,   allocatable   ::   tair(:,:)
   real,   allocatable   ::   pres(:,:)
   real,   allocatable   ::   relh(:,:)
   real,   allocatable   ::   wspd(:,:)
   real,   allocatable   ::   rsi(:,:)
   real,   allocatable   ::   rsr(:,:)
   real,   allocatable   ::   precip(:,:)
   
   common /ephstc001/ maskmap,topomap,aspectmap,slopmap,soilmap,      &
                      sdepthmap,vegmap,vegfrac,ndvimap,laimap,        &
!                      wtdepthmap,tair,pres,relh,wspd,rsi,rsr,precip
   
   
   
   
