!
! SUMMARY:      data.h - header file with data structures
! USAGE:        Part of DHSVM
!
! AUTHOR:       Bart Nijssen
! ORG:          University of Washington, Department of Civil Engineering
! E-MAIL:       nijssen@u.washington.edu
! ORIG-DATE:    Apr-1996
! ADAPTATION:   Afonso Augusto M. de Araujo
! ADAP-DATE:    Dec-2003
!
!
! DESCRIPTION:  header file with data structures
! DESCRIP-END.
! FUNCTIONS:    
! COMMENTS:
! $Id: data.h,v 1.8 2002/11/19 17:08:00 nijssen Exp $     
!

!DEC$ IF .not. DEFINED (GLOBAL_INC)
!DEC$ DEFINE DATA_INC

character(128)   ::   runname
character(128)   ::   initime
character(128)   ::   endtime
integer          ::   dtbig
integer          ::   dtsfc
integer          ::   dtforc
integer          ::   dtroutsfc
integer          ::   tstart
integer          ::   tstop
integer          ::   nsfcst

common   /ephst_001/ runname,initime,endtime,dtbig,dtsfc,dtforc,      &
                     dtroutsfc,tstart,tstop,nsfcst
                     
real(4)          ::   DX
real(4)          ::   DY
real(4)          ::   Xorig
real(4)          ::   Yorig
real(4)          ::   OffSetX
real(4)          ::   OffSetY
real(4)          ::   DX_Forc
real(4)          ::   DY_Forc

common   /ephst_002/ DX,DY,Xorig,Yorig,OffSetX,OffSetY,DX_Forc,DY_Forc

integer          ::   forcflopt
character(128)   ::   tairforcfl
character(128)   ::   airpressfl
character(128)   ::   relhforcfl
character(128)   ::   wspdforcfl
character(128)   ::   rsiforcfl
character(128)   ::   rsrforcfl
character(128)   ::   precipforcfl

common   /ephst_003/ forcflopt,tairforcfl,airpressfl,relhforcfl,  &
                     wspdforcfl,rsiforcfl,rsrforcfl,precipforcfl
                     
character(128)   ::  stablefl
character(128)   ::  vtablefl

common   /ephst_004/ stablefl, vtablefl


character(128)   ::   maskfl
character(128)   ::   topofl
character(128)   ::   aspectfl
character(128)   ::   slopefl
character(128)   ::   soilfl
character(128)   ::   sdepthfl
character(128)   ::   vegfl
character(128)   ::   ndvifl

common   /ephst_005/ maskfl,topofl,aspectfl,slopefl,soilfl,sdepthfl,   &
                     vegfl,ndvifl
                     
character(128)   ::   wtdepthfl

common   /ephst_006/ wtdepthfl                     


integer          ::   slopeopt

common   /ephst_007/ slopeopt


integer          ::   ho

common   /ephst_008/ ho


real(4)          ::   sigma
real(4)          ::   emissa
real(4)          ::   emissg
real(4)          ::   roa
real(4)          ::   row
real(4)          ::   pi
integer(4)       ::   tau
real(4)          ::   L
real(4)          ::   Cp
real(4)          ::   Cv
real(4)          ::   Wl
real(4)          ::   Cdhlnd
real(4)          ::   Cdhwtr
real(4)          ::   Cdqlnd
real(4)          ::   Cdqwtr
real(4)          ::   C4lref
real(4)          ::   C4b
real             ::   d1
real             ::   d2
real(4)          ::   d3
real(4)          ::   dr

common   /ephst_009/ sigma,emissa,emissg,roa,row,pi,tau,L,Cp,Cv,Wl,   &
                     Cdhlnd,Cdhwtr,Cdqlnd,Cdqwtr,C4lref,C4b,d1,d2,d3, &
                     dr

integer          ::   ephstoutopt
character(128)   ::   ephstoutfl
character(128)   ::   ephstchkfl

common   /ephst_010/ ephstoutopt,ephstoutfl, ephstchkfl

!DEC$ ENDIF
