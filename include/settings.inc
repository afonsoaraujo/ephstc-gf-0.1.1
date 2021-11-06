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


!DEC$ IF .not. DEFINED (SETTINGS_INC)
!DEC$ DEFINE SETTINGS_INC

integer   ::   BUFSIZE
parameter (BUFSIZE = 255)

integer   ::   NDIRSMAX    ! número máximo de direções 
parameter (NDIRSMAX = 8)

integer   ::   NDIRS       ! número de direções consideradas
parameter (NDIRS = 4)

integer   ::   NX,NY,NZ
!PARAMETER( NX = 1417, NY = 1003, NZ = 1)
PARAMETER( NX = 425, NY = 301, NZ = 1)

integer   ::   NX_Forc, NY_Forc
PARAMETER(NX_Forc = 43, NY_Forc=30)

!integer   ::   DX,DY
!PARAMETER( DX = 1,  DY = 1 )

integer   ::   NSTypes
PARAMETER (NSTypes = 13)

integer   ::   NSLayers
PARAMETER (NSLayers = 3)

integer   ::   NVTypes
PARAMETER (NVTypes = 14)

integer   ::   NVLayers
PARAMETER (NVLayers = 3)

integer   ::   NFiles
PARAMETER (NFiles = 20)

integer   ::   neq
PARAMETER (neq = 6)

!DEC$ ENDIF
