




!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE TRNCNTL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE trncntl(nx,ny, ternfn, x,y)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS terrain data description file for GrADS display
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  10/05/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    x        The x-coord. of the physical and computational grid.
!             Defined at u-point.
!    y        The y-coord. of the physical and computational grid.
!             Defined at v-point.
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny             ! Number of grid points in 3 directions

  CHARACTER (LEN=*) :: ternfn  ! Terrain data file name

  REAL :: x(nx)                ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y(ny)                ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.

  REAL :: temxy1(nx,ny)        ! Temporary array
  REAL :: temxy2(nx,ny)        ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: trnctlfl
  CHARACTER (LEN=15) :: chrstr, chr1

  INTEGER :: nchout0

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/

  CHARACTER (LEN=2) :: dtunit
  INTEGER :: ntm,tinc

  REAL :: lonmin,latmin,lonmax,latmax
  REAL :: xbgn,ybgn
  REAL :: xinc,yinc
  REAL :: lat11,lat12,lon11,lon12,lat21,lat22,lon21,lon22
  REAL :: latinc,loninc

  INTEGER :: ctllen, fnlen
  INTEGER :: i,k, is
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Open the GrADS data control file
!
!-----------------------------------------------------------------------
!
  fnlen = LEN( ternfn )
  CALL strlnth( ternfn, fnlen )

  ctllen = fnlen + 4
  trnctlfl(1:ctllen) = ternfn(1:fnlen)//'.ctl'

  CALL fnversn( trnctlfl, ctllen )
  CALL getunit (nchout0)
  WRITE (6,'(a)') 'The GrADS data control file is '                     &
                    //trnctlfl(1:ctllen)

  OPEN (nchout0, FILE = trnctlfl(1:ctllen), STATUS = 'unknown')

  xbgn = 0.5 * (x(1) + x(2))
  ybgn = 0.5 * (y(1) + y(2))

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  CALL xytoll(nx,ny,x,y,temxy1,temxy2)

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  CALL a3dmax0(temxy1,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               latmax,latmin)
  CALL a3dmax0(temxy2,1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,                    &
               lonmax,lonmin)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'latmin:latmax:latinc = ',                                   &
            latmin,':',latmax,':',latinc
  WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                                 &
           'lonmin:lonmax:loninc = ',                                   &
           lonmin,':',lonmax,':',loninc

  IF(month <= 0) month = 1
  WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                         &
      hour,':',minute,'Z',day,monnam(month),year

  ntm = 1
  tinc = 1
  dtunit = 'MN'

  WRITE (nchout0,'(a/a)')                                               &
      'TITLE   ARPS terrain data control for '                          &
      //runname(1:lfnkey),'*'

  WRITE (nchout0,'(a,a)')                                               &
      'DSET    ',ternfn(1:fnlen)

  WRITE (nchout0,'(a)')                                                 &
      'OPTIONS sequential big_endian'

  WRITE (nchout0,'(a,i10)')                                             &
      'FILEHEADER ',192        !!! The number depends on the file
                               !!! structure of ternfn. See iolib3d.f90

  WRITE (nchout0,'(a/a)')                                               &
      'UNDEF   -9.e+33','*'

  IF ( mapproj == 2 .OR. mapproj == -2) THEN
!    WRITE (nchout0,'(a)')                                               &
!        '* For lat-lon-lev display, umcomment the following 4 lines.'

    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
              trulat1,trulat2,trulon,xinc,yinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',lonmin,loninc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',latmin,latinc
  ELSE
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx, '  LINEAR  ',xbgn,xinc

    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny, '  LINEAR  ',ybgn,yinc
  END IF

  WRITE (nchout0,'(a,1x,i8,a,i1)')                                      &
      'ZDEF',1,'  LEVELS  ',0

  WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                           &
      'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

  WRITE (nchout0,'(a,1x,i3)')                                           &
      'VARS',1

  WRITE (nchout0,'(a)')                                                 &
      'trn      0   99   ARPS terrain (m)'

  WRITE (nchout0,'(a)')                                                 &
      'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

  RETURN
END SUBROUTINE trncntl

