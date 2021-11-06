!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                      PHYCST.INC                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'phycst.inc' for ARPS

!  This file initializes values of some physical constants
!  that remain unchanged throughout a model run
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/1991
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  2/24/97 (Jinxing Zong, Ming Xue and Yuhe Liu)
!  Added five pwr arrays and a common block TAOICECST for look-up
!  tables to replace fractional power calculations in Tao ice scheme.
!
!  02/15/2002 (J. Brotzge)
!  Added variables for air, water, soil, heat capacities
!  Added tau, tauinv, cnbeta, cpinv 
!
!-----------------------------------------------------------------------
!

!!DEC$ IF .not. DEFINED (PHYCST_INC)
!!DEC$ DEFINE PHYCST_INC

  REAL :: rd        ! Gas constant for dry air  (m**2/(s**2*K))
  PARAMETER( rd     = 287.0 )

  REAL :: rv        ! Gas constant for water vapor  (m**2/(s**2*K)).
  PARAMETER( rv     = 461.0 )

  REAL :: cv        ! Specific heat of dry air at constant volume
                    ! (m**2/(s**2*K)).
  PARAMETER( cv     = 717.0 )

  REAL :: cp        ! Specific heat of dry air at constant pressure
                    ! (m**2/(s**2*K)).
  PARAMETER( cp     = 1004.0 )

  REAL :: cpinv
  PARAMETER (cpinv = 1.0/cp )

  REAL :: cpdcv     ! cp/cv used in sound wave speed formula
  PARAMETER( cpdcv  = cp/cv )

  REAL :: rddcp
  PARAMETER( rddcp  = rd/cp )

  REAL :: rvdrd
  PARAMETER( rvdrd  = rv/rd )

  REAL :: rddrv
  PARAMETER( rddrv  = rd/rv )

  REAL :: lathv     ! Latent heat of vaporization at 0 degrees Celsius
                    ! (m**2/(s**2)).
  PARAMETER( lathv  = 2.500780E6 )

  REAL :: latheatv(151)
  REAL :: pwr1364(0:1001), pwr875 (0:1001)
  REAL :: pwr525 (0:1001), pwr2046(0:1001)
  COMMON /microcst/ latheatv, pwr1364, pwr525, pwr875, pwr2046

  REAL :: pwr81(0:151), pwr1666(0:10001)
  REAL :: pwr2(0:10001), pwr0625(0:10001), pwr15625(0:10001)
  COMMON /taoicecst/ pwr81, pwr1666, pwr2, pwr0625, pwr15625

! EMK 8 November 2002                     
  REAL :: pwr105(0:10001), pwr1596(0:10001), pwrlam195ratio(0:10001) ! EMK
  COMMON /ferriericecst/ pwr105, pwr1596, pwrlam195ratio ! EMK

  REAL :: lathf     ! Latent heat of fusion at 0 degrees Celsius
                    ! (m**2/(s**2)).
  PARAMETER( lathf  = 3.4E5 )

  REAL :: g         ! Acceleration due to gravity at the earth surface.
                    ! (m/(s**2))
  PARAMETER( g      = 9.8 )

  REAL :: p0        ! Surface reference pressure, is 100000 Pascal.
  PARAMETER( p0     = 1.0E5 )

  REAL :: rhow      ! Liquid water reference density (kg/m**3)
  PARAMETER( rhow = 1.0E3 )
  
!Ephst - Afonso Araujo - 30/12/2004
  real :: roa       ! Near sfc air density
  PARAMETER( roa = 1.2)

  REAL :: omega     ! The angular velocity of earth rotation (1/s)
  PARAMETER( omega=7.292E-5 )

  REAL :: solarc    ! Solar constant
  PARAMETER ( solarc = 1353.0 )

  REAL :: sbcst     ! Stefen-Boltzmann constant
  PARAMETER ( sbcst = 5.67E-8 )

  REAL :: emissg    ! Emissivity of the ground
  PARAMETER ( emissg = 0.995 )


  REAL :: emissa    ! Emissivity of the atmosphere
  PARAMETER ( emissa = 0.725 ) 

  REAL :: tau                        ! Seconds of a day = 24. * 3600.
  PARAMETER (tau = 86400.)           ! One day

  REAL :: tauinv                     ! Seconds of a day = 24. * 3600.
  PARAMETER (tauinv = 1.0/tau)       ! One day

  REAL :: cnbeta
  PARAMETER (cnbeta = 0.5)           ! Crank-Nicholson coefficient

  REAL :: cair
  PARAMETER (cair = 1004.0)           ! Volumetric heat capacity for air

  REAL :: cwater
  PARAMETER (cwater = 4200000.0)      ! Volumetric heat capacity for water

  REAL :: csoil
  PARAMETER (csoil =  1260000.0)      ! Volumetric heat capacity for soil

  REAL :: cice
  PARAMETER (cice = 2106000.0)        ! Volumetric heat capacity for ice

!!DEC$ ENDIF
