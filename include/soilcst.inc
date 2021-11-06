!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    SOILCST.INC                       ######
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
!  Include file 'soilcst.inc' for ARPS soil model.
!
!  This file defines the parameters and common blocks to be used
!  in the soil-vegetation energy and moisture model
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  09/25/1995
!
!  MODIFICATION HISTORY:
!
!  02/25/02 (J. Brotzge)
!  Added several new variables to specify soil moisure/temp
!    profiles as a function of soil type.
!    Data taken from Chen and Dudhia, MWR, 2001;
!    Cosby et al. (1984), Wat Res Res.
!    and Rawls et al (1983), J. of Hydraulic Eng.
!
!-----------------------------------------------------------------------
!

!!DEC$ IF .not. DEFINED (SOILCST_INC)
!!DEC$ DEFINE SOILCST_INC

  REAL :: bslope(13)         ! Slope of the retention curve
  DATA bslope/  4.05,  4.38,  4.90,  5.30,  5.39,  7.12,                &
                7.75,  8.52, 10.40, 10.40, 11.40, 11.55, 0.00/

!Arps
!  REAL :: cgsat(13)          ! Thermal coefficient for bare ground
!                             ! at saturation
!  DATA cgsat/ 3.222E-6, 3.057E-6, 3.560E-6, 4.418E-6,                   &
!              4.111E-6, 3.670E-6, 3.593E-6, 3.995E-6,                   &
!              3.058E-6, 3.729E-6, 3.600E-6, 1.00E-30, 8.0116E-6/

!Ephst - Afonso Araujo - 30/12/2004
  real :: cgsat(13)           ! Thermal coefficient for bare ground
                              ! at saturation
  data cgsat/ 3.222E-6, 3.057E-6, 3.560E-6, 4.418E-6,                 &
              4.111E-6, 3.670E-6, 3.593E-6, 3.995E-6,                 &
              3.058E-6, 3.729E-6, 3.600E-6, 1.00E-30, 8.0116E-6/

  REAL :: psisat(13)         ! Saturated matric potential (m)
  DATA psisat/ .121, .090, .218, .786, .478, .299,                    &
               .356, .630, .153, .490, .405, .036, .000/

  REAL :: kns(13)            ! Saturated hydraulic conductivity (m/s)
  DATA kns/ 1.07E-6, 1.41E-5, 5.23E-6, 2.81E-6, 3.38E-6, 4.45E-6,     &
              2.04E-6, 2.45E-6, 7.22E-6, 1.34E-6, 9.74E-7, 1.34E-6,   &
              0.0/
!Arps
!  INTEGER :: pwgeq(13)       ! Exponent in Wgeq formula
!  DATA pwgeq/  4,  4,  4,  6,  6,  6,                                   &
!               8, 10,  8, 10, 12,  0,  0/

!Ephst - Afonso Araujo - 30/12/2004
  integer :: pwgeq(13)       ! Exponent in Wgeq formula
  data pwgeq/  4,  4,  4,  6,  6,  6, 8, 10,  8, 10, 12,  0,  0/

!Arps
!  REAL :: awgeq(13)          ! Coefficient in Wgeq formula
!  DATA awgeq/ .387, .404, .219, .105, .148, .135,                       &
!              .127, .084, .139, .075, .083, .000, .000/

!Ephst - Afonso Araujo - 30/12/2004
  real :: awgeq(13)          ! Coefficient in Wgeq formula
  data awgeq/ .387, .404, .219, .105, .148, .135,                       &
              .127, .084, .139, .075, .083, .000, .000/

!Arps
!  REAL :: c1sat(13)          ! Value of c1 at saturation (when d1=10cm)
!  DATA c1sat/ .082, .098, .132, .153, .191, .213,                    &
!              .385, .227, .421, .375, .342, .000, .000/

!Ephst - Afonso Araujo - 30/12/2004
  real :: c1sat(13)          ! Value of c1 at saturation (when d1=10cm)
  data c1sat/ .082, .098, .132, .153, .191, .213,                     &
              .385, .227, .421, .375, .342, .000, .000/

! Arps
!  REAL :: c2ref(13)          ! Value of c2 at wetdp = 0.5*wsat
!  DATA c2ref/ 3.9, 3.7, 1.8,  .8,  .8,  .8,                          &
!               .4,  .6,  .3,  .3,  .3,  .0,  .0/

!Ephst - Afonso Araujo - 30/12/2004
  real :: c2ref(13)
  data c2ref/  3.9, 3.7, 1.8,  .8,  .8,  .8, .4,  .6,  .3,            &
                .3,  .3,  .0,  .0/
		
!Ephst - Afonso Araujo - 30/12/2004
  real :: c1wlt(13)
  data c1wlt/ 15.00, 15.00, 15.00, 15.00, 15.00, 15.00, 15.00,        &
              15.00, 15.00, 15.00, 15.00, 15.00, 15.00/
	      
!Ephst - Afonso Araujo - 30/12/2004
  real :: C4lref(13)
  data c4lref/ 25096., 19485., 4623.,  736., 518., 3835., 546.,       &
                2257., 18525., 3059., 6089.,   0.,    0./
		
!Ephst - Afonso Araujo - 30/12/2004
  real :: C4b(13)
  data C4b/ 5.800,  6.000,  6.400, 6.900, 6.400, 7.800, 8.300, 9.300, &		
           10.600, 11.100, 11.500, 0.000, 0.000/

!Ephst - Afonso Araujo - 30/12/2004
! Boone et al. (1999)
   real :: c4(13)
   data c4/ 0.0185,  0.0470, 0.0166, 0.0023, 0.0096, 0.0250, 0.0040,  &
            0.0040, 0.0158, 0.0022, 0.0074, 0.0000, 0.0000/
           

!  REAL :: wsat(13)           ! Saturated volumetric moisture content.(ETA)
!  DATA wsat / .421, .421, .434, .434, .439, .404,                     &
!              .464, .465, .406, .406, .468, 1.0E-20, 1.00/

!Ephst - Afonso Araujo - 30/12/2004
   real :: wsat(13)           ! Saturated volumetric moisture content.
   data wsat/ 0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, &
              0.426, 0.482, 0.482, 1.0E-30, 1.00/

! Arps
!  REAL :: wfc(13)            ! Field capacity moisture. (ETA) 
!  DATA wfc  / .283, .283, .312, .312, .329, .315,                    &
!             .387, .382, .338, .338, .412, 1.0E-25, 1.00/

!Ephst - Afonso Araujo - 30/12/2004
  real :: wfc(13)             ! Field capacity moisture.
  data wfc/ 0.135, 0.150, 0.195, 0.255, 0.240, 0.255, 0.322, 0.325,   &
            0.310, 0.370, 0.367, 1.0E-25, 1.00/

! Arps
!  REAL :: wwlt(13)           ! Wilting volumetric moisture content. (ETA)
!  DATA wwlt / .029, .029, .047, .047, .066, .069,                       &
!              .119, .103, .020, .020, .139, 1.0E-30, 1.00/

!Ephst - Afonso Araujo - 30/12/2004
   real :: wwlt(13)           ! Wilting volumetric moisture content.
   data wwlt/ 0.068, 0.075, 0.114, 0.179, 0.155, 0.175, 0.218, 0.250, &
              0.219, 0.283, 0.286, 1.0E-30, 1.00/

!  REAL :: wsat(13)           ! Saturated volumetric moisture content.
!  DATA wsat / .395, .410, .435, .485, .451, .420,                       &
!              .477, .476, .426, .492, .482, 1.0E-20, 1.00/
!
!  REAL :: wfc(13)            ! Field capacity moisture.
!  DATA wfc  / .135, .150, .195, .255, .240, .255,                       &
!             .322, .325, .310, .370, .367, 1.0E-25, 1.00/
!
!  REAL :: wwlt(13)           ! Wilting volumetric moisture content.
!  DATA wwlt / .068, .075, .114, .179, .155, .175,                       &
!              .218, .250, .219, .283, .286, 1.0E-30, 1.00/

  REAL :: rgl(14)            ! Species-dependent threshold value of
                             ! incoming solar radiation
  DATA rgl  / 100.0, 100.0, 100.0,  30.0,  30.0,  30.0,  30.0,          &
               30.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0/

  REAL :: rsmin(14)          ! Minimum of surface resistance
  DATA rsmin/ 999.0, 500.0,  40.0,  80.0,  80.0, 100.0, 150.0,          &
              100.0, 999.0,  40.0,   5.0,  80.0, 300.0, 1.0E-30/

  REAL :: alpha(7)    ! Polynomial fit used for estimating des/dT
  DATA alpha/ 610.7799961, 44.36518521,1.428945805,            &
           2.650648471E-2, 3.031240396E-4,2.034080948E-6,      &
           6.136820929E-9/

  REAL :: quartz(13)         !Quartz content
  DATA quartz/ 0.92, 0.82, 0.60, 0.25, 0.40, 0.60, 0.10,       &
               0.35, 0.52, 0.10, 0.25, 0.00, 0.00/
!Arps
!  REAL :: porosity(13)       ! Effective porosity (cm3/cm3)
!  DATA porosity/ 0.417, 0.421, 0.434, 0.486, 0.439, 0.404,     &
!                 0.464, 0.465, 0.406, 0.423, 0.468, 0.000,     &
!                 0.000/

!Ephst - Afonso Araujo - 30/12/2004
   real :: porosity(13)       ! Effective porosity (cm3/cm3)
   data porosity/ 0.430, 0.420, 0.400, 0.460, 0.430, 0.390,     &
                  0.480, 0.460, 0.410, 0.490, 0.470, 1.000,     &
		  1.000/

!Ephst - Afonso Araujo - 30/12/2004
   real :: press(13)          ! Suction pressure
   data press/ 0.121, 0.090, 0.218, 0.786, 0.478, 0.299, 0.356, &
               0.630, 0.153, 0.490, 0.405, 1.000, 1.000/
		  
		  
!Ephst - Afonso Araujo - 30/12/2004
   real :: poredist(13)       ! Pore distribution
   data poredist/ 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,     &
                  1.000, 1.000, 1.000, 1.000, 1.000, 1.000,     &
		  1.000/


  REAL :: rootzone(14)       ! Effective root zone depth (m) 
  DATA rootzone/ 1.0, 0.40, 1.0, 1.0, 1.0, 2.0, 2.0,           &
                 2.0, 0.40, 1.0, 1.0, 1.0, 1.0, 2.0/
                 
!Ephst - Afonso Araujo - 30/12/2004
  real :: albdveg(14)
!  data albdveg /0.55,0.45,0.40,0.35,0.37,0.28,0.26, 0.21,      &
!                  0.55, 0.40, 0.40, 0.40, 0.33, 0.17/

  data albdveg /0.47,0.19,0.16,0.16,0.12,0.13,0.12, 0.13,      &
                0.47, 0.18, 0.18, 0.18, 0.16, 0.12/

!  data albdveg /0.45,0.31,0.27,0.25,0.26,0.18,0.16, 0.11,      &
!                  0.45, 0.19, 0.19, 0.25, 0.23, 0.19/

  REAL :: snup(14)           ! Maximum snow depth threshold at max albedo
  DATA snup/ 0.040, 0.025, 0.040, 0.040, 0.040, 0.080, 0.080,  &
             0.080, 0.025, 0.040, 0.040, 0.040, 0.025, 0.000/ 

  REAL :: hsf2(14)           ! Variable used for F2 in Rc estimate 
  DATA hsf2/ 0.00, 42.0, 36.35, 42.0, 54.53, 54.53, 41.69,     &
             60.0, 999.0, 36.35, 60.0, 42.0, 200.0, 0.0/   

  REAL :: cgv                ! Thermal coefficient for vegetation
!  DATA cgv/1.5E-05/
  DATA cgv/2.4E-05/
!  DATA cgv/1.0E-04/

  REAL :: dirfc              ! Directional factor of Rayleigh
                             ! scattering
                             ! and absorption transmission function,
                             ! above cloud: dirfc = dirf (not for sfc.)
  DATA dirfc/1.6666667/   ! below & within cloud: dirfc = 5/3

  REAL :: d1                 ! Top layer depth of soil column: 10 cm
  DATA d1/ 0.10/             ! Chaning d1 requires changing C1sat. The
                             ! latter is propotional to the inverse of d1

  REAL :: d2                 ! Deep layer depth of soil column: 1 m
  DATA d2/ 0.6/              ! Deep layer depth of soil column: 1 m
  
!Ephst - Afonso Araujo - 30/12/2004
  real :: d3                 !
  data  d3 /1.1/
  
!Ephst - Afonso Araujo - 30/12/2004
  real :: dr                 !
  data  dr /0.5/  

  REAL :: rsmax              ! Maximum surface resistance in s/m
  DATA rsmax/5000./       ! Maximum surface resistence

  REAL :: wetsml             ! NP, Eq. 18, small value to limit
  DATA wetsml/1.0E-30/    ! Small value to limit c2wg at saturation
                          ! w2 = wsat

!  DATA roufns/ 0.002, 0.020, 0.010, 0.100, 0.150, 0.750, 1.000,      &
!               2.000, 0.005, 0.010, 0.020, 0.060, 0.040, 0.001/

!  DATA veg  / 0.10, 0.10, 0.60, 0.40, 0.40, 0.90, 0.99,              &
!                 0.99, 0.01, 0.30, 0.99, 0.40, 0.20, 0.00/

!Ephst - Afonso Araujo - 30/12/2004
  real  :: vegf(14)
  DATA vegf  / 0.10, 0.10, 0.60, 0.40, 0.40, 0.90, 0.99,               &
               0.99, 0.01, 0.30, 0.99, 0.40, 0.20, 0.00/

  REAL :: snow_albedo        ! Albedo of snow covered ground (was .65, now .6)
  REAL :: snowdepth_crit     ! Depth (m) at which snow cover is considered to
                             ! to be 100%.
  REAL :: snowflxfac         ! Factor by which to modify snow fluxes
  REAL :: cg_snow            ! heat capacity of snow
  PARAMETER (snow_albedo=0.60,snowdepth_crit=0.025,                   &
             cg_snow=1.0/60610.0,snowflxfac=1.0)

!Ephst - Afonso Araujo - 30/12/2004
  real :: ks(13)              ! Hydraulic conductivy  [m/s]
  data ks/ 0.0001760, 0.0001563, 0.0000341, 0.0000072, 0.0000070,     &
           0.0000063, 0.0000017, 0.0000025, 0.0000022, 0.0000010,     &
           0.0000013, 0.0000000, 0.0000000/
           
!Ephst - Afonso Araujo - 30/12/2004
!  real :: ks(13)              ! Hydraulic conductivy  [m/h]
!  data ks/ 0.6336000, 0.5508000, 0.1227600, 0.0259200, 0.0025200,     &
!           0.0022680, 0.0061200, 0.0090000, 0.0079200, 0.0036000,     &
!           0.0046800, 0.0000000, 0.0000000/           

!  data ks/ 1.490, 1.520, 1.570, 1.420, 1.490, 1.600, 1.380, 1.430,    &
!           1.570, 1.350, 1.390, 0.000, 0.000/
	   
!Ephst - Afonso Araujo - 30/12/2004
  real :: kslat(13)           ! Lateral hydraulic conductivy [m/s]
  data kslat/ 0.010560, 0.002990, 0.001440, 0.001089, 0.000542,       &
              0.000660, 0.001256, 0.000487, 0.000327, 0.000811,       &
              0.000874, 0.000000, 0.000000/

!  real :: kslat(13)           ! Lateral hydraulic conductivy [m/h]
!  data kslat/  38.016000, 10.764000, 5.1840000, 3.9204000, 1.9512000,     &
!                2.376000,  4.521600, 1.7532000, 1.1772000, 2.9196000,     &
!                3.146400,  0.0000000, 0.0000000/
              
	      
!Ephst - Afonso Araujo - 30/12/2004
  real :: kslatexp(13)            ! Lateral hydraulic conductivity expoent
  data kslatexp/ 4.100,  3.990, 4.840,  3.790, 5.300, 8.660, 7.480,   &
                 8.020, 13.000, 9.760, 12.280, 1.000, 1.000/

!Ephst - Afonso Araujo - 30/12/2004
  real :: maxinfiltrate(13)       ! Maximum infiltration rate
  data maxinfiltrate/ 500.00, 400.00, 400.00, 350.00, 200.00,         &
                      200.00, 300.00, 200.00, 200.00, 250.00,         &
		                280.00,    0.00,  0.00/

!Ephst - Afonso Araujo - 30/12/2004
  integer :: nlayers(13)
  data nlayers/ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1/

!!DEC$ ENDIF
