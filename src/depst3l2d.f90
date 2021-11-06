!
!  ##################################################################
!  ##################################################################
!  ######                                                      ######
!  ######                   MODULE EPST3L2D                    ######
!  ######                                                      ######
!  ######                     Developed by                     ######
!  ######               Laboratorio de Hidrologia              ######
!  ######      Universidade Federal do Rio de Janeiro          ######
!  ######                      COPPE/UFRJ                      ######
!  ######                                                      ######
!  ##################################################################
!  ##################################################################
!
SUBROUTINE Depst3l2d(nx,ny,nz,soiltyp,vegtyp,lai,veg,                   &
                     tsfc,tdeep,wetsfc,wetdp2,wetdp3,wetcanp,           &
                     snowdpth,qvsfc,windsp,psfc,rhoa,precip,            &
                     tair,qvair,cdha,cdqa,radsw,                        &
                     rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,evaprr,   &
                     qvsat,qvsata,f34,                                  &
                     SatFlow,RunOffei,RunOffCnp,wtdepthmap,sdepthmap)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Predict the soil surface temperature and moisture contents by solving
!  the surface energy and moisture budget equtions:
!
!  1. the ground surface temperature, Ts -- tsfc
!  2. the deep ground temperature,    T2 -- tdeep
!  3. the surface soil moisture,      wg -- wetsfc
!  4. the deep soil moisture,         w2 -- wetdp
!  5. the canopy moisture,            wr -- wetcanp
!
!-----------------------------------------------------------------------
!
!  The equations are listed as follows.
!
!
!       d Ts                        2PI
!      ------ = Ct (Rn - H - LE) - ----- (Ts - T2)
!       d t                         Tau
!
!
!       d T2      1
!      ------ = ----- (Ts - T2)
!       d t      Tau
!
!
!       d Wg      C1                 C2
!      ------ = ------ (Pg - Eg) - ----- (Wg - Wgeq)
!       d t     ROw d1              Tau
!
!
!       d W2      1
!      ------ = ------ (Pg - Eg - Etr)
!       d t     ROw d2
!
!
!       d Wr
!      ------ = Veg P - Er
!       d t
!
!
!  where
!
!      Tau    -- 1 day in seconds = 86400 seconds
!      PI     -- number of PI = 3.141592654
!      ROw    -- Density of liquid water
!      d1     -- Top layer depth of soil column, 0.01 m
!      d2     -- Deep layer depth of soil column, 1 m
!      Veg    -- Vegetation fraction
!      Ct     -- Thermal capacity
!      Rn     -- Radiation flux, rnflx
!      H      -- Sensible heat flux, shflx
!      LE     -- Latent heat flux, lhflx = latent*(Eg + Ev)
!      Eg     -- Evaporation from ground
!      Ev     -- Evapotranspiration from vegetation, Ev = Etr + Er
!      Etr    -- Transpiration of the remaining part of the leaves
!      Er     -- Evaporation directly from the foliage covered by
!                intercepted water
!      P      -- Precipitation rates
!      Pg     -- Precipitation reaching the ground,
!                Pg = (1 - Veg) P
!      Wgeq   -- Surface volumetric moisture
!      C1, C2 -- Coefficients
!
!  For detailed information about the surface energy budget model,
!  see the articles in the reference list.
!
!  The second-order Rouge-Kutta time integration scheme is used,
!  which is described below.
!
!  Assume a equation in the form of
!
!      d X
!     ----- = F(X, t)
!      d t
!
!  In the forward scheme, we have
!
!      X(1) = X(0) + dt * F[X(0), t0]
!
!  We split one time step into two halves, dt2 = dt/2, and use the
!  forward scheme to calculate the first half step X(1/2).
!
!      X(1/2)  = X(0) + dt2 * F[X(0), t0]
!
!  Then we can calculate the Right Hand Side (RHS) of the equation at
!  the half step, F[X(1/2), t(1/2)]. Finally, we calculate the one
!  step prediction, X(t1), by use of the average of F[X(0), t0] and
!  F[X(1/2), t(1/2)].
!
!      X(1) = X(0) + dt * 0.5 * { F[X(0),t0] + F[X(1/2), t(1/2)] }
!  REFERENCES:
!
!  Jacquemin, B. and J. Noilhan, 1990: Sensitivity Study and
!       Validation of a Land Surface Parameterization Using the
!       HAPEX-MOBILHP Data Set, Boundary-Layer Meteorology, 52,
!       93-134, (JN).
!
!  Noilhan, J. and S. Planton, 1989: A Simple Parameterization of
!       Land Surface Processes for Meteorological Model, Mon. Wea.
!       Rev., 117, 536-549, (NP).
!
!  Pleim, J. E. and A. Xiu, 1993: Development and Testing of a
!       Land-Surface and PBL Model with Explicit Soil Moisture
!       Parameterization, Preprints, Conf. Hydroclimat., AMS, 45-51,
!       (PX).
!
!  Bougeault, P., et al., 1991: An Experiment with an Advanced
!       Surface Parameterization in a Mesobeta-Scale Model. Part I:
!       Implementation, Monthly Weather Review, 119, 2358-2373.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  11/16/93
!
!  MODIFICATION HISTORY:
!
!  10/30/94 (Y. Liu)
!  Fixed a bug reported by Richard Carpenter.
!
!  11/01/1994 (Y. Liu)
!  Subroutine name of soil model were changed from SFCEBM to SOILEBM
!  and the arguments passed were changed to 2-D arrays instead of 3-D
!  temporary arrays.
!
!  12/12/1994 (Y. Liu)
!  Fixed a bug for the final calcultaion of qvsfc.
!
!  12/14/1994 (Y. Liu)
!  Fixed a bug in phycst.inc for the water density value which was
!  previously mistakingly set to 1 kg/m**3. The correct value should
!  be 1000 kg/m**3. This bug largely influenced the integration of Wg
!  and W2.
!
!  12/23/1994 (Y. Liu)
!  Added the runoff calculation for Wr.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the soil model and
!  at the same time delete the table data array veg(14).
!
!  03/27/1995 (Yuhe Liu)
!  Changed the solor radiation used in the calculation of surface
!  resistence factor F1 from the one at the top of atmosphere to the
!  one at the surface.
!
!  Changed the formula of calculating the surface resistence factor
!  F3 to F3=1, instead of varying with qvsat(Tair) and qvair.
!
!  12/8/1998 (Donghai Wang and Vince Wong)
!  Added a new 2-D permanent array, snowcvr(nx,ny), for snow cover.
!  We just used a simple scheme to consider the snow cover process.
!
!  2000/01/10 (Gene Bassett)
!  Snow cover (0 or 1) changed to snow depth (snowdpth).  For simplicity
!  a fractional value for snow cover is not used (simply say the grid
!  point is completely covered with snow if snowdpth > snowdepth_crit,
!  otherwise no snow).
!
!  2000/02/04 (Gene Bassett, Yang Kun)
!  Fixed an error in tsoil integration (rhst2) causing tsoil to change
!  by only 50% of what it should.
!
!  2001/12/07 (Diandong Ren, Ming Xue) 
!  Re-structured the code, moved the calculations of the right hand
!  side terms of the soil-vegetation model into subroutine depst3l2d_frc.
!
!  Soil seasonal temperature trend to be added.
!
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!  2002/06/9 (Ming Xue)
!  Revoked some modifications to the soil moisture related caculations
!  that Diandong Ren put in since IHOP_2 - the mods need more testing.
!
!  2002/12/13 (Jerry Brotzge) 
!  Updated code to match recommendations by Pleim and Xiu (1995) and 
!  Xiu and Pleim (2001 - JAM).  
!
!-----------------------------------------------------------------------
! 
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (sfc/top)
!
!    soiltyp  Soil type at the horizontal grid points
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    windsp   Wind speed just above the surface (m/s)
!    psfc     Surface pressure (Pascal)
!    rhoa     Near sfc air density
!    precip   Precipitation flux reaching the surface (m/s)
!
!    cdha     Array for cdh, surface drag coefficient for heat
!    cdqa     Array for cdq, surface drag coefficient for moisture
!
!    pres     3-dimensional pressure
!    temp     3-dimensional temperature
!    qv       3-dimensional specific humidity
!
!  INPUT/OUTPUT: 
!
!    tsfc     Temperature at ground surface (K)
!    tdeep    Deep soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Canopy water amount
!
!  OUTPUT:
!
!    qvsfc    Effective S. H. at sfc.
!
!  Local automatic work arrays for storing the right hand forcing terms 
!  of the soil model equations and for storing intermediate values of the
!  soil state variables. 
!
!    frc_tsfc   Temporary array
!    frc_tdeep  Temporary array
!    frc_wsfc   Temporary array
!    frc_wdp    Temporary array
!    frc_wcnp   Temporary array
!
!    tsfcn      Temporary array
!    tdeepn     Temporary array
!    wgn        Temporary array
!    w2n        Temporary array
!    wrn        Temporary array
!
!-----------------------------------------------------------------------
!
   IMPLICIT NONE

   INTEGER :: nx,ny,nz

   REAL :: windsp(nx,ny)         ! Wind speed
   REAL :: psfc  (nx,ny)         ! Surface pressure
   REAL :: rhoa  (nx,ny)         ! Air density near the surface
   REAL :: precip(nx,ny)         ! Precipitation rate at the surface

   INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
   INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

   REAL :: lai    (nx,ny)        ! Leaf Area Index
   REAL :: veg    (nx,ny)        ! Vegetation fraction

   REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
   REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
   REAL :: wetsfc (nx,ny)        ! Surface soil moisture
   REAL :: wetdp2 (nx,ny)        ! Deep soil moisture of layer 2
   REAL :: wetdp3 (nx,ny)        ! Deep soil moisture of layer 3
   REAL :: wetcanp(nx,ny)        ! Canopy water amount
   REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface 

   REAL :: snowdpth(nx,ny)       ! Snow depth (m)

   REAL :: tair   (nx,ny)        ! Surface air temperature (K)
   REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
   REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
   REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture  

   REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
   REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
   REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
   REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
   REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
   REAL :: ct     (nx,ny)        ! Soil thermal coefficient

   REAL :: evaprg (nx,ny)        ! Evaporation
   REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
   REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
   REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
   REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
   REAL :: qvsat  (nx,ny)        !
   
   real :: SatFlow(nx,ny)
   real :: RunOffei(nx,ny)
   real :: RunOffCnp(nx,ny)
   real :: wtdepthmap(nx,ny)
   real :: sdepthmap(nx,ny)

   REAL,allocatable :: frc_tsfc(:,:)       ! Right hand side forcing for tsfc eq.
   REAL,allocatable :: frc_tdeep(:,:)      ! Right hand side forcing for tdeep eq.
   REAL,allocatable :: frc_wsfc(:,:)       ! Right hand side forcing for wetsfc eq.
   REAL,allocatable :: frc_wdp2(:,:)       ! Right hand side forcing for wetsp eq.
   REAL,allocatable :: frc_wdp3(:,:)       ! Right hand side forcing for wetsp eq.
   REAL,allocatable :: frc_wcnp(:,:)       ! Right hand side forcing for wetcanp eq.

   REAL,allocatable :: tsfcn(:,:)          ! Temporary array, tsn
   REAL,allocatable :: tdeepn(:,:)         ! Temporary array, t2n
   REAL,allocatable :: wgn(:,:)            ! Temporary array, wetsfcNEW
   REAL,allocatable :: w2n(:,:)            ! Temporary array, w2n
   REAL,allocatable :: w3n(:,:)            ! Temporary array, w3n
   REAL,allocatable :: wrn(:,:)            ! Temporary array, wrn
   REAL             :: relief              ! Difference between seasonal average skin and deep soil 

!-----------------------------------------------------------------------
!  Include files: globcst.inc and phycst.inc
!-----------------------------------------------------------------------
!
!  Parameters and variables are defined in globcst.inc:
!
!    dtsfc      Surface model time step
!    nsfcst     # of surface model time steps
!
!    moist      Moist flag
!
!    year       Reference year
!    month      Reference month
!    day        Reference day
!    jday       Reference Julian day
!    hour       Hour of reference time
!    minute     Minute of reference time
!    second     Second of reference time
!
!    latitud    Latitude at the domain center
!    longitud   Longitude at the domain center
!
!    curtim     Current model time
!    dtbig      Length of big time step
!
!    bslope     Slope of the retention curve
!    cgsat      Soil thermal coefficient for bare ground at saturation
!    cgv        Soil thermal coef. for totally shielded ground by veg.
!    pwgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    awgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    c1sat      Value of C1 at saturation. NP, Tab. 2
!    c2ref      Value of C2 for W2 = .5 * Wsat. NP, Tab. 2
!    wsat       Saturated volumetric moisture content. JN, Tab. 1
!    wfc        Field capacity moisture. JN, Tab. 1
!    wwlt       Wilting volumetric moisture content. JN, Tab. 1
!
!  Parameters and variables are defined in phycst.inc:
!
!    solarc     Solar constant (W/m**2)
!    emissg     Emissivity of the ground
!    emissa     Emissivity of the atmosphere
!    sbcst      Stefen-Boltzmann constant
!
!    rhow       Liquid water reference density (kg/m**3)
!    rd         Gas constant for dry air (kg/(m s**2))
!    cp         Gas heat capacity at constant pressure
!    cv         Gas heat capacity at constant volume
!
!-----------------------------------------------------------------------
!
   INCLUDE 'globcst.inc'
   INCLUDE 'phycst.inc'
   INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
   REAL :: pi                         ! Pi
   PARAMETER (pi = 3.141592654)

!   REAL :: tau                        ! Seconds of a day = 24. * 3600.
!   PARAMETER (tau = 86400.)

   REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.

   REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun

   REAL :: cg          ! Soil thermal coefficient for bare ground

   REAL,allocatable :: rhsts(:,:)       ! Right hand side of Eq. for Ts at current time
   REAL,allocatable :: rhst2(:,:)       ! Right hand side of Eq. for T2 at current time
   REAL,allocatable :: rhswg(:,:)       ! Right hand side of Eq. for Wg at current time
   REAL,allocatable :: rhsw2(:,:)       ! Right hand side of Eq. for W2 at current time
   REAL,allocatable :: rhsw3(:,:)       ! Right hand side of Eq. for W3 at current time
   REAL,allocatable :: rhswr(:,:)       ! Right hand side of Eq. for Wr at current time
   REAL,allocatable :: wrmax(:,:)       ! Maximum value for canopy moisture, wetcanp
   REAL :: c1wg               ! Coefficient in the surface moisture Eq. of Wg
   REAL :: c2wg               ! Coefficient in the surface moisture Eq. of Wg 

   REAL :: wgeq               ! Surface moisture when gravity balances the capillarity

   REAL :: wr2max             ! Tendency to reach the maximum wrmax
   REAL :: pnet               ! Residual of precip. and evap.
   REAL :: vegp               ! Precip. intercepted by vegetation

   INTEGER :: i,j,k,it
   integer :: istatus   

   REAL :: tema, temb,EBCal

   LOGICAL :: firstcall        ! First call flag of this subroutine

!   SAVE firstcall, log100, dtsfc2, tauinv
   SAVE firstcall, log100, dtsfc2
   DATA firstcall/.true./

   INTEGER :: jday_min         ! offset value from Jan 01. 
!---------------------------------------------------------------------  
! Alocação dinâmica de memória
!---------------------------------------------------------------------
  if (.not. allocated(frc_tsfc)) ALLOCATE(frc_tsfc    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_tsfc")
  frc_tsfc = 0

  if (.not. allocated(frc_tdeep)) ALLOCATE(frc_tdeep    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_tdeep")
  frc_tdeep = 0

  if (.not. allocated(frc_wsfc)) ALLOCATE(frc_wsfc    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_wsfc")
  frc_wsfc = 0
  
  if (.not. allocated(frc_wdp2)) ALLOCATE(frc_wdp2    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_wdp2")
  frc_wdp2 = 0

  if (.not. allocated(frc_wdp3)) ALLOCATE(frc_wdp3    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_wdp3")
  frc_wdp3 = 0

  if (.not. allocated(frc_wcnp)) ALLOCATE(frc_wcnp    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:frc_wcnp")
  frc_wcnp = 0
  
  if (.not. allocated(tsfcn)) ALLOCATE(tsfcn    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:tsfcn")
  tsfcn = 0
  
  if (.not. allocated(tdeepn)) ALLOCATE(tdeepn    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:tdeepn")
  tdeepn = 0
  
  if (.not. allocated(wgn)) ALLOCATE(wgn    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:wgn")
  wgn = 0
  
  if (.not. allocated(w2n)) ALLOCATE(w2n    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:w2n")
  w2n = 0
  
  if (.not. allocated(w3n)) ALLOCATE(w3n    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:w3n")
  w3n = 0
  
  if (.not. allocated(wrn)) ALLOCATE(wrn    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:wrn")
  wrn = 0
  
  if (.not. allocated(rhsts)) ALLOCATE(rhsts    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhsts")
  rhsts = 0
  
  if (.not. allocated(rhst2)) ALLOCATE(rhst2    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhst2")
  rhst2 = 0
  
  if (.not. allocated(rhswg)) ALLOCATE(rhswg    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhswg")
  rhswg = 0
  
  if (.not. allocated(rhsw2)) ALLOCATE(rhsw2    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhsw2")
  rhsw2 = 0
  
  if (.not. allocated(rhsw3)) ALLOCATE(rhsw3    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhsw3")
  rhsw3 = 0
  
  if (.not. allocated(rhswr)) ALLOCATE(rhswr    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:rhswr")
  rhswr = 0
  
  if (.not. allocated(wrmax)) ALLOCATE(wrmax    (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "ephst:wrmax")
  wrmax = 0  
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   IF (firstcall) THEN
      log100    = ALOG(100.0)
      dtsfc2    = dtsfc/2.0
      firstcall = .false.
   END IF
!====================================================================== 
! Modelo de solo
!======================================================================
   IF ( moist /= 0 ) THEN
!----------------------------------------------------------------------
!   Calculate saturated specific humidity near the surface, qvsata.
!----------------------------------------------------------------------
      CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tair,qvsata)

      DO j = 1, ny-1
      DO i = 1, nx-1
!      f34(i,j) = MAX( 0., 1.0 - 0.0016 * (298.0-tair(i,j))**2 )

         IF (tair(i,j) >= 302.15) THEN                    !(XP2001) 
            f34(i,j) = 1.0 / (1.0 + EXP( 1.18 * (tair(i,j)-314.00)))
         ELSE
            f34(i,j) = 1.0 / (1.0 + EXP( -0.41* (tair(i,j)-282.05)))
         ENDIF 
!         if (i == sx .and. j == sy) then
!            write (6,'(2f18.10)') tair(i,j), f34(i,j)
!         endif	 

      END DO
      END DO
   END IF

   jday_min = 61       ! may change with latitude
   relief=tsoil_offset_amplitude*sin((jday-jday_min)/365.0*2.0*PI)

!=======================================================================
!  Start time integration loop
!=======================================================================
   DO it = 1, nsfcst               
  
      if (epelopt == 1 .or. epelopt == 2) then
         do j = 1, ny
         do i = 1, nx
            wetdp2(i,j) = wetdp2(i,j) - (SatFlow(i,j)/(d2 *      &
	                    porosity(soiltyp(i,j)) ))

            wetdp3(i,j) = wetdp3(i,j) - (SatFlow(i,j)/((d3-d2) *      &
	                    porosity(soiltyp(i,j)) ))
!            if (i == sx .and. j == sy) then
!               write (6,'(3f18.10)') wetdp2(i,j),SatFlow(i,j),wetdp3(i,j)
!            endif
         enddo
         enddo
      endif

      CALL depst3l2d_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,             &
                         tsfc,tdeep,wetsfc,wetdp2,wetdp3,wetcanp,     &
                         snowdpth,qvsfc,windsp,psfc,rhoa,precip,      &
                         tair,qvair,cdha,cdqa,radsw,                  &
                         rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,    &
                         evaprr,qvsat,qvsata,f34,                     &
                         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp2,        &
                         frc_wdp3,frc_wcnp,wrmax,relief,              &
                         SatFlow,RunOffei,RunOffCnp,wtdepthmap,       &
                         sdepthmap)  
                        
      DO j = 1, ny-1
      DO i = 1, nx-1
         tsfcn (i,j) = tsfc(i,j) + dtsfc2 * frc_tsfc(i,j)
         tdeepn(i,j) = tdeep(i,j)+ dtsfc2 * frc_tdeep(i,j)
      END DO
      END DO

      IF ( moist /= 0 ) THEN    

      DO j = 1, ny-1
      DO i = 1, nx-1
         wgn(i,j) = wetsfc(i,j) + dtsfc2 * frc_wsfc(i,j)
         wgn(i,j) = MAX(wgn(i,j), 0.0 )
         wgn(i,j) = MIN(wgn(i,j), wsat(soiltyp(i,j)) )

         w2n(i,j) = wetdp2(i,j) + dtsfc2 * frc_wdp2(i,j)
         w2n(i,j) = MAX( w2n(i,j), 0.0 )
         w2n(i,j) = MIN(w2n(i,j), wsat(soiltyp(i,j)) )
         
         w3n(i,j) = wetdp3(i,j) + dtsfc2 * frc_wdp3(i,j)
         w3n(i,j) = MAX( w3n(i,j), 0.0 )
         w3n(i,j) = MIN(w3n(i,j), wsat(soiltyp(i,j)) )         

         wrn(i,j) = wetcanp(i,j) + dtsfc2 * frc_wcnp(i,j)
         wrn(i,j) = MAX(wrn(i,j), 0.0 )
         wrn(i,j) = MIN(wrn(i,j), wrmax(i,j) )
      END DO
      END DO

      ELSE

      DO j = 1, ny-1
      DO i = 1,nx-1
         w2n(i,j) = wetdp2(i,j)
         w3n(i,j) = wetdp3(i,j)
      END DO
      END DO

      END IF
      
      DO j = 1, ny-1
      DO i = 1, nx-1
         IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
            rhsts(i,j) = 0.0
            rhst2(i,j) = 0.0
         ELSE
            rhsts(i,j)=frc_tsfc(i,j)
            rhst2(i,j)=frc_tdeep(i,j)
         END IF
      END DO
      END DO

      IF ( moist /= 0 ) THEN    

         DO j = 1, ny-1
         DO i = 1, nx-1
            rhswg(i,j)=frc_wsfc(i,j)
            rhsw2(i,j)=frc_wdp2(i,j)
            rhsw3(i,j)=frc_wdp3(i,j)
            rhswr(i,j)=frc_wcnp(i,j)
         END DO
         END DO
      END IF

      CALL depst3l2d_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,             &
                         tsfcn,tdeepn,wgn,w2n,w3n,wrn,                &
                         snowdpth,qvsfc,windsp,psfc,rhoa,precip,      &
                         tair,qvair,cdha,cdqa,radsw,                  &
                         rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,    &
                         evaprr,qvsat,qvsata,f34,                     &
                         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp2,        &
                         frc_wdp3,frc_wcnp,wrmax,relief,              &
                         SatFlow,RunOffei,RunOffCnp,wtdepthmap,       &
                         sdepthmap)
                         
!  Integration for Ts and T2 at one time step.

      DO j = 1, ny-1
      DO i = 1, nx-1
         IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
            rhsts(i,j) = 0.0
            rhst2(i,j) = 0.0
         ELSE
            rhsts(i,j) = 0.5 * (frc_tsfc(i,j) +rhsts(i,j))
            rhst2(i,j) = 0.5 * (frc_tdeep(i,j)+rhst2(i,j)) 
         END IF
         tsfc(i,j)  = tsfc(i,j) + dtsfc * rhsts(i,j)
         tdeep(i,j) = tdeep(i,j)+ dtsfc * rhst2(i,j)
      END DO
      END DO
      
 
      IF ( moist /= 0 ) THEN
         DO j = 1, ny-1
         DO i = 1, nx-1

            IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.         &
                 snowdpth(i,j) >= snowdepth_crit ) THEN
               rhswg(i,j) = 0.0
               rhsw2(i,j) = 0.0
               rhsw3(i,j) = 0.0
               rhswr(i,j) = 0.0
            ELSE
               rhswg(i,j) = 0.5 * (frc_wsfc(i,j)+rhswg(i,j))
               rhsw2(i,j) = 0.5 * (frc_wdp2(i,j)+rhsw2(i,j))
               rhsw3(i,j) = 0.5 * (frc_wdp3(i,j)+rhsw3(i,j))
               rhswr(i,j) = 0.5 * (frc_wcnp(i,j)+rhswr(i,j))
            END IF
            
            wetsfc(i,j) = wetsfc(i,j) + dtsfc * rhswg(i,j)
            wetsfc(i,j) = MAX( wetsfc(i,j), 0.0 )
            wetsfc(i,j) = MIN( wetsfc(i,j), wsat(soiltyp(i,j)) )

            wetdp2(i,j) = wetdp2(i,j) + dtsfc * rhsw2(i,j)
            wetdp2(i,j) = MAX( wetdp2(i,j), 0.0 )
            wetdp2(i,j) = MIN( wetdp2(i,j), wsat(soiltyp(i,j)) )
          
            wetdp3(i,j) = wetdp3(i,j) + dtsfc * rhsw3(i,j)
            wetdp3(i,j) = MAX( wetdp3(i,j), 0.0 )
            wetdp3(i,j) = MIN( wetdp3(i,j), wsat(soiltyp(i,j)) )

            wetcanp(i,j) = wetcanp(i,j) + dtsfc * rhswr(i,j)
            wetcanp(i,j) = MAX( wetcanp(i,j), 0.0 )
            wetcanp(i,j) = MIN( wetcanp(i,j), wrmax(i,j) )

	 END DO
         END DO

      END IF

   END DO  ! TIME INTEGRATION


   DO j = 1, ny-1                               !SOIL
   DO i = 1, nx-1
      tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

      IF (snowdpth(i,j) >= snowdepth_crit) THEN
         ct(i,j) = cg_snow
         gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)             &
                  * snowflxfac/(tau*ct(i,j))
                                        ! Snow cover
      ELSE
         cg = cgsat(soiltyp(i,j))                                     &
              * ( wsat(soiltyp(i,j))/tema )                           &
              **( bslope(soiltyp(i,j))/log100 )

        ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv                  &
                           + veg(i,j) * cg )       ! NP, Eq. 8

!         ct(i,j) = cg                                  !(PX1995) 

         gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))         
                                        ! Ground diffusive heat flux
      END IF

      shflx(i,j) =rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)            &
         * ( tsfc(i,j) - tair(i,j) *(psfc(i,j)/1.0E5)**rddcp )         !RDDRDD

   END DO
   END DO

   CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1,psfc,tsfc,qvsat)             !HYDROLOGY

   CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,          &
        soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp2,wetcanp,            &
        snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)

   DO j=1, ny-1
   DO i=1, nx-1
      qvsfc(i,j) = lhflx(i,j) + qvair(i,j)
      evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
      evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
      evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
      lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)           !MOISTURE
      lhflx  (i,j) = lhflx(i,j) * lathv                                   !QE 
   END DO
   END DO

   EBCal = rnflx(sx,sy) - lhflx(sx,sy) - shflx(sx,sy) - gflx(sx,sy)   

   RETURN
END SUBROUTINE Depst3l2d

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE depst3l2d_frc              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE depst3l2d_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,         &
                         tsfc,tdeep,wetsfc,wetdp2,wetdp3,wetcanp, &
                         snowdpth,qvsfc,windsp,psfc,rhoa,precip,  &
                         tair,qvair,cdha,cdqa,radsw,              &
                         rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,&
                         evaprr,qvsat,qvsata,f34,                 &
                         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp2,    &
                         frc_wdp3,frc_wcnp,wrmax,relief,          &  
                         SatFlow,RunOffei,RunOffCnp,wtdepthmap,   &
                         sdepthmap)
!------------------------------------------------------------------
!                                                                  
!  AUTHOR: Diandong Ren (dd_ren@rossby.metr.ou.edu)                
!          Ming Xue     (mxue@ou.edu)
!  12/08/2001                                                      
!                                                                  
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!------------------------------------------------------------------
!                                                                  
!  PURPOSE:                                                        
!  To calculate the right hand side forcing terms in the 
!  soil-vegetation model.
!------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER :: nx,ny,nz
!   REAL :: eta_fac,wmax,c1max,sigma_sqr

   REAL :: windsp(nx,ny)         ! Wind speed
   REAL :: psfc  (nx,ny)         ! Surface pressure
   REAL :: rhoa  (nx,ny)         ! Air density near the surface
   REAL :: precip(nx,ny)         ! Precipitation rate at the surface

   INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
   INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

   REAL :: lai    (nx,ny)        ! Leaf Area Index
   REAL :: veg    (nx,ny)        ! Vegetation fraction

   REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
   REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
   REAL :: wetsfc (nx,ny)        ! Surface soil moisture
   REAL :: wetdp2 (nx,ny)       ! Deep soil moisture
   REAL :: wetdp3 (nx,ny)
   REAL :: wetcanp(nx,ny)        ! Canopy water amount
   REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

   REAL :: snowdpth(nx,ny)       ! Snow depth (m)

   REAL :: tair   (nx,ny)        ! Surface air temperature (K)
   REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
   REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
   REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

   REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
   REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
   REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
   REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
   REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
   REAL :: ct     (nx,ny)        ! Soil thermal coefficient

   REAL :: evaprg (nx,ny)        ! Evaporation
   REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
   REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
   REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
   REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
   REAL :: qvsat  (nx,ny)        !

   REAL :: frc_tsfc(nx,ny)       ! Right hand side forcing for tsfc eq.
   REAL :: frc_tdeep(nx,ny)      ! Right hand side forcing for tsoil eq.
   REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq. 
   REAL :: frc_wdp2(nx,ny)       ! Right hand side forcing for wetsp eq. 
   REAL :: frc_wdp3(nx,ny)       ! Right hand side forcing for wetsp eq. 
   REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.
   REAL :: wrmax(nx, ny)         ! Maximum value for canopy moisture, wetcanp
  
   REAL :: relief                ! Difference between seasonal average skin and deep soil 
                                 ! temperature  (t_skin - t_deeplayer).
   real :: SatFlow(nx,ny)
   real :: RunOffei(nx,ny)
   real :: RunOffCnp(nx,ny)
   real :: wtdepthmap(nx,ny)
   real :: sdepthmap(nx,ny)                        

   REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
   REAL :: c2wg        ! Coefficient in the surface moisture Eq. of W2
   real :: w3s         ! Coefficient in the surface moisture Eq. of W3
   real :: c3          ! Coefficient in the surface moisture Eq. of W3
   real :: D12         !
   real :: D23         !
   real :: K2          !
   real :: K3          !

   REAL :: pi          ! Pi
   PARAMETER (pi = 3.141592654)

!   REAL :: tau         ! Seconds of a day = 24. * 3600.
!   PARAMETER (tau = 86400.)

   REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                       ! dtsfc2 = dtsfc/2.
   REAL :: log100      ! Constant: alog(100)
                       ! dependent distance from the earth to the sun
   REAL :: cg          ! Soil thermal coefficient for bare ground
   REAL :: wgeq        ! Surface moisture when gravity balances the capillarity
   REAL :: wr2max      ! Tendency to reach the maximum wrmax
   REAL :: pnet        ! Residual of precip. and evap.
   real :: pgrd        ! Precip. that reaches the ground
   real :: pref        ! precipitação de referência (Habets et.al. (1999) jh217-75
   REAL :: vegp        ! Precip. intercepted by vegetation
!-----------------------------------------------------------------------
!	real               ::   NRootLayers
!       real,allocatable   ::   LPorosity(:)
!	real,allocatable   ::   LFCap(:)
!	real,allocatable   ::   RootDepth(:)
!	real,allocatable   ::   Moisture(:)
!       real(4)            ::   CalcWTableDepth     
!
!-----------------------------------------------------------------------
!  Variáveis auxiliares - Auxiliary variables
!-----------------------------------------------------------------------
   INTEGER :: i,j,k

   REAL :: tema,temb
   REAL :: eta, c1max, wmax, sig2 
   integer :: istatus
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
   INCLUDE 'globcst.inc'
   INCLUDE 'phycst.inc'
   INCLUDE 'soilcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   write(6,'(a)') ' Entrando em depst3l2d_frc!'
   log100    = ALOG(100.0)
   dtsfc2    = dtsfc/2.0

!-----------------------------------------------------------------------  
   DO j = 1, ny-1
   DO i = 1, nx-1
      tema = MAX( wetdp2(i,j), wwlt(soiltyp(i,j)) )
  
      IF (snowdpth(i,j) >= snowdepth_crit) THEN !snow cover
        ct(i,j) = cg_snow
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                     &
                          *snowflxfac/(tau*ct(i,j))
      ELSE
        cg = cgsat(soiltyp(i,j))                                      &
            * ( wsat(soiltyp(i,j))/tema )                             &
            **( bslope(soiltyp(i,j))/log100 )
!        write(6,'(a,4f18.8,i4)') 'cg=', cg,cgsat(soiltyp(i,j)),wsat(soiltyp(i,j)),bslope(soiltyp(i,j)),soiltyp(i,j)
        ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv + veg(i,j) * cg )     
!        if (i==sx .and. j==sy) then
!           write(6,'(a,5f18.8)') 'ct =', ct(i,j),cg,cgv,veg(i,j),bslope(soiltyp(i,j))
!        endif
	
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
	
!        if (i==sx .and. j==sy) then		
!           write(6,'(a)') '    '
!	   write(6,'(a,2i4,5f18.8)') 'gflx =', i,j,gflx(i,j),tsfc(i,j),tdeep(i,j),relief,ct(i,j)
!	endif
      END IF
  
      shflx(i,j) = rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)           &
                  * ( tsfc(i,j) - tair(i,j) )  
!      if (i==sx .and. j==sy) then
!         write(6,'(a,2i4,6f18.8)')' shflx=',i,j,shflx(i,j),           &
!	       cp,cdha(i,j),windsp(i,j),tsfc(i,j),tair(i,j)
!      endif		  
   END DO
   END DO

   IF ( moist == 0 ) THEN

      DO j = 1, ny-1
      DO i = 1, nx-1
         evaprg(i,j) = 0.0    
         evaprtr(i,j) = 0.0   
         evaprr(i,j) = 0.0  
         lhflx(i,j) = 0.0   
      END DO
      END DO

   ELSE

      CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)

      CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,       &
                   soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp2,wetcanp,  &
                   snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)
   END IF

   DO j = 1, ny-1
   DO i = 1, nx-1
      evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
      evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
      evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
      lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)
      
      lhflx(i,j) = lhflx(i,j) * lathv       ! Latent heat flux
   END DO
   END DO

   DO j = 1, ny-1
   DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        frc_tsfc (i,j) = 0.
        frc_tdeep(i,j) = 0.
      ELSE

      
         frc_tsfc(i,j) = ct(i,j)*(rnflx(i,j)-shflx(i,j)-lhflx(i,j)-gflx(i,j))

         
         IF ( snowdpth(i,j) >= snowdepth_crit ) THEN
            frc_tdeep(i,j)= 1.03* (tsfc(i,j)-tdeep(i,j)-relief)*tauinv*snowflxfac ! Snow cover
         ELSE
!            frc_tdeep(i,j)= 1.03*(tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
	    frc_tdeep(i,j)= (tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
        END IF
      END IF
   END DO
   END DO
   
!---------------------------------------------------------------------  
! Modelo de solo - Soil model
!---------------------------------------------------------------------   
   IF ( moist /= 0 ) THEN   

      DO j = 1, ny-1                                       !HYDROLOGY
      DO i = 1, nx-1
         wrmax(i,j) = (0.2 * 1.e-3 * veg(i,j) * lai(i,j))/dr     ! meter
         IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.        &
              snowdpth(i,j) >= snowdepth_crit ) THEN
            frc_wsfc(i,j) = 0.
            frc_wdp2(i,j) = 0.
            frc_wdp3(i,j) = 0.
            frc_wcnp(i,j) = 0.
         ELSE
!-----------------------------------------------------------------------
!  Parcelas de Precipitacao - Precipitations parcels
!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!  Precipitação sobre o solo

            if (wetcanp(i,j) > wrmax(i,j)) then 
               pgrd = (1-veg(i,j)) * precip(i,j) + (wetcanp(i,j) -    &
                  wrmax(i,j))*(rhow*dr/dtsfc)
		  
            else
               pgrd = (1-veg(i,j))*precip(i,j)
            end if
            
!---------------------------------------------------------------------    
! Precipitação sobre a vegetação (vegp), precipitação líquida (pnet) e
! Geração de escoamento superficial RunOffCnp
     
            wr2max = ( wrmax(i,j) - wetcanp(i,j) ) / dtsfc2    
            vegp = veg(i,j) * precip(i,j)

            pnet = vegp - evaprr(i,j)

            tema = pnet - wr2max * (rhow)   !  Eu acrescentei *dr ?(rhow*dr)??  Rio 24/01/06

            RunOffCnp(i,j) = MAX( tema, 0.0 )   ! formulação do SvatArps(soilebm3d)            
            
!---------------------------------------------------------------------    
! Geração de RunOff de acordo com Habets.ea:1999
!---------------------------------------------------------------------    
            tema = (1 - ((wetdp2(i,j)-wwlt(soiltyp(i,j))) /           &
	           (wsat(soiltyp(i,j))-wwlt(soiltyp(i,j)))))**        &
		   (1/(Bvic+1))
		   
            pref = (1+Bvic)*(wsat(soiltyp(i,j))-wwlt(soiltyp(i,j))) * &
	            tema
	     
            temb = (d2*(wsat(soiltyp(i,j))-wetdp2(i,j)))              &
	            + ((d2*(wsat(soiltyp(i,j))-wwlt(soiltyp(i,j)))) * &
		    (pref-(pgrd*(dtsfc/rhow*d2)/((wsat(soiltyp(i,j))- &
		    wwlt(soiltyp(i,j)))*(Bvic+1))))**(1+Bvic))


            if (pgrd*(dtsfc/rhow*d2) > pref) then
	       
	       RunOffei(i,j) = pgrd - (d2*(wsat(soiltyp(i,j))-        &
                                wetdp2(i,j))) * (rhow/dtsfc)
               
            else if (pgrd*(dtsfc/rhow*d2) .le. pref) then
	    
               RunOffei(i,j) = pgrd - temb * (rhow*d2/dtsfc)
				
            endif
	    
	    if (RunOffei(i,j) < 0.0 )  then
	       RunOffei(i,j) = 0.0
	    endif

!-----------------------------------------------------------------------
! Coeficientes para umidade - Soil moisture coefficients
!-----------------------------------------------------------------------
            tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

            c1wg = 0.4* c1sat(soiltyp(i,j))                           &
                   * ( wsat(soiltyp(i,j)) / tema )                    &
                   **( bslope(soiltyp(i,j)) / 2.0 + 1.0)

!-----------------------------------------------------------------------
! Replacement Cl to improve dry soils (NM1996 - A.3)   (JAB)
!-----------------------------------------------------------------------
            IF (wetsfc(i,j) < wwlt(soiltyp(i,j)) ) THEN

               eta = (-0.01815 * tsfc(i,j) + 6.41 ) * wwlt(soiltyp(i,j)) + &
                     ( 0.0065 * tsfc(i,j) - 1.4 )
      
               c1max = (1.19 * wwlt(soiltyp(i,j)) - 5.09)*0.01*tsfc(i,j)   &
                      +( 1.464 * wwlt(soiltyp(i,j)) + 17.86 )
      
               wmax = eta * wwlt(soiltyp(i,j))
      
               sig2 = - ( (2 * alog ( 0.01 / c1max )) / ( wmax*wmax) )
      
               c1wg = c1max * EXP ( -0.5* ( (wetsfc(i,j)-wmax)*(      &
                      wetsfc(i,j)-wmax)*sig2 ) )
            ENDIF

!----------------------------------------------------------------------------
            c2wg = c2ref(soiltyp(i,j)) * wetdp2(i,j)                  &
                   / ( wsat(soiltyp(i,j)) - wetdp2(i,j) + wetsml )
                   
            wgeq = wetdp2(i,j) - wsat(soiltyp(i,j))                   &
                   * awgeq(soiltyp(i,j))                              &
                   * ( wetdp2(i,j) / wsat(soiltyp(i,j)) )             &
                   ** pwgeq(soiltyp(i,j))                             &
                   * ( 1.0 - ( wetdp2(i,j) / wsat(soiltyp(i,j)) )     &
                   ** ( 8 * pwgeq(soiltyp(i,j)) ) )
                   
!----------------------------------------------------------------------------
                  
            w3s = wfc(soiltyp(i,j)) + (wsat(soiltyp(i,j)) -           &
	               wfc(soiltyp(i,j))) / exp(1.0)  
                    
            c3 = (tau * (2*bslope(soiltyp(i,j)) + 2) *                &
                  ks(soiltyp(i,j))) / (d3 * ((W3s /                   &
                  wsat(soiltyp(i,j))) ** (-2 * bslope(soiltyp(i,j))   &
                  -  2) - 1))
!-----------------------------------------------------------------------
! Difusão vertical subsuperficial - Subsurface vertival diffusion
!-----------------------------------------------------------------------                            
            D12 = (c2wg/tau) * (wetsfc(i,j) - wgeq)
          
            D23 = (c4(soiltyp(i,j))/tau) * (wetdp2(i,j) - wetdp3(i,j))

!-----------------------------------------------------------------------
! Drenagem gravitacional - Gravitational dreinage
!-----------------------------------------------------------------------
            K2 = (c3/tau) * (d3/d2) * MAX(0.0,(wetdp2(i,j) -          &
                  wfc(soiltyp(i,j))))
  
            K3 = (c3/tau) * (d3/(d3-d2)) * MAX(0.0,(wetdp3(i,j) -     &
                  wfc(soiltyp(i,j))))

!-----------------------------------------------------------------------  
!  Equações Básicas do Modelo de água no solo
!-----------------------------------------------------------------------
            frc_wsfc(i,j) = c1wg * ( pgrd - RunOffei(i,j) -          &
                             evaprg(i,j) ) / (rhow * d1) - D12


            frc_wdp2(i,j) = ((pgrd - RunOffei(i,j) - evaprg(i,j) -    &
	                    evaprtr(i,j) ) / (rhow * d2)) - K2 - D23 


            frc_wdp3(i,j) =  (d2/(d3-d2)) * (K2 + D23) - K3

         
            frc_wcnp(i,j) = ( vegp - evaprr(i,j) - RunOffCnp(i,j) ) / (rhow*dr)
            
         END IF

      END DO
      END DO           
      

   END IF

   RETURN
END SUBROUTINE depst3l2d_frc

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SOILEBM                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE soilebm(nx,ny,nz,soiltyp,vegtyp,lai,veg,                     &
           tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,                    &
           qvsfc,windsp,psfc,rhoa,precip,                               &
           tair,qvair,cdha,cdqa,radsw,rnflx,shflx,lhflx,gflx,ct,        &
           evaprg,evaprtr,evaprr,qvsat,qvsata,f34)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Predict the soil surface temperature and moisture contents by solving
!  the surface energy and moisture budget equtions:
!
!  1. the ground surface temperature, Ts -- tsfc
!  2. the deep ground temperature,    T2 -- tdeep
!  3. the surface soil moisture,      wg -- wetsfc
!  4. the deep soil moisture,         w2 -- wetdp
!  5. the canopy moisture,            wr -- wetcanp
!
!-----------------------------------------------------------------------
!
!  The equations are listed as follows.
!
!
!       d Ts                        2PI
!      ------ = Ct (Rn - H - LE) - ----- (Ts - T2)
!       d t                         Tau
!
!
!       d T2      1
!      ------ = ----- (Ts - T2)
!       d t      Tau
!
!
!       d Wg      C1                 C2
!      ------ = ------ (Pg - Eg) - ----- (Wg - Wgeq)
!       d t     ROw d1              Tau
!
!
!       d W2      1
!      ------ = ------ (Pg - Eg - Etr)
!       d t     ROw d2
!
!
!       d Wr
!      ------ = Veg P - Er
!       d t
!
!
!  where
!
!      Tau    -- 1 day in seconds = 86400 seconds
!      PI     -- number of PI = 3.141592654
!      ROw    -- Density of liquid water
!      d1     -- Top layer depth of soil column, 0.01 m
!      d2     -- Deep layer depth of soil column, 1 m
!      Veg    -- Vegetation fraction
!      Ct     -- Thermal capacity
!      Rn     -- Radiation flux, rnflx
!      H      -- Sensible heat flux, shflx
!      LE     -- Latent heat flux, lhflx = latent*(Eg + Ev)
!      Eg     -- Evaporation from ground
!      Ev     -- Evapotranspiration from vegetation, Ev = Etr + Er
!      Etr    -- Transpiration of the remaining part of the leaves
!      Er     -- Evaporation directly from the foliage covered by
!                intercepted water
!      P      -- Precipitation rates
!      Pg     -- Precipitation reaching the ground,
!                Pg = (1 - Veg) P
!      Wgeq   -- Surface volumetric moisture
!      C1, C2 -- Coefficients
!
!  For detailed information about the surface energy budget model,
!  see the articles in the reference list.
!
!  The second-order Rouge-Kutta time integration scheme is used,
!  which is described below.
!
!  Assume a equation in the form of
!
!      d X
!     ----- = F(X, t)
!      d t
!
!  In the forward scheme, we have
!
!      X(1) = X(0) + dt * F[X(0), t0]
!
!  We split one time step into two halves, dt2 = dt/2, and use the
!  forward scheme to calculate the first half step X(1/2).
!
!      X(1/2)  = X(0) + dt2 * F[X(0), t0]
!
!  Then we can calculate the Right Hand Side (RHS) of the equation at
!  the half step, F[X(1/2), t(1/2)]. Finally, we calculate the one
!  step prediction, X(t1), by use of the average of F[X(0), t0] and
!  F[X(1/2), t(1/2)].
!
!      X(1) = X(0) + dt * 0.5 * { F[X(0),t0] + F[X(1/2), t(1/2)] }
!  REFERENCES:
!
!  Jacquemin, B. and J. Noilhan, 1990: Sensitivity Study and
!       Validation of a Land Surface Parameterization Using the
!       HAPEX-MOBILHP Data Set, Boundary-Layer Meteorology, 52,
!       93-134, (JN).
!
!  Noilhan, J. and S. Planton, 1989: A Simple Parameterization of
!       Land Surface Processes for Meteorological Model, Mon. Wea.
!       Rev., 117, 536-549, (NP).
!
!  Pleim, J. E. and A. Xiu, 1993: Development and Testing of a
!       Land-Surface and PBL Model with Explicit Soil Moisture
!       Parameterization, Preprints, Conf. Hydroclimat., AMS, 45-51,
!       (PX).
!
!  Bougeault, P., et al., 1991: An Experiment with an Advanced
!       Surface Parameterization in a Mesobeta-Scale Model. Part I:
!       Implementation, Monthly Weather Review, 119, 2358-2373.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  11/16/93
!
!  MODIFICATION HISTORY:
!
!  10/30/94 (Y. Liu)
!  Fixed a bug reported by Richard Carpenter.
!
!  11/01/1994 (Y. Liu)
!  Subroutine name of soil model were changed from SFCEBM to SOILEBM
!  and the arguments passed were changed to 2-D arrays instead of 3-D
!  temporary arrays.
!
!  12/12/1994 (Y. Liu)
!  Fixed a bug for the final calcultaion of qvsfc.
!
!  12/14/1994 (Y. Liu)
!  Fixed a bug in phycst.inc for the water density value which was
!  previously mistakingly set to 1 kg/m**3. The correct value should
!  be 1000 kg/m**3. This bug largely influenced the integration of Wg
!  and W2.
!
!  12/23/1994 (Y. Liu)
!  Added the runoff calculation for Wr.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the soil model and
!  at the same time delete the table data array veg(14).
!
!  03/27/1995 (Yuhe Liu)
!  Changed the solor radiation used in the calculation of surface
!  resistence factor F1 from the one at the top of atmosphere to the
!  one at the surface.
!
!  Changed the formula of calculating the surface resistence factor
!  F3 to F3=1, instead of varying with qvsat(Tair) and qvair.
!
!  12/8/1998 (Donghai Wang and Vince Wong)
!  Added a new 2-D permanent array, snowcvr(nx,ny), for snow cover.
!  We just used a simple scheme to consider the snow cover process.
!
!  2000/01/10 (Gene Bassett)
!  Snow cover (0 or 1) changed to snow depth (snowdpth).  For simplicity
!  a fractional value for snow cover is not used (simply say the grid
!  point is completely covered with snow if snowdpth > snowdepth_crit,
!  otherwise no snow).
!
!  2000/02/04 (Gene Bassett, Yang Kun)
!  Fixed an error in tsoil integration (rhst2) causing tsoil to change
!  by only 50% of what it should.
!
!  2001/12/07 (Diandong Ren, Ming Xue) 
!  Re-structured the code, moved the calculations of the right hand
!  side terms of the soil-vegetation model into subroutine soilebm_frc.
!
!  Soil seasonal temperature trend to be added.
!
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!  2002/06/9 (Ming Xue)
!  Revoked some modifications to the soil moisture related caculations
!  that Diandong Ren put in since IHOP_2 - the mods need more testing.
!
!  2002/12/13 (Jerry Brotzge) 
!  Updated code to match recommendations by Pleim and Xiu (1995) and 
!  Xiu and Pleim (2001 - JAM).  
!
!-----------------------------------------------------------------------
! 
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (sfc/top)
!
!    soiltyp  Soil type at the horizontal grid points
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    windsp   Wind speed just above the surface (m/s)
!    psfc     Surface pressure (Pascal)
!    rhoa     Near sfc air density
!    precip   Precipitation flux reaching the surface (m/s)
!
!    cdha     Array for cdh, surface drag coefficient for heat
!    cdqa     Array for cdq, surface drag coefficient for moisture
!
!    pres     3-dimensional pressure
!    temp     3-dimensional temperature
!    qv       3-dimensional specific humidity
!
!  INPUT/OUTPUT: 
!
!    tsfc     Temperature at ground surface (K)
!    tdeep    Deep soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Canopy water amount
!
!  OUTPUT:
!
!    qvsfc    Effective S. H. at sfc.
!
!  Local automatic work arrays for storing the right hand forcing terms 
!  of the soil model equations and for storing intermediate values of the
!  soil state variables. 
!
!    frc_tsfc   Temporary array
!    frc_tdeep  Temporary array
!    frc_wsfc   Temporary array
!    frc_wdp    Temporary array
!    frc_wcnp   Temporary array
!
!    tsfcn      Temporary array
!    tdeepn     Temporary array
!    wgn        Temporary array
!    w2n        Temporary array
!    wrn        Temporary array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny)        ! Surface soil moisture
  REAL :: wetdp  (nx,ny)        ! Deep soil moisture
  REAL :: wetcanp(nx,ny)        ! Canopy water amount
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !


  REAL :: frc_tsfc(nx,ny)       ! Right hand side forcing for tsfc eq.
  REAL :: frc_tdeep(nx,ny)      ! Right hand side forcing for tdeep eq.
  REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq.
  REAL :: frc_wdp(nx,ny)        ! Right hand side forcing for wetsp eq.
  REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.

  REAL :: tsfcn(nx,ny)          ! Temporary array, tsn
  REAL :: tdeepn(nx,ny)         ! Temporary array, t2n
  REAL :: wgn(nx,ny)            ! Temporary array, wetsfcNEW
  REAL :: w2n(nx,ny)            ! Temporary array, w2n
  REAL :: wrn(nx,ny)            ! Temporary array, wrn
  REAL :: relief          ! Difference between seasonal average skin and deep soil 
!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
!  Parameters and variables are defined in globcst.inc:
!
!    dtsfc      Surface model time step
!    nsfcst     # of surface model time steps
!
!    moist      Moist flag
!
!    year       Reference year
!    month      Reference month
!    day        Reference day
!    jday       Reference Julian day
!    hour       Hour of reference time
!    minute     Minute of reference time
!    second     Second of reference time
!
!    latitud    Latitude at the domain center
!    longitud   Longitude at the domain center
!
!    curtim     Current model time
!    dtbig      Length of big time step
!
!    bslope     Slope of the retention curve
!    cgsat      Soil thermal coefficient for bare ground at saturation
!    cgv        Soil thermal coef. for totally shielded ground by veg.
!    pwgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    awgeq      Coefficient of Wgeq formula. NP, Tab. 2
!    c1sat      Value of C1 at saturation. NP, Tab. 2
!    c2ref      Value of C2 for W2 = .5 * Wsat. NP, Tab. 2
!    wsat       Saturated volumetric moisture content. JN, Tab. 1
!    wfc        Field capacity moisture. JN, Tab. 1
!    wwlt       Wilting volumetric moisture content. JN, Tab. 1
!
!  Parameters and variables are defined in phycst.inc:
!
!    solarc     Solar constant (W/m**2)
!    emissg     Emissivity of the ground
!    emissa     Emissivity of the atmosphere
!    sbcst      Stefen-Boltzmann constant
!
!    rhow       Liquid water reference density (kg/m**3)
!    rd         Gas constant for dry air (kg/(m s**2))
!    cp         Gas heat capacity at constant pressure
!    cv         Gas heat capacity at constant volume
!
!-----------------------------------------------------------------------
!
   INCLUDE 'globcst.inc'
   INCLUDE 'phycst.inc'
   INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
   REAL :: pi                         ! Pi
   PARAMETER (pi = 3.141592654)

!  REAL :: tau                        ! Seconds of a day = 24. * 3600.
!  PARAMETER (tau = 86400.)

  REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.

  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun

  REAL :: cg          ! Soil thermal coefficient for bare ground

  REAL :: rhsts(nx,ny)       ! Right hand side of Eq. for Ts at current time
  REAL :: rhst2(nx,ny)       ! Right hand side of Eq. for T2 at current time
  REAL :: rhswg(nx,ny)       ! Right hand side of Eq. for Wg at current time
  REAL :: rhsw2(nx,ny)       ! Right hand side of Eq. for W2 at current time
  REAL :: rhswr(nx,ny)       ! Right hand side of Eq. for Wr at current time
  REAL :: wrmax(nx,ny)       ! Maximum value for canopy moisture, wetcanp
  REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
  REAL :: c2wg        ! Coefficient in the surface moisture Eq. of Wg

  REAL :: wgeq        ! Surface moisture when gravity balances the capillarity

  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation

  INTEGER :: i, j, it

  REAL :: tema
  real :: EBCal

  LOGICAL :: firstcall        ! First call flag of this subroutine

!  SAVE firstcall, log100, dtsfc2, tauinv
  SAVE firstcall, log100, dtsfc2
  DATA firstcall/.true./

  INTEGER :: jday_min         ! offset value from Jan 01. 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (firstcall) THEN
    log100    = ALOG(100.0)
    dtsfc2    = dtsfc/2.0
    firstcall = .false.
  END IF

  IF ( moist /= 0 ) THEN

!  Calculate saturated specific humidity near the surface, qvsata.

    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tair,qvsata)


    DO j = 1, ny-1
    DO i = 1, nx-1
!      f34(i,j) = MAX( 0., 1.0 - 0.0016 * (298.0-tair(i,j))**2 )

       IF (tair(i,j) >= 302.15) THEN                    !(XP2001) 
         f34(i,j) = 1.0 / (1.0 + EXP( 1.18 * (tair(i,j)-314.00)))
       ELSE
         f34(i,j) = 1.0 / (1.0 + EXP( -0.41* (tair(i,j)-282.05)))
       ENDIF 

    END DO
    END DO
  END IF

  jday_min = 61       ! may change with latitude
  relief=tsoil_offset_amplitude*sin((jday-jday_min)/365.0*2.0*PI)

!
!-----------------------------------------------------------------------
!
!  Start time integration loop
!
!-----------------------------------------------------------------------
!
  DO it = 1, nsfcst               

    CALL soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,                   &
         tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,                      &
         qvsfc,windsp,psfc,rhoa,precip,tair,qvair,cdha,cdqa,            &
         radsw,rnflx,shflx,lhflx,gflx,ct,evaprg,evaprtr,                &
         evaprr,qvsat,qvsata,f34,                                       &
         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

    DO j = 1, ny-1
    DO i = 1, nx-1
      tsfcn (i,j) = tsfc(i,j) + dtsfc2 * frc_tsfc(i,j)
      tdeepn(i,j) = tdeep(i,j)+ dtsfc2 * frc_tdeep(i,j)
    END DO
    END DO

    IF ( moist /= 0 ) THEN    

      DO j = 1, ny-1
      DO i = 1, nx-1
        wgn(i,j) = wetsfc(i,j) + dtsfc2 * frc_wsfc(i,j)
        wgn(i,j) = MAX(wgn(i,j), 0.0 )
        wgn(i,j) = MIN(wgn(i,j), wsat(soiltyp(i,j)) )

        w2n(i,j) = wetdp(i,j) + dtsfc2 * frc_wdp(i,j)
        w2n(i,j) = MAX( w2n(i,j), 0.0 )
        w2n(i,j) = MIN(w2n(i,j), wsat(soiltyp(i,j)) )

        wrn(i,j) = wetcanp(i,j) + dtsfc2 * frc_wcnp(i,j)
        wrn(i,j) = MAX(wrn(i,j), 0.0 )
        wrn(i,j) = MIN(wrn(i,j), wrmax(i,j) )
      END DO
      END DO

    ELSE

      DO j = 1, ny-1
        DO i = 1,nx-1
          w2n(i,j) = wetdp(i,j)
        END DO
      END DO

    END IF

    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        rhsts(i,j) = 0.0
        rhst2(i,j) = 0.0
      ELSE
        rhsts(i,j)=frc_tsfc(i,j)
        rhst2(i,j)=frc_tdeep(i,j)
      END IF
    END DO
    END DO

    IF ( moist /= 0 ) THEN    
      DO j = 1, ny-1
      DO i = 1, nx-1
        rhswg(i,j)=frc_wsfc(i,j)
        rhsw2(i,j)=frc_wdp(i,j)
        rhswr(i,j)=frc_wcnp(i,j)
      END DO
      END DO
    END IF

    CALL soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,                   &
         tsfcn,tdeepn,wgn,w2n,wrn,snowdpth,                             &
         qvsfc,windsp,psfc,rhoa,precip,tair,qvair,cdha,cdqa,            &
         radsw,rnflx,shflx,lhflx,gflx,ct,                               &
         evaprg,evaprtr,evaprr,qvsat,qvsata,f34,                        &
         frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

!  Integration for Ts and T2 at one time step.

    DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        rhsts(i,j) = 0.0
        rhst2(i,j) = 0.0
      ELSE
        rhsts(i,j) = 0.5 * (frc_tsfc(i,j) +rhsts(i,j))
        rhst2(i,j) = 0.5 * (frc_tdeep(i,j)+rhst2(i,j)) 
      END IF
        tsfc(i,j)  = tsfc(i,j) + dtsfc * rhsts(i,j)
        tdeep(i,j) = tdeep(i,j)+ dtsfc * rhst2(i,j)
    END DO
    END DO

    IF ( moist /= 0 ) THEN
      DO j = 1, ny-1
      DO i = 1, nx-1
          IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.         &
            snowdpth(i,j) >= snowdepth_crit ) THEN
            rhswg(i,j) = 0.0
            rhsw2(i,j) = 0.0
            rhswr(i,j) = 0.0
          ELSE
            rhswg(i,j) = 0.5 * (frc_wsfc(i,j)+rhswg(i,j))
            rhsw2(i,j) = 0.5 * (frc_wdp (i,j)+rhsw2(i,j))
            rhswr(i,j) = 0.5 * (frc_wcnp(i,j)+rhswr(i,j))
          END IF
          wetsfc(i,j) = wetsfc(i,j) + dtsfc * rhswg(i,j)
          wetsfc(i,j) = MAX( wetsfc(i,j), 0.0 )
          wetsfc(i,j) = MIN( wetsfc(i,j), wsat(soiltyp(i,j)) )

          wetdp(i,j) = wetdp(i,j) + dtsfc * rhsw2(i,j)
          wetdp(i,j) = MAX( wetdp(i,j), 0.0 )
          wetdp(i,j) = MIN( wetdp(i,j), wsat(soiltyp(i,j)) )

          wetcanp(i,j) = wetcanp(i,j) + dtsfc * rhswr(i,j)
          wetcanp(i,j) = MAX( wetcanp(i,j), 0.0 )
          wetcanp(i,j) = MIN( wetcanp(i,j), wrmax(i,j) )
       END DO
      END DO
    END IF

  END DO  ! TIME INTEGRATION

  DO j = 1, ny-1                               !SOIL
  DO i = 1, nx-1
    tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

    IF (snowdpth(i,j) >= snowdepth_crit) THEN
      ct(i,j) = cg_snow
      gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                 &
                  *snowflxfac/(tau*ct(i,j))
                                        ! Snow cover
    ELSE
      cg = cgsat(soiltyp(i,j))                                        &
          * ( wsat(soiltyp(i,j))/tema )                               &
          **( bslope(soiltyp(i,j))/log100 )

      ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv                     &
                           + veg(i,j) * cg )       ! NP, Eq. 8

!      ct(i,j) = cg                                  !(PX1995) 

      gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
                                        ! Ground diffusive heat flux
    END IF

      shflx(i,j) =rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)             &
         * ( tsfc(i,j) - tair(i,j) *(psfc(i,j)/1.0E5)**rddcp )         !RDDRDD


  END DO
  END DO

  CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)             !HYDROLOGY

  CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,            &
       soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp,wetcanp,               &
       snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)

  DO j=1, ny-1
  DO i=1, nx-1
    qvsfc(i,j) = lhflx(i,j) + qvair(i,j)
    evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
    evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
    evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
    lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)           !MOISTURE
    lhflx  (i,j) = lhflx(i,j) * lathv                                   !QE 
  END DO
  END DO

   EBCal = rnflx(sx,sy) - lhflx(sx,sy) - shflx(sx,sy) - gflx(sx,sy)   

  RETURN
END SUBROUTINE soilebm


!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SOILEBM_FRC               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE soilebm_frc(nx,ny,nz,soiltyp,vegtyp,lai,veg,           &
     tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,qvsfc,windsp,psfc,  &
     rhoa,precip,tair,qvair, cdha,cdqa,radsw, rnflx,shflx,lhflx,  &
     gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34,              &
     frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief)  

!------------------------------------------------------------------
!                                                                  
!  AUTHOR: Diandong Ren (dd_ren@rossby.metr.ou.edu)                
!          Ming Xue     (mxue@ou.edu)
!  12/08/2001                                                      
!                                                                  
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!------------------------------------------------------------------
!                                                                  
!  PURPOSE:                                                        
!  To calculate the right hand side forcing terms in the 
!  soil-vegetation model.
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz
!  REAL :: eta_fac,wmax,c1max,sigma_sqr

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny)        ! Surface soil moisture
  REAL :: wetdp  (nx,ny)        ! Deep soil moisture
  REAL :: wetcanp(nx,ny)        ! Canopy water amount
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Surface air specific humidity (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !

  REAL :: frc_tsfc(nx,ny)      ! Right hand side forcing for tsfc eq.
  REAL :: frc_tdeep(nx,ny)       ! Right hand side forcing for tsoil eq.
  REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq. 
  REAL :: frc_wdp(nx,ny)        ! Right hand side forcing for wetsp eq. 
  REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.
  REAL :: wrmax(nx, ny)         ! Maximum value for canopy moisture, wetcanp

  REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
  REAL :: c2wg        ! Coefficient in the surface moisture Eq. of Wg

  REAL :: pi          ! Pi

  PARAMETER (pi = 3.141592654)

!  REAL :: tau         ! Seconds of a day = 24. * 3600.
!  PARAMETER (tau = 86400.)

  REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.
  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun
  REAL :: cg          ! Soil thermal coefficient for bare ground
  REAL :: wgeq        ! Surface moisture when gravity balances the capillarity
  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation

  INTEGER :: i, j

  REAL :: tema
  REAL :: eta, c1max, wmax, sig2 


  REAL :: relief       ! Difference between seasonal average skin and deep soil 
                      ! temperature  (t_skin - t_deeplayer).
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  log100    = ALOG(100.0)
  dtsfc2    = dtsfc/2.0

  DO j = 1, ny-1
    DO i = 1, nx-1
      tema = MAX( wetdp(i,j), wwlt(soiltyp(i,j)) )
  
      IF (snowdpth(i,j) >= snowdepth_crit) THEN !snow cover
        ct(i,j) = cg_snow
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                     &
                          *snowflxfac/(tau*ct(i,j))
      ELSE
        cg = cgsat(soiltyp(i,j))                                      &
            * ( wsat(soiltyp(i,j))/tema )                             &
            **( bslope(soiltyp(i,j))/log100 )
        ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv + veg(i,j) * cg )     
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
      END IF
  
      shflx(i,j) = rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)           &
                  * ( tsfc(i,j) - tair(i,j) )  
    END DO
  END DO

  IF ( moist == 0 ) THEN

    DO j = 1, ny-1
      DO i = 1, nx-1
         evaprg(i,j) = 0.0    
         evaprtr(i,j) = 0.0   
         evaprr(i,j) = 0.0  
         lhflx(i,j) = 0.0   
      END DO
    END DO

  ELSE

    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)

    CALL evapflx(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,         &
                   soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp,wetcanp,    &
                   snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat)
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
      evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
      evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
      evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
      lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)
      lhflx(i,j) = lhflx(i,j) * lathv       ! Latent heat flux
    END DO
  END DO

  DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        frc_tsfc (i,j) = 0.
        frc_tdeep(i,j) = 0.
      ELSE
        frc_tsfc(i,j) = ct(i,j)*(rnflx(i,j)-shflx(i,j)-lhflx(i,j)-gflx(i,j))
        IF ( snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_tdeep(i,j)= 1.03* (tsfc(i,j)-tdeep(i,j)-relief)*tauinv*snowflxfac ! Snow cover
        ELSE
!          frc_tdeep(i,j)= 1.03*(tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
           frc_tdeep(i,j)= (tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
        END IF
      END IF
    END DO
  END DO

  IF ( moist /= 0 ) THEN   

    DO j = 1, ny-1                                       !HYDROLOGY
      DO i = 1, nx-1
        wrmax(i,j) = (0.2 * 1.e-3 * veg(i,j) * lai(i,j))/dr     ! meter
        IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.           &
             snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_wsfc(i,j) = 0.
          frc_wdp(i,j)  = 0.
          frc_wcnp(i,j) = 0.
        ELSE

          tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

          c1wg = 0.4* c1sat(soiltyp(i,j))                               &
                 * ( wsat(soiltyp(i,j)) / tema )                        &
                 **( bslope(soiltyp(i,j)) / 2.0 + 1.0)
!--------------------------------------------------------------------------
! Replacement Cl to improve dry soils (NM1996 - A.3)   (JAB)
!--------------------------------------------------------------------------

          IF (wetsfc(i,j) < wwlt(soiltyp(i,j)) ) THEN

            eta = (-0.01815 * tsfc(i,j) + 6.41 ) * wwlt(soiltyp(i,j)) + &
                  ( 0.0065 * tsfc(i,j) - 1.4 )
      
            c1max = (1.19 * wwlt(soiltyp(i,j)) - 5.09)*0.01*tsfc(i,j)   &
                   +( 1.464 * wwlt(soiltyp(i,j)) + 17.86 )
      
            wmax = eta * wwlt(soiltyp(i,j))
      
            sig2 = - ( (2 * alog ( 0.01 / c1max )) / ( wmax*wmax) )
      
            c1wg = c1max * EXP ( -0.5* ( (wetsfc(i,j)-wmax)*(wetsfc(i,j)-wmax)*sig2 ) )
          ENDIF
!          if(i==sx .and. j==sy) then                 
!             write(6,'(a,f16.6)') 'c1wg =',c1wg
!          endif             

!----------------------------------------------------------------------------

          c2wg = c2ref(soiltyp(i,j)) * wetdp(i,j)                       &
                / ( wsat(soiltyp(i,j)) - wetdp(i,j) + wetsml )
          wgeq = wetdp(i,j) - wsat(soiltyp(i,j))                        &
                  * awgeq(soiltyp(i,j))                                 &
                  * ( wetdp(i,j) / wsat(soiltyp(i,j)) )                 &
                  ** pwgeq(soiltyp(i,j))                                &
                  * ( 1.0 - ( wetdp(i,j) / wsat(soiltyp(i,j)) )         &
                  ** ( 8 * pwgeq(soiltyp(i,j)) ) )

          frc_wsfc(i,j) = c1wg * ( ( 1.0 - veg(i,j) )                   &
                                  * precip(i,j) - evaprg(i,j) )         &
                          /(rhow * d1)                                  &
                         - c2wg * ( wetsfc(i,j) - wgeq ) * tauinv

 
          frc_wdp(i,j) = ( ( 1.0 - veg(i,j) ) * precip(i,j)             &
                        - evaprg(i,j) - evaprtr(i,j) )                  &
                      / ( rhow * d2 )     

          wr2max = ( wrmax(i,j) - wetcanp(i,j) ) / dtsfc2    
          vegp = veg(i,j) * precip(i,j)
          pnet = vegp - evaprr(i,j)
          tema = pnet - wr2max * rhow    ! aqui nao teveria ser (rhow*dr)????
          runoff = MAX( tema, 0.0 )
          vegp = vegp - runoff

          frc_wcnp(i,j) = ( vegp - evaprr(i,j) ) / rhow 

        END IF
!        if(i==sx .and. j==sy) then                 
!           write(6,'(a,f16.6,a,f16.6)') 'c2wg  =',c2wg,'wgeq =',wgeq
!           write(6,'(a,f16.6,a,f16.6)') 'runoff=',runoff,'wgeq= ',wgeq
!        endif             
        

      END DO
    END DO           

  END IF

  RETURN
END SUBROUTINE soilebm_frc


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE EVAPFLX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE evapflx(nx,ny, radsw, f34, cdqa,                             &
           windsp,psfc,rhoa,qvair,                                      &
           soiltyp,vegtyp,lai,veg,                                      &
           tsfc,wetsfc,wetdp,wetcanp,snowdpth,                          &
           evaprg,evaprtr,evaprr,qvflx,qvsat)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate:
!
!  1. Evaporation from ground surface,
!
!        evaprg = rhoa * cdq * windsp * evaprg'
!
!  where
!
!        evaprg' = (1.-veg) * (rhgs * qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  2. Direct evaporation from the fraction delta of the foliage covered
!     by intercepted water.
!
!        Er = rhoa * cdq * windsp * Er'
!
!        Er' = delta * veg * (qvsats-qvair)
!
!  3. Transpiration of the remaining part (1-delta) of leaves,
!
!        Etr = rhoa * cdq * windsp * Etr'
!
!  where
!
!        Etr' = veg * (1-delta) * Ra/(Ra+Rs) * (qvsats-qvair )
!
!  and Ra is aerodynamic resistance and Rs is the surface resistance
!
!                     1
!        Ra = ----------------
!               cdq * windsp
!
!                   Rsmin
!        Rs = ------------------
!              LAI*F1*F2*F3*F4
!
!
!               f + Rsmin/Rsmax
!        F1 = -------------------
!                   f + 1
!
!                   Rg    2
!        f = 0.55 ----- -----
!                  Rgl   LAI
!
!             -  1,                             Wfc < W2
!             |
!             |    W2 - Wwlt
!        F2 = -  ------------,              Wwlt <= W2 <= Wfc
!             |   Wfc - Wwlt
!             |
!             -  0,                             W2 < Wwlt
!
!
!               1-0.06*(qvsats-qvair),   qvsats-qvair <= 12.5 g/kg
!        F3 = {
!               0.25,                      otherwise
!
!
!        F4 = 1 - 0.0016 * (298-tair)**2
!
!  4. Water vapor flux, qvflx,
!
!        qvflx = rhoa * cdq * windsp * qvflx'
!
!        qvflx' = (Eg' + Etr' + Er') = (qvsfc - qvair)
!
!     where qvsfc is the effective surface specific humidity
!
!        (qvsfc - qvair) = (Eg' + Etr' + Er')
!
!  This subroutine will solve Eg', Etr', Er', and qvflx'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  4/20/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    radsw    Incoming solar radiation
!
!    f34      f3*f4;
!             f3 = Fractional conductance of atmospheric vapor pressure
!             f4 = Fractional conductance of air temperature
!
!    cdqa     Array for surface drag coefficient for water vapor
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    tsfc     Surface soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Vegetation moisture
!
!    psfc     Surface pressure
!    qvair    Specific humidity near the surface
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!
!  OUTPUT:
!
!    evaprp   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!    qvflx    Water vapor flux
!
!  WORK ARRAY:
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny
!
  REAL :: radsw  (nx,ny)

  REAL :: f34    (nx,ny)

  REAL :: cdqa   (nx,ny)

  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp (nx,ny)
  REAL :: lai    (nx,ny)
  REAL :: veg    (nx,ny)

  REAL :: tsfc   (nx,ny)
  REAL :: wetsfc (nx,ny)
  REAL :: wetdp  (nx,ny)
  REAL :: wetcanp(nx,ny)
  REAL :: snowdpth(nx,ny)
!
  REAL :: psfc   (nx,ny)
  REAL :: qvair  (nx,ny)
  REAL :: windsp (nx,ny)
  REAL :: rhoa   (nx,ny)

  REAL :: evaprg (nx,ny)
  REAL :: evaprtr(nx,ny)
  REAL :: evaprr (nx,ny)
  REAL :: qvflx  (nx,ny)
  REAL :: qvsat  (nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654 )
!
  INTEGER :: i, j

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp
  REAL :: rstcoef     ! Coefficient of resistance
  REAL :: delta
  REAL :: rhgs        ! R.H. at ground surface
!
  REAL :: tema
  REAL :: temb

  REAL :: pterm
  REAL :: mterm
  REAL :: waf         ! Available soil moisture fraction
  REAL :: bw          ! Half point of the available soil mstr frctn curve
  REAL :: ps          ! Shelter factor 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j = 1, ny-1
    DO i = 1, nx-1

      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
!
!-----------------------------------------------------------------------
!
!  Over water and ice, qvsfc should be saturated
!
!-----------------------------------------------------------------------
!
        qvflx(i,j) = qvsat(i,j) - qvair(i,j)

        evaprg (i,j) = qvflx(i,j)
        evaprtr(i,j) = 0.0
        evaprr (i,j) = 0.0

      ELSE                                  ! over land

!        wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! in meter - original do soilebm
	wrmax = 0.2 * veg(i,j) * lai(i,j) / rhow
!        if (i==sx .and. j==sy) then
!           write(6,'(a,3f22.12,2i4)') 'wrmax=', wrmax,veg(i,j),lai(i,j),vegtyp(i,j),soiltyp(i,j)
!	endif
!
!-----------------------------------------------------------------------
!
!  In order to calculate the qv flux at the surface, we need to
!  calculate some parameters like the resistance coefficient,
!
!      rstcoef = rsta / (rsts + rsta)
!
!  where rsta is the aerodynamic resistance and rsts is the surface
!  resistances.
!
!      rsta = 1 / ( cdq * Va )
!
!      rsts = rsmin(vtyp) / ( lai(vtyp) * f1 * f2 * f3 * f4 )
!
!  f3 * f4 is time-independent and has been calculated previously
!  and stored in f34(i,j)
!
!-----------------------------------------------------------------------
!
!  Calculate f1.
!
!        f = 0.55*(radsw/rgl(vtyp))*(2./lai(vtyp))       ! NP, Eq. 34
!        f1 = ( rsmin(vtyp)/rsmax + f ) / ( 1. + f )     ! NP, Eq. 34
!
!  Note: the incoming solar radiation radsw is stored in radsw(i,j).
!
!-----------------------------------------------------------------------
!
        IF ( lai(i,j) == 0. ) THEN        ! No vegetation, I/W
          rstcoef = 1.
        ELSE
!          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
!               * ( 2.0 / lai(i,j) )

          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
               * ( 2.0 )       !(XP2001 - Eq.8)   (JAB) 

          rstcoef = ( rsmin(vegtyp(i,j)) / rsmax + temb )               &
                  / ( 1.0 + temb )
        END IF
!
!-----------------------------------------------------------------------
!
!  Calculate f2 and f1*f2.
!
!-----------------------------------------------------------------------
!
!        pterm = .5 + SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!        mterm = SIGN( .5, wetdp(i,j) - wwlt(soiltyp(i,j)) )              &
!              - SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!
!        rstcoef = rstcoef * ( pterm + mterm                             &
!                * ( wetdp(i,j)-wwlt(soiltyp(i,j)) )                     &
!                / ( wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) ) )

        waf = ( wetdp(i,j)-wwlt(soiltyp(i,j)) )            &
                / ( wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )

        bw = wwlt(soiltyp(i,j)) + (wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) )  &
                / 3.0

        rstcoef = rstcoef / ( 1.0 + EXP( -5.0 * (waf - bw) ) )
!             (XP2001 - Eq.9)                                       (JAB)
!
!-----------------------------------------------------------------------
!
!  Calculate lai*f1*f2*f3*f4 where f3*f4 is stored in f34(i,j).
!
!-----------------------------------------------------------------------
!
!        rstcoef = lai(i,j)*rstcoef*f34(i,j)      ! lai*f1*f2*f3*f4

        ps = 0.3 * lai(i,j) + 0.7                 ! XP2001             (JAB)
        rstcoef = lai(i,j)*rstcoef*f34(i,j)/ps    ! lai*f1*f2*f3*f4/ps (JAB)
!
!-----------------------------------------------------------------------
!
!  Calculate the resistance coefficient, rsta/(rsts+rsta)
!
!        rsts = rsmin(vtyp)/(lai(i,j)*f1*f2*f3*f4) ! Sfc. resistance
!        rsta = 1./(cdh*va)                 ! NP, between Eq. 32 & 33
!        rstcoef = rsta/(rsta+rsts)
!                = 1/(1+rsts/rsta)
!
!-----------------------------------------------------------------------
!
        tema = rsmin(vegtyp(i,j)) * cdqa(i,j) * windsp(i,j)

        IF ( ABS(rstcoef) > 1.0E-30 ) THEN
          rstcoef = 1.0 / (1.0 + tema/rstcoef)
        END IF
!
!-----------------------------------------------------------------------
!
!  1. evaprg'
!
!        evaprg' = (1.-veg) * (rhgs*qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  evaprg will be stored for current and future use to
!  calculate the latent heat flux and soil moisture transports.
!
!-----------------------------------------------------------------------
!
        pterm = .5 + SIGN( .5, wetsfc(i,j)-1.1*wfc(soiltyp(i,j)) )

        rhgs = pterm                                                    &
             + (1.-pterm) * ( 0.25 * ( 1.0 - COS( wetsfc(i,j)           &
                              * pi / (1.1*wfc(soiltyp(i,j))))) ** 2 )

!      IF (snowdpth(i,j) .ge. snowdepth_crit) rhgs=1.0 !Snow cover

!        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
!                    * ( rhgs * qvsat(i,j) - qvair(i,j) )

        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
                    *  rhgs * ( qvsat(i,j) - qvair(i,j) )  !XP2001 (JAB) 

!
!-----------------------------------------------------------------------
!
!  2. Transpiration of the remaining part (1-delta) of leaves, Etr',
!
!        Etr' = (1-delta) * veg
!             * Ra/(Ra+Rs) * ( qvsats - qvair )
!
!  3. Direct evaporation from the fraction delta, Er'
!
!        Er' = delta * veg * ( qvsats - qvair )
!
!  Er' and Etr' are stored for future use to calculate the latent heat
!  flux and soil moisture transports.
!
!-----------------------------------------------------------------------
!
!        IF ( wrmax == 0.0 ) THEN
        IF ( wrmax .lt. 0.000000001 ) THEN
          delta = 0.0
!          write(6,'(a,f18.8)') 'delta0=', delta
        ELSE
          delta = ( wetcanp(i,j) / wrmax ) ** 0.66666667
!          if (i==sx .and. j==sy) then
!             write(6,'(a,3f18.14)') 'delta1=', delta, wetcanp(i,j), wrmax
!	  endif
        END IF

        pterm = .5 + SIGN( .5, qvsat(i,j) - qvair(i,j) )
        delta = pterm * delta + ( 1. - pterm )

        tema = veg(i,j) * ( qvsat(i,j) - qvair(i,j) )

        evaprtr(i,j) = ( 1.0 - delta ) * rstcoef * tema
        evaprr (i,j) = delta * tema
!
!-----------------------------------------------------------------------
!
!  4. Water vapor flux, qvflx',
!
!        qvflx' = evaprg' + evaprtr' + evaprr'
!
!  qvflx will be saved for the future use to calculate the latent
!  heat flux.
!
!-----------------------------------------------------------------------
!
        qvflx(i,j) = evaprg(i,j) + evaprtr(i,j) + evaprr(i,j)
      END IF
                                                   ! NP, expl. Eq. 26-27
    END DO
  END DO

  RETURN
END SUBROUTINE evapflx  



!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRIDIAG2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tridiag2(nx,ny,nz,ibgn,iend,jbgn,jend,                       &
           kbgn,kend,a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the solution for a tridiagonal equations.
!  Note:
!       a: left of main diagonal, useful range [kbgn+1,kend];
!       b: main diagonal, useful range [kbgn,kend];
!       c: right of main diagonal, useful range [kbgn,kend-1];
!       d: right hand side of equations.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: X. Song, Donghai Wang and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where operation begins.
!    iend     i-index where operation ends.
!    jbgn     j-index where operation begins.
!    jend     j-index where operation ends.
!    kbgn     k-index where operation begins.
!    kend     k-index where operation ends.
!
!    a        left of main diagonal, useful range [kbgn+1,kend]
!    b        main diagonal, useful range [kbgn,kend]
!    c        right of main diagonal, useful range [kbgn,kend-1]
!    d        right hand side of equations
!
!  OUTPUT:
!
!    d        The solution
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: c(nx,ny,nz)       ! Right of main diagonal
  REAL :: a(nx,ny,nz)       ! Left of main diagonal
  REAL :: d(nx,ny,nz)       ! Right hand side of equations
  REAL :: b(nx,ny,nz)       ! Main diagonal

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: i,j,k
  REAL :: r
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn+1,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        r=a(i,j,k)/b(i,j,k-1)
        b(i,j,k)=b(i,j,k)-r*c(i,j,k-1)
        d(i,j,k)=d(i,j,k)-r*d(i,j,k-1)
      END DO
    END DO
  END DO

  DO j=jbgn,jend
    DO i=ibgn,iend
      d(i,j,kend)=d(i,j,kend)/b(i,j,kend)
    END DO
  END DO

  DO k=kend-1,kbgn,-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        d(i,j,k)=(d(i,j,k)-c(i,j,k)*d(i,j,k+1))/b(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tridiag2
