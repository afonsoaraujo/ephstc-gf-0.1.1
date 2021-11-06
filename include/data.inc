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



!DEC$ IF .not. DEFINED (DATA_INC)
!DEC$ DEFINE DATA_INC

   include 'settings.inc'

type OPTIONS
   integer          ::   SlopeAspect
   integer          ::   teste1
end type OPTIONS


type CALENDAR
   integer(4)       ::   ano
   integer(4)       ::   mes
   integer(4)       ::   diam
   integer(4)       ::   diaj     ! Dia juliano
   integer(4)       ::   hora
   integer(4)       ::   min
end type CALENDAR


type FILENAME
   character(128)   ::   Maskfl      ! Arquivo contendo a máscara da bacia
   character(128)   ::   AirTempfl   ! Arquivo contendo forçante de temperatura do ar
   character(128)   ::   AirPressfl  ! Arquivo contendo forçante de pressão do ar
   character(128)   ::   Relhfl      ! Arquivo contendo forçante de umidade relativa do ar
   character(128)   ::   Wspdfl      ! Arquivo contendo forçante de velocidade do vento
   character(128)   ::   Rsifl       ! Arquivo contendo forçante de radiação solar insidente
   character(128)   ::   Rsrfl       ! Arquivo contendo forçante de radiação solar refletida
   character(128)   ::   Precipfl    ! Arquivo contendo forçante de precipitação
   character(128)   ::   Topofl      ! Arquivo de topografia
   character(128)   ::   Slopefl     ! Arquivo de Declividades (Slope)
   character(128)   ::   Aspectfl    ! Arquivo de Aspecto (Aspect)
   character(128)   ::   Soilfl      ! Arquivo de Tipos de solos
   character(128)   ::   SDepthfl    ! Arquivo com as profundidades do solo
   character(128)   ::   WTDepthfl   ! Arquivo com a profundidade do lençol
   character(128)   ::   STblfl      ! Tabela das propriedades do solo
   character(128)   ::   Vegfl       ! Arquivo com tipos de vegetação
   character(128)   ::   NDVIfl      ! Arquivo contendo o índice de vegetação (NDVI)
   character(128)   ::   VTblfl      ! Tabela das propriedades da vegetação
end type FILENAME


type MAPSIZE
   character   ::   System(BUFSIZE)   ! Coordinate system
   real        ::   Xorig             ! X coordinate of Northwest corner
   real        ::   Yorig             ! Y coordinate of Northwest corner
   integer     ::   X                 ! Current x position
   integer     ::   Y                 ! Current y position
   integer     ::   NX                ! Number of pixels in x direction
   integer     ::   NY                ! Number of pixels in y direction
   real        ::   DX                ! Pixel spacing in x-direction
   real        ::   DY                ! Pixel spacing in y-direction
   real        ::   DXY               ! Pixel spacing in diagonal
   integer     ::   OffsetX           ! Offset in x-direction compared to basemap
   integer     ::   OffsetY           ! Offset in y-direction compared to basemap
end type MAPSIZE


type FORCPIX
!   integer     ::  NX_Forc           !
!   integer     ::  NY_Forc           !
   real        ::   AirTemp           ! Air Temperature at surface      -  [C]
   real        ::   Press             ! Air Pressure                    -  [hPa]
   real        ::   Relh              ! Relative humidity               -  [%]
   real        ::   Wspd              ! Wind speed                      -  [m/s]
   real        ::   Rsi               ! Insident solar radiation (sw)   -  [W/m2]  
   real        ::   Rsr               ! Reflected solar readiation (lw) -  [W/m2]  
   real        ::   Rsn               ! Net radiation at surface        -  [W/m2]
   real        ::   Precip            ! Precipitation                   -  [mm/h]
end type FORCPIX


type TOPOPIX
   real        ::   Dem               ! Elevations
   integer     ::   Mask              ! Mask for modeled area (0 - out; 1 - in)
   real        ::   Travel            ! Travel time
   real        ::   Grad              ! Sum of downslope slope-width products
   real        ::   Slope             ! Land surface slope
   real        ::   Aspect            ! Land surface slope direction
   real        ::   FlowGrad          ! Magnitude of subsurface flow gradient slope * width
   real        ::   Dir(NDIRS)        ! Fraction of flux moving in each direction
   integer     ::   TotalDir          ! Sum of Dir array
end type TOPOPIX


type VEGTABLE
   character          ::   Desc(BUFSIZE)        ! Vegetation type
   integer            ::   Index                !
   integer            ::   NVegLayers           ! Number of vegetation layers
   integer            ::   NSoilLayers          ! Number of soil layers
!  unsigned char      ::   OverStory            ! TRUE if there is an overstory
!  unsigned char      ::   UnderStory           ! TRUE if there is an understory
   real,allocatable   ::   Adjust(:)
   real,allocatable   ::   Height(:)            ! Height of vegetation (in m)
   real,allocatable   ::   VegFrac(:)           ! Fractional coverage
   real,allocatable   ::   HemiFrac(:)          ! used to calculated longwave radiation balance
   real               ::   LAI                  ! One Sided Leaf Area Index
   real,allocatable   ::   Roufns(:)            ! Surface roughness 
   real               ::   Wrmax                ! Maximum value for canopy moisture, wetcanp
   real,allocatable   ::   MaxInt(:)            ! Maximum interception storage (m)
   real,allocatable   ::   RsMax(:)             ! Maximum stomatal resistance 
   real,allocatable   ::   RsMin(:)             ! Minimum stomatal resistance
   real,allocatable   ::   Rgl(:)               ! Calculo de f(2)
   real,allocatable   ::   MoistThres(:)        ! Soil moisture threshold above which soil moisture 
                                                !    does not restrict transpiration
                                                !   stomatal closure occurs (Pa)
   real,allocatable   ::   RootFract(:,:)       ! Fraction of roots in each soil layer
   real,allocatable   ::   RootDepth(:)         ! Depth of root zones
                                                !    when the "canopy radiation attenuation" 
                                                !    option is set to fixed
   real               ::   TotalDepth           ! total depth of all root zones
   real               ::   LeafAngleA           ! parameter describing the Leaf Angle Distribution
   real               ::   LeafAngleB           ! parameter describing the leaf Angle Distribution
   real               ::   Scat                 ! scattering parameter (between 0.7 and 0.85)
   real,allocatable   ::   Albedo(:)            ! Albedo for each vegetation layer
   real               ::   Cn                   ! Canopy attenuation coefficient for wind profile
   real               ::   MaxSnowInt           ! Maximum snow interception capacity for the overstory
   real               ::   ImpervFrac           ! fraction of pixel that is impervious 
   real               ::   Ra(2)                ! Aerodynamic resistance in the absence of snow
   real               ::   RaSnow               ! Aerodynamic resistance for the lower boundary in the presence of snow
   real               ::   Trunk                ! Fraction of overstory height that identifies the top of the trunk space
   real               ::   U(2)                 ! Wind speed profile (m/s)
   real               ::   USnow                ! wind speed 2, above snow surface (m/s)
end type VEGTABLE


type VEGPIX
   integer            ::   Veg                   ! Vegetation type
   real               ::   Tcanopy               ! Canopy temperature (C)
   real               ::   NDVI                  ! Índice de vegetação normalizado 
end type VEGPIX


type SOILTABLE
   character          ::   Desc(BUFSIZE)         ! Soil type
   integer            ::   Index
   integer            ::   NLayers               ! Number of soil layers 
   real               ::   Albedo                ! Albedo of the soil surface 
   real,allocatable   ::   Porosity(:)           ! Soil porosity for each layer
   real,allocatable   ::   PoreDist(:)           ! Pore size distribution for each layer
   real,allocatable   ::   Press(:)              ! Soil bubling pressure for each layer
   real,allocatable   ::   Wsat(:)               ! Saturated volumetric moisture content. JN, Tab. 1
   real,allocatable   ::   Wfcap(:)              ! Field capacity for each layer
   real,allocatable   ::   Wwilt(:)              ! Wilting point for each layer
   real,allocatable   ::   Dens(:)               ! Soil density (kg/m^3)
   real,allocatable   ::   Ks(:)                 ! Saturated hydraulic conductivity (vertical) for each layer
   real               ::   KsLat                 ! Saturated hydraulic conductivity (lateral)
   real               ::   KsLatExp              ! Exponent for vertical change of KsLat
   real,allocatable   ::   KhDry(:)              ! Thermal conductivity for dry soil     (W/(m*K))
   real,allocatable   ::   KhSol(:)              ! Effective solids thermal conductivity (W/(M*K))
   real,allocatable   ::   Ch(:)                 ! Heat capacity for soil medium
   real               ::   Cgsat                 ! Coeficiente térmico de saturação para o solo
   real               ::   C1sat                 ! Coeficiente C1 solo na saturação (W=Wsat)
   real               ::   C2ref                 ! Coeficiente C2 solo semi-saturado (W=0.5 Wsat)
   real               ::   C1wilt                ! Coeficiente C1 solo ponto de murchamento (W=Wwilt) 
   real               ::   C4b                   ! Parâmetro de textura do solo
   real               ::   C4lref                ! Coeficiente de ajuste de grid em função do tipo de solo
   real               ::   a                     !
   real               ::   p                     !  
   real               ::   MaxInfiltrationRate   ! Maximum infiltration rate for upper layer (m/s)
end type SOILTABLE


type SOILPIX
   integer            ::   Soil                  ! Soil type
   real               ::   Depth                 ! Depth of total soil zone, including all root
                                                 ! zone layers, and the saturated zone
   real,allocatable   ::   Moist(:)              ! Soil moisture content in layers
   real,allocatable   ::   Perc(:)               ! Percolation from layers
   real,allocatable   ::   Temp(:)               ! Temperature in each layer (C)
   real               ::   TableDepth            ! Depth of water table below ground surface (m)
   real               ::   WaterLevel            ! Absolute height of the watertable above datum (m), i.e. 
                                                 ! corrected for terrain elevation
   real               ::   SatFlow               ! amount of saturated flow generated
   real               ::   Runoff                ! amount of surface runoff generated from HOF and Return flow
   real               ::   ChannelInt            ! amount of subsurface flow intercepted by the channel
   real               ::   RoadInt               ! amount of water intercepted by the road
   real               ::   TSurf                 ! Soil surface temperature
   real               ::   Qnet                  ! Net radiation exchange at surface
   real               ::   Qrest                 ! Rest term for energy balance (should be 0)
   real               ::   Qs                    ! Sensible heat exchange
   real               ::   Qe                    ! Latent heat exchange 
   real               ::   Qg                    ! Ground heat exchange
   real               ::   Qst                   ! Ground heat storage
   real               ::   Ra                    ! Soil surface aerodynamic resistance (s/m)
   real               ::   SurfaceWater          ! used in the impervious calculations (m)
end type SOILPIX


type SVATPIX
   real           ::   qvsatta
   real           ::   qvsatts
   real           ::   qva
   real           ::   ea
   real           ::   esatta
   real           ::   esatts
   real           ::   hu
   real           ::   Rg
   real           ::   LE
   real           ::   H
   real           ::   Rn
   real           ::   F1
   real           ::   F2
   real           ::   F3
   real           ::   F4
   real           ::   f
   real           ::   Pg
   real           ::   Pr
   real           ::   Eg
   real           ::   Etr
   real           ::   Er
   real           ::   Ra
   real           ::   Rs
   real           ::   Fw
   real           ::   G
   real           ::   LEg
   real           ::   LEr
   real           ::   LEtr
end type SVATPIX


type SVATABLE
   real           ::   aw
   real           ::   bw
   real           ::   Ct         ! Capacidade de calor da superfície
   real           ::   Cg         ! Coeficiente térmico para o solo
   real           ::   Cv         ! Coeficiente térmico para a vegetação
   real(4)        ::   Cdh        ! Coeficiente de transferência de calor 
   real(4)        ::   Cdq        ! Coeficiente de transferência de umidade
   real           ::   C1         ! Parâmetro de propriedade hidráulica
   real           ::   C1max      ! Valor máximo para C1
   real           ::   C2         ! Parâmetro de equilíbrio do perfil de umidade
   real           ::   C3         ! Parâmetro de restauração newtoniana
   real           ::   C4         ! Parâmetro de equilíbrio na interface
   real           ::   C4l        ! Combinação de C4b e C4lref
   real           ::   Wl         ! Limite numérico de W para C2
   real           ::   Wmax       ! Limite máximo para W
   real           ::   Wgeq       ! Umidade de equilíbrio das forças gravitacionais [m3 m-3]
   real           ::   W3s        ! Umidade média de restauração newtoniana
   real           ::   W23        ! Umidade na interface da zona de raízes
   real           ::   d1         ! Profundidade da camada superficial do solo
   real           ::   d2         ! Profundidade da zona de raízes
   real           ::   d3         ! Profundidade da zona sub-raízes
   real           ::   dr         !
   real           ::   K2         !
   real           ::   K3         !
   real           ::   D12        !
   real           ::   D23        !
   real           ::   tau        ! Constante de restauração diária
   real           ::   sigma
   real           ::   sigma2
   real           ::   emissa
   real           ::   emissg
   real           ::   roa
   real           ::   row
   real           ::   pi
   real           ::   L
   real           ::   Cp
   real           ::   a          !
   real           ::   p          !
end type SVATABLE

!DEC$ ENDIF
