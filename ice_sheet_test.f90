!
!----------------------------------------------------------
!   The Globally Resolved Energy Balance (GREB) Model
!----------------------------------------------------------
!
!   Authors: Dietmar Dommenget, Janine Flöter, Tobias Bayr and Christian Stassen
!
!   References:	- Conceptual Understanding of Climate Change with a
!              	Globally Resolved Energy Balance Model
!               by Dietmar Dommenget and Janine Flöter, J. Clim Dyn (2011) 37: 2143.
!               doi:10.1007/s00382-011-1026-0

!               - A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0
!               by Christian Stassen, Dietmar Dommenget & Nicholas Loveday.
!               Geosci. Model Dev., 12, 425-440, https://doi.org/10.5194/gmd-12-425-2019, 2019.
!
!		            - The Monash Simple Climate Model Experiments: An interactive database
!		            of the mean climate, climate change and scenarios simulations
!		            by Dietmar Dommenget, Kerry Nice, Tobias Bayr, Dieter Kasang, Christian Stassen
!		            and Mike Rezny, submitted to Geoscientific Model Development
!
!
!  input fields: The GREB model needs the following fields to be specified before
!                the main subroutine greb_model is called:
!
!    Tclim(xdim,ydim,nstep_yr):   mean Tsurf                       [K]
!    uclim(xdim,ydim,nstep_yr):   mean zonal wind speed            [m/s]
!    vclim(xdim,ydim,nstep_yr):   mean meridional wind speed       [m/s]
!    qclim(xdim,ydim,nstep_yr):   mean atmospheric humidity        [kg/kg]
!  cldclim(xdim,ydim,nstep_yr):   total cloud cover                [0-1]
! swetclim(xdim,ydim,nstep_yr):   soil wetness, fraction of total  [0-1]
!   Toclim(xdim,ydim,nstep_yr):   mean deep ocean temperature      [K]
!  mldclim(xdim,ydim,nstep_yr):   mean ocean mixed layer depth     [m]
!   z_topo(xdim,ydim):            topography (<0 are ocean points) [m]
!  glacier(xdim,ydim):            glacier mask ( >0.5 are glacier points )
! sw_solar(ydim,nstep_yr):        24hrs mean solar radiation       [W/m^2]
!
!
! possible experiments:
!
!  log_exp = 1  deconstruct mean climate
!  		you can switch on or off climate components as given in the
!		module physics in 'deconstruct mean state switches' section
!
!  log_exp = 10  deconstruct 2xCO2 response
!  		 you can switch on or off climate components as given in the
!		 module physics in 'deconstruct 2xco2 switches' section
!
!  log_exp = 20  2x   CO2
!  log_exp = 21  4x   CO2
!  log_exp = 22  10x  CO2
!  log_exp = 23  0.5x CO2
!  log_exp = 24  0x   CO2
!
!  log_exp = 25  CO2-wave 30yrs-period
!  log_exp = 26  2xCO2 30yrs followed by 70yrs CO2-ctrl
!  log_exp = 27  solar constant +27W/m2 (~2xCO2 warming)
!  log_exp = 28  11yrs solar cycle
!
!  log_exp = 30  paleo solar 231Kyr before present & CO2=200ppm
!  log_exp = 31  paleo solar 231Kyr before present
!  log_exp = 32  paleo CO2=200ppm 231Kyr before present
!
!  log_exp = 35  solar radiation obliquity changes
!  log_exp = 36  solar radiation eccentricity changes
!  log_exp = 37  solar radiation radius changes
!
!  log_exp = 40  partial 2xCO2 Northern hemisphere
!  log_exp = 41  partial 2xCO2 Southern hemisphere
!  log_exp = 42  partial 2xCO2 Tropics
!  log_exp = 43  partial 2xCO2 Extratropics
!  log_exp = 44  partial 2xCO2 Ocean
!  log_exp = 45  partial 2xCO2 Land
!  log_exp = 46  partial 2xCO2 Boreal Winter
!  log_exp = 47  partial 2xCO2 Boreal Summer
!
!  log_exp = 95  IPCC A1B scenario
!  log_exp = 96  IPCC RCP26 scenario
!  log_exp = 97  IPCC RCP45 scenario
!  log_exp = 98  IPCC RCP60 scenario
!  log_exp = 99  IPCC RCP85 scenario
!
!  log_exp = 230 run a climate change experiment with forced boundary conditions
!            (surface temperature, hodrizontal winds and omega) of the CMIP5
!            rcp85 ensemble mean response
!
!  log_exp = 240 & 241 run a El Nino (La Nina) experiment with forced boundary conditions
!            (surface temperature, hodrizontal winds and omega) of the ERA-Interim
!            composite mean response
!
!  log_exp = 310 GREB with coupled ice sheet model
!
!  log_exp = 100 run model with your own CO2 input file
!
!+++++++++++++++++++++++++++++++++++++++
module mo_numerics
!+++++++++++++++++++++++++++++++++++++++

! numerical parameter
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: dt        = 12*3600           ! time step [s]
  integer, parameter :: dt_crcl   = 0.5*3600          ! time step circulation [s]
  integer, parameter :: ndt_days  = 24*3600/dt        ! number of timesteps per day
  integer, parameter :: nstep_yr  = ndays_yr*ndt_days ! number of timesteps per year
  integer            :: time_flux = 0                 ! length of integration for flux correction [yrs]
  integer            :: time_ctrl = 0                 ! length of integration for control run  [yrs]
  integer            :: time_scnr = 0                 ! length of integration for scenario run [yrs]
  integer            :: ipx       = 1                 ! points for diagonstic print outs
  integer            :: ipy       = 1                 ! points for diagonstic print outs
  integer, parameter, dimension(12) :: jday_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat

  integer            :: ireal     = 1         	      ! record length for IO (machine dependent)
! 						        ireal = 4 for Mac Book Pro and Ubuntu Linux

  namelist / numerics / time_flux, time_ctrl, time_scnr

end module mo_numerics


!+++++++++++++++++++++++++++++++++++++++
module mo_physics
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  integer  :: log_exp   = 0              ! process control logics for expiments (see header)
! deconstruct mean state (dmc) switches
  integer  :: log_cloud_dmc   = 1              ! process control clouds
  integer  :: log_ocean_dmc   = 1              ! process control ocean
  integer  :: log_atmos_dmc   = 1              ! process control Atmosphere
  integer  :: log_co2_dmc     = 1              ! process control CO2
  integer  :: log_hydro_dmc   = 1              ! process control hydro
  integer  :: log_qflux_dmc   = 1              ! process control qflux corrections
! deconstruct 2xco2 (drsp) switches
  integer  :: log_topo_drsp   = 1              ! process control for topo
  integer  :: log_cloud_drsp  = 1              ! process control for clouds
  integer  :: log_humid_drsp  = 1              ! process control for humidity clim
  integer  :: log_ocean_drsp  = 1              ! process control for ocean
  integer  :: log_hydro_drsp  = 1              ! process control for hydro
! switches that are the same for both deconstructions
  integer  :: log_ice         = 1              ! process control ice-albedo
  integer  :: log_ice_sheet   = 1              ! process control for ice sheet
  integer  :: log_snow        = 1              ! process control snow-albedo (ice sheet)
  integer  :: log_hdif        = 1              ! process control Diffusion of heat
  integer  :: log_hadv        = 1              ! process control Advection of heat
  integer  :: log_vdif        = 1              ! process control Diffusion of vapor
  integer  :: log_vadv        = 1              ! process control Advection of vapor
! switches for the hydrological cycle
  integer  :: log_rain        = 0              ! process control precipitation parameterisation
  integer  :: log_eva         = 0              ! process control evaporation parameterisation
  integer  :: log_conv        = 0              ! process control advection parameterisation
  integer  :: log_clim        = 0              ! process control for reference climatology
! switches for external forcing files
  integer  :: log_tsurf_ext   = 0              ! process control evaporation parameterisation
  integer  :: log_hwind_ext   = 0              ! process control advection parameterisation
  integer  :: log_omega_ext   = 0              ! process control for reference climatology

! parameters for scenarios
  real     :: dradius   = 0.		 ! deviations from actual earth radius in %

! physical parameter (natural constants)
  parameter( pi        = 3.1416 )
  parameter( sig       = 5.6704e-8 )     ! stefan-boltzmann constant [W/m^2/K^4]
  parameter( rho_ocean = 999.1 )         ! density of water at T=15C [kg/m^2]
  parameter( rho_land  = 2600. )         ! density of solid rock [kg/m^2]
  parameter( rho_air   = 1.2 )           ! density of air at 20C at NN
  parameter( rho_ice   = 910. )           ! density of incompressible polycrystalline ice [kg/m^3]
  parameter( cp_ocean  = 4186. )         ! specific heat capacity of water at T=15C [J/kg/K]
  parameter( cp_land   = cp_ocean/4.5 )  ! specific heat capacity of dry land [J/kg/K]
  parameter( cp_air    = 1005. )         ! specific heat capacity of air      [J/kg/K]
  parameter( cp_ice    = 2009. )         ! specific heat capacity of ice [J/kg/K]
  parameter( eps       = 1. )            ! emissivity for IR
  real :: S0_var       = 100.            ! variation of solar constant   [%]

! physical parameter (model values)
  parameter( d_ocean   = 50. )                     ! depth of ocean column [m]
  parameter( d_land    = 2. )                      ! depth of land column  [m]
  parameter( d_air     = 5000. )                   ! depth of air column   [m]
  parameter( d_ice_max = 1. )                      ! maximum shortwave radiation penetarting depth of ice column   [m]
  parameter( cap_ocean = cp_ocean*rho_ocean )      ! heat capacity 1m ocean  [J/K/m^2]
  parameter( cap_land  = cp_land*rho_land*d_land ) ! heat capacity land   [J/K/m^2]
  parameter( cap_air   = cp_air*rho_air*d_air )    ! heat capacity air    [J/K/m^2]
  real :: ice_svlm  = 1./rho_ice                   ! specific volume of ice/snow [m3/kg]
  real :: ct_sens   = 22.5                         ! coupling for sensible heat
  real :: da_ice    = 0.25                         ! albedo diff for ice covered points
  real :: a_no_ice  = 0.1                          ! albedo for non-ice covered points
  real :: a_cloud   = 0.35                         ! albedo for clouds
  real :: a_snow    = 0.8                          ! process control snow-albedo (ice sheet)
  real :: Tl_ice1   = 273.15-10.                   ! temperature range of land snow-albedo feedback
  real :: Tl_ice2   = 273.15                       ! temperature range of land snow-albedo feedback
  real :: To_ice1   = 273.15-7.                    ! temperature range of ocean ice-albedo feedback
  real :: To_ice2   = 273.15-1.7                   ! temperature range of ocean ice-albedo feedback
  real :: co_turb   = 5.0                          ! turbolent mixing to deep ocean [W/K/m^2]
  real :: kappa     = 8e5                          ! atmos. diffusion coefficient [m^2/s]
  parameter( ce        = 2e-3  )                   ! laten heat transfer coefficient for ocean
  parameter( cq_latent = 2.257e6 )                 ! latent heat of condensation/evapoartion f water [J/kg]
  parameter( ci_latent = 3.335e5 )                 ! latent heat of condensation/fusion f ice [J/kg]
  parameter( cq_rain   = -0.1/24./3600. )          ! decrease in air water vapor due to rain [1/s]
  parameter( z_air     = 8400. )                   ! scaling height atmos. heat, CO2
  parameter( z_vapor   = 5000. )                   ! scaling height atmos. water vapor diffusion
  parameter( grav      = 9.81  )                   ! gravitational acceleration [m/s^2]
  real :: r_qviwv   = 2.6736e3                     ! regres. factor between viwv and q_air  [kg/m^3]

  ! physical paramter (rainfall)
  real :: c_q, c_rq, c_omega, c_omegastd

! parameter emisivity
  real, dimension(10) :: p_emi = (/9.0721, 106.7252, 61.5562, 0.0179, 0.0028,     &
&                                             0.0570, 0.3462, 2.3406, 0.7032, 1.0662/)

! declare climate fields
  real, dimension(xdim,ydim)          ::  z_topo, glacier,z_ocean
  real, dimension(xdim,ydim,nstep_yr) ::  Tclim, uclim, vclim, omegaclim, omegastdclim, wsclim
  real, dimension(xdim,ydim,nstep_yr) ::  qclim, mldclim, Toclim, cldclim
  real, dimension(xdim,ydim,nstep_yr) ::  TF_correct, qF_correct, ToF_correct, swetclim, dTrad
  real, dimension(ydim,nstep_yr)      ::  sw_solar, sw_solar_ctrl, sw_solar_scnr
  real, dimension(xdim,ydim)          ::  co2_part      = 1.0
  real, dimension(xdim,ydim)          ::  co2_part_scn  = 1.0

! declare anomaly fields for enso and climate change
  real, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_enso     = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_enso = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_enso    = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   Tclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   uclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   vclim_anom_cc       = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   omegaclim_anom_cc   = 0.
  real, dimension(xdim,ydim,nstep_yr) ::   wsclim_anom_cc      = 0.

! declare constant fields
  real, dimension(xdim,ydim)          ::  cap_surf
! ice sheet : capacity of snow  
  real, dimension(xdim,ydim)          ::  cap_ice
  integer jday, ityr

! Mike: declare some program constants
  real, dimension(xdim, ydim)         :: wz_air, wz_vapor
  real, dimension(xdim,ydim,nstep_yr) :: uclim_m, uclim_p
  real, dimension(xdim,ydim,nstep_yr) :: vclim_m, vclim_p

  real :: t0, t1, t2

  namelist / physics / log_exp, ct_sens, da_ice, a_no_ice, a_cloud, co_turb, kappa, 	&
&                      p_emi, Tl_ice1, Tl_ice2, To_ice1, To_ice2, r_qviwv,          	&
&		                   log_cloud_dmc, log_ocean_dmc, log_atmos_dmc, log_co2_dmc,        &
&                      log_hydro_dmc, log_qflux_dmc, 					                          &
&                      log_topo_drsp, log_cloud_drsp, log_humid_drsp, log_hydro_drsp,   &
&                      log_ocean_drsp, log_ice, log_hdif, log_hadv, log_vdif, log_vadv, &
& 		                 S0_var, dradius, log_rain, log_eva, log_conv, log_clim,          &
&                      log_tsurf_ext, log_hwind_ext, log_omega_ext, log_ice_sheet

end module mo_physics

!+++++++++++++++++++++++++++++++++++++++
module mo_diagnostics
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim

 ! declare diagnostic fields
  real, dimension(xdim,ydim)          :: Tsmn, Tamn, qmn, swmn, lwmn, qlatmn, qsensmn, &
&                                        ftmn, fqmn, icmn, Tomn, ice_Hmn, ice_Tsmn

! declare output fields
  real, dimension(xdim,ydim)          :: Tmm, Tamm, Tomm, qmm, icmm, prmm, evamm, qcrclmm, &
&                                        ice_Hmm, ice_Tsmm
  real, dimension(xdim,ydim,12)       :: Tmn_ctrl, Tamn_ctrl, Tomn_ctrl, ice_Hmn_ctrl, ice_Tsmn_ctrl, ice_mask_ctrl
  real, dimension(xdim,ydim,12)       :: qmn_ctrl, icmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl 

end module mo_diagnostics

!+++++++++++++++++++++++++++++++++++++++
subroutine greb_model
!+++++++++++++++++++++++++++++++++++++++
!   climate model main loop

  use mo_numerics
  use mo_physics
  use mo_diagnostics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts0, Ts1, Ta0, Ta1, To0, To1, q0, q1,       &
&                               ts_ini, ta_ini, q_ini, to_ini,              &
&                               iceH_ini,ice_Ts0, ice_Ts1, ice_H1, ice_H0

! open output files
  open(101,file='control.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(102,file='scenario.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(103,file='scenario.gmean.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal)
! ice sheet files
  open(201,file='control_ice_sheet.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(202,file='ice_sheet_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)


  dTrad = -0.16*Tclim -5. ! offset Tatmos-rad
! set ocean depth
  z_ocean=0
  do i=1,nstep_yr
     where(mldclim(:,:,i).gt.z_ocean) z_ocean = mldclim(:,:,i)
  end do
  z_ocean = 3.0*z_ocean

! decon mean state switch
  if (log_cloud_dmc ==  0) cldclim = 0.0
  if( log_hydro_dmc ==  0) qclim   = 0.0

! decon2xco2 switch
  if (log_topo_drsp  ==  0) where(z_topo > 1.) z_topo = 1.0     ! sens. exp. constant topo
  if (log_cloud_drsp ==  0) cldclim = 0.7                       ! sens. exp. constant cloud cover
  if (log_humid_drsp ==  0) qclim   = 0.0052                    ! sens. exp. constant water vapor
  if (log_ocean_drsp ==  0) mldclim = d_ocean                   ! sens. exp. no deep ocean

! heat capacity global [J/K/m^2]
  where (z_topo  > 0.) cap_surf = cap_land
  where (z_topo <= 0.) cap_surf = cap_ocean*mldclim(:,:,1)

! decon mean state switch
  if (log_ocean_dmc == 0) cap_surf = cap_land
! ice sheet : heat capacity of snow  
!  if (log_ice_sheet == 1) cap_ice  = cap_surf

! initialize fields
!  Ts_ini   = Tclim(:,:,nstep_yr)                          ! initial value temp. surf
  Ts_ini   = 273.15
!  Ta_ini   = Ts_ini                                       ! initial value atm. temp.
  Ta_ini   = Ts_ini
  To_ini   = Toclim(:,:,nstep_yr)                         ! initial value temp. surf
  q_ini    = qclim(:,:,nstep_yr)                          ! initial value atmos water vapor
  iceH_ini = 0.                                           ! initial ice sheet thickness

  CO2_ctrl = 340.0
! decon mean state switch
  if (log_co2_dmc == 0) CO2_ctrl = 0.

  if (log_exp .ge. 95 .and. log_exp .le. 100 )  CO2_ctrl = 280.  ! IPCC scenarios

  sw_solar = sw_solar_ctrl

! define some program constants
  wz_air   = exp(-z_topo/z_air)
  wz_vapor = exp(-z_topo/z_vapor)
  where (uclim(:,:,:) >= 0.0)
     uclim_m = uclim
     uclim_p = 0.0
  elsewhere
     uclim_m = 0.0
     uclim_p = uclim
  end where
  where (vclim(:,:,:) >= 0.0)
     vclim_m = vclim
     vclim_p = 0.0
  elsewhere
     vclim_m = 0.0
     vclim_p = vclim
  end where

  ! initialize the rainfall parameterisation
  select case( log_rain )
  case(-1) ! Original GREB
      c_q=1.; c_rq= 0.; c_omega=0.; c_omegastd=0.
    case(1) ! Adding relative humidity (rq)
      c_q=-1.391649; c_rq=3.018774; c_omega= 0.; c_omegastd=0.
    case(2) ! Adding omega
      c_q=0.862162; c_rq=0.; c_omega=-29.02096; c_omegastd=0.
    case(3) ! Adding rq and omega
      c_q=-0.2685845; c_rq=1.4591853; c_omega=-26.9858807; c_omegastd=0.
    case(0) ! Best GREB
        c_q=-1.88; c_rq=2.25; c_omega=-17.69; c_omegastd=59.07 ! Rainfall parameters (ERA-Interim)
      if (log_clim == 1) then
        c_q=-1.27; c_rq=1.99; c_omega=-16.54; c_omegastd=21.15 ! Rainfall parameters (NCEP)
      end if
  end select


! control run
  print*,'% CONTROL RUN CO2=',CO2_ctrl,'  time=', time_ctrl,'yr'
  Ts1 = Ts_ini; Ta1 = Ta_ini; To1 = To_ini; q1 = q_ini; ice_Ts1 = Ts1; ice_H1 = iceH_ini                  ! initialize fields / ice sheet 
  year=1970; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  do it=1, time_ctrl*nstep_yr                                             ! main time loop
    call time_loop(it, isrec, year, CO2_ctrl, irec, mon, 101, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, ice_Ts0, ice_H0, ice_Ts1, ice_H1 )
    Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0; ice_Ts1=ice_Ts0; ice_H1=ice_H0
    if (log_exp .eq. 1 .and. mod(it,nstep_yr) .eq. 0) year=year+1
  end do



! scenario run
if ( log_exp .ne. 1 .or. time_scnr .ne. 0 ) then
  if( log_exp .eq. 30  ) sw_solar = sw_solar_scnr ! paleo 231 kyr bp
  if( log_exp .eq. 31  ) sw_solar = sw_solar_scnr ! paleo 231 kyr bp
  if( log_exp .eq. 35  ) sw_solar = sw_solar_scnr ! change obliquity
  if( log_exp .eq. 36  ) sw_solar = sw_solar_scnr ! change eccentricity
  if( log_ep .eq. 37  ) then ! change solar constant as function of radius
     radius = 1+0.01*(dradius)
     print*,'Solar radius [AU] = ', radius
     rS0 = (1/radius)**2
     sw_solar = rS0*sw_solar
  end if
  if ( log_exp .eq. 230 ) then ! ice sheet setting ; change boundary conditions for Climate Change forcing
     Tclim      = Tclim + Tclim_anom_cc
     uclim      = uclim + uclim_anom_cc
     vclim      = vclim + vclim_anom_cc
     omegaclim  = omegaclim + omegaclim_anom_cc
     wsclim     = wsclim + wsclim_anom_cc
  end if
  if ( log_exp .eq. 240 .or. log_exp .eq. 241 ) then ! change boundary conditions for ENSO forcing
     Tclim      = Tclim + Tclim_anom_enso
     uclim      = uclim + uclim_anom_enso
     vclim      = vclim + vclim_anom_enso
     omegaclim  = omegaclim + omegaclim_anom_enso
     wsclim     = wsclim + wsclim_anom_enso
  end if

  print*,'% SCENARIO EXP: ',log_exp,'  time=', time_scnr,'yr'
  Ts1 = Ts_ini; Ta1 = Ta_ini; q1 = q_ini; To1 = To_ini; ice_H1 = ice_ini; ice_Ts1 = Ts1                     ! initialize field / ice sheet
  year=1950.; CO2=340.0; mon=1; irec=0; Tmm=0.; Tamm=0.; qmm=0.; apmm=0.;
  if (log_exp .ge. 35 .and. log_exp .le. 37) year=1.

  do it=1, time_scnr*nstep_yr                                              ! main time loop
     call forcing(it, year, CO2, Ts1)
     call time_loop(it,isrec, year, CO2, irec, mon, 102, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, ice_Ts0, ice_H0, ice_Ts1, ice_H1 )
     Ts1=Ts0; Ta1=Ta0; q1=q0; To1=To0; ice_Ts1=Ts_ini; ice_H1=ice_H0
     if (mod(it,nstep_yr) == 0) year=year+1
  end do

end if !( log_exp .ne. 1 )

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine time_loop(it, isrec, year, CO2, irec, mon, ionum, Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, ice_Ts0, ice_H0, ice_Ts1, ice_H1)
!+++++++++++++++++++++++++++++++++++++++
! main time loop

  use mo_numerics
  use mo_physics

  real, dimension(xdim,ydim):: Ts1, Ta1, q1, To1, Ts0,Ta0, q0, To0, sw,       &
&                              ice_cover, Q_sens, Q_lat, Q_lat_air, dq_eva,   &
&                              dq_rain, dTa_crcl, dq_crcl, dq, dT_ocean, dTo, &
&                              LW_surf, LWair_down, LWair_up, em, Fn_surf,    &
&                              ice_dH, ice_H1, ice_H0,  ice_Ts0, ice_Ts1

  jday = mod((it-1)/ndt_days,ndays_yr)+1  ! current calendar day in year
  ityr = mod((it-1),nstep_yr)+1           ! time step in year

  call tendencies(CO2, Ts1, Ta1, To1, q1, ice_cover, SW, LW_surf, Q_lat,      &
&                    Q_sens, Q_lat_air, dq_eva, dq_rain, dq_crcl,             &
&                    dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em,       &
&                    ice_dH)

  Tmin_limit = 40 ! no very low Tsurf/Tatmoss;  numerical stability


  ! surface temperature
  Fn_surf = 3.9 
  Ts0  = Ts1 + dt*Fn_surf / cap_surf
  where(Ts0 .le. Tmin_limit )     Ts0 = Tmin_limit ! no very low Tsurf;  numerical stability
  Ta0  = 273.15

  ! sea ice heat capacity
  ! ice sheet : spread and ablation
  if(log_ice_sheet == 1) call ice_sheet(it, ionum, irec, mon, ice_Ts0, ice_H0, ice_Ts1, ice_H1, ice_dH, Fn_surf, dT_ocean, z_topo)
  ! diagnostics: annual means plots
  where(ice_H0>0.) Ts0 = ice_Ts0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for test 
end subroutine time_loop

!+++++++++++++++++++++++++++++++++++++++
subroutine tendencies(CO2, Ts1, Ta1, To1, q1, ice_cover, SW, LW_surf, Q_lat, Q_sens, Q_lat_air,  &
&                     dq_eva, dq_rain, dq_crcl, dTa_crcl, dT_ocean, dTo, LWair_down, LWair_up, em, &
&                     ice_dH)  ! ice sheet 
!+++++++++++++++++++++++++++++++++++++++

  use mo_numerics
  use mo_physics

! declare temporary fields
  real, dimension(xdim,ydim) :: Ts1, Ta1, To1, q1, ice_cover, sw, LWair_up,    &
&                               LWair_down, em, Q_sens, Q_lat, Q_lat_air,      &
&                               dq_eva, dq_rain, dTa_crcl, dq_crcl, LW_surf,   &
&                               dT_ocean, dTo, ice_dH

!$omp parallel sections
!$omp section
    ! SW radiation model
    call SWradiation(Ts1, sw, ice_cover)
!$omp section
    ! LW radiation model
    call LWradiation(Ts1, Ta1, q1, CO2, LW_surf, LWair_up, LWair_down, em)
    ! sensible heat flux
    Q_sens = ct_sens*(Ta1-Ts1)
!$omp section
! decon mean state switch
    if (log_atmos_dmc == 0) Q_sens = 0.

    ! hydro. model
    ! ice sheet : ice accumulation, thickness increase by snowfall
    if (log_ice_sheet == 1) call ice_accumulation(ice_dH, Ts1, Ta1, dq_rain, wz_vapor) 
    ! atmos. circulation

end subroutine tendencies


!+++++++++++++++++++++++++++++++++++++++
subroutine SWradiation(Tsurf, sw, ice_cover)
!+++++++++++++++++++++++++++++++++++++++
!    SW radiation model

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: ityr, sw_solar,da_ice, a_no_ice, a_cloud, z_topo  &
&                         , Tl_ice1, Tl_ice2, To_ice1, To_ice2, glacier       &
&                         , cldclim, log_exp, log_atmos_dmc, log_ice, S0_var

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, sw, albedo, ice_cover, a_surf, a_atmos

! atmos albedo
  a_atmos=cldclim(:,:,ityr)*a_cloud

! isnow/ice diagnostic only
  ice_cover=0.0
  where(z_topo >= 0. .and. Tsurf <= Tl_ice1) ice_cover = 1.0
  where(z_topo >= 0. .and. Tsurf >  Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&       ice_cover = (1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))
  where(z_topo < 0.  .and. Tsurf <= To_ice1) ice_cover = 1.0
  where(z_topo < 0.  .and. Tsurf >  To_ice1 .and. Tsurf < To_ice2 ) &
&       ice_cover = (1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

! surface albedo
! Land:  ice -> albedo linear function of T_surf
   where(z_topo >= 0. .and. Tsurf <= Tl_ice1) a_surf = a_no_ice+da_ice   ! ice
   where(z_topo >= 0. .and. Tsurf >= Tl_ice2) a_surf = a_no_ice          ! no ice
   where(z_topo >= 0. .and. Tsurf > Tl_ice1 .and. Tsurf < Tl_ice2 ) &
&       a_surf = a_no_ice +da_ice*(1-(Tsurf-Tl_ice1)/(Tl_ice2-Tl_ice1))
! Ocean: ice -> albedo/heat capacity linear function of T_surf
  where(z_topo < 0. .and. Tsurf <= To_ice1) a_surf = a_no_ice+da_ice      ! ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) a_surf = a_no_ice             ! no ice
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       a_surf = a_no_ice+da_ice*(1-(Tsurf-To_ice1)/(To_ice2-To_ice1))

! glacier -> no albedo changes
  where(glacier > 0.5) a_surf = a_no_ice+da_ice

! dmc & decon2xco2 switch
  if (log_ice  ==    0) a_surf = a_no_ice

! SW flux
  albedo=a_surf+a_atmos-a_surf*a_atmos
  forall (i=1:xdim)
     sw(i,:)=0.01*S0_var*SW_solar(:,ityr)*(1-albedo(i,:))
  end forall

end subroutine SWradiation


!+++++++++++++++++++++++++++++++++++++++
subroutine LWradiation(Tsurf, Tair, q, CO2, LWsurf, LWair_up, LWair_down, em)
!+++++++++++++++++++++++++++++++++++++++
! new approach with LW atmos

  USE mo_numerics,    ONLY: xdim, ydim
  USE mo_physics,     ONLY: sig, eps, qclim, cldclim, z_topo, jday, ityr,         &
&                           r_qviwv, z_air, z_vapor, dTrad, p_emi, log_exp,       &
&                           log_atmos_dmc, co2_part

! declare temporary fields
  real, dimension(xdim,ydim)  :: Tsurf, Tair, q, LWsurf, LWair, e_co2, e_cloud,   &
&                                LWair_up, LWair_down, e_vapor, em


  e_co2   = exp(-z_topo/z_air)*co2_part*CO2   ! CO2
  e_vapor = exp(-z_topo/z_air)*r_qviwv*q      ! water vapor
  e_cloud = cldclim(:,:,ityr)                 ! clouds

! total
  em      = p_emi(4)*log( p_emi(1)*e_co2 +p_emi(2)*e_vapor +p_emi(3) ) +p_emi(7)   &
&          +p_emi(5)*log( p_emi(1)*e_co2   +p_emi(3) )                             &
&          +p_emi(6)*log( p_emi(2)*e_vapor +p_emi(3) )
  em      = (p_emi(8)-e_cloud)/p_emi(9)*(em-p_emi(10))+p_emi(10)

  LWsurf      = -sig*Tsurf**4
  LWair_down  = -em*sig*(Tair+dTrad(:,:,ityr))**4
  LWair_up    = LWair_down

! decon mean state switch
  if( log_atmos_dmc == 0) LWair_down = 0

end subroutine LWradiation





subroutine forcing(it, year, CO2, Tsurf)
!+++++++++++++++++++++++++++++++++++++++

  USE mo_numerics,    ONLY: xdim, ydim, ndays_yr, ndt_days, nstep_yr
  USE mo_physics,     ONLY: log_exp, sw_solar, sw_solar_ctrl, sw_solar_scnr,     &
&                           co2_part, co2_part_scn, z_topo, ityr, Tclim
  USE mo_diagnostics,  ONLY: icmn_ctrl

  ! input fields
  real, dimension(xdim,ydim)  :: Tsurf

  ! declare variables
  real, dimension(xdim,ydim) 	:: icmn_ctrl1

  ! calculate annual mean ice cover (for log_exp 44 & 45)
  icmn_ctrl1 = 0.0
  do i=1,12
  icmn_ctrl1 = icmn_ctrl1 + icmn_ctrl(:,:,i)
  end do
  icmn_ctrl1 = icmn_ctrl1/12.


  if( log_exp .eq. 10 ) CO2 = 2*340.

  if( log_exp .eq. 20 ) CO2 = 2*340.
  if( log_exp .eq. 21 ) CO2 = 4*340.
  if( log_exp .eq. 22 ) CO2 = 10*340.
  if( log_exp .eq. 23 ) CO2 = 0.5*340.
  if( log_exp .eq. 24 ) CO2 = 0*340.

  if( log_exp .eq. 25 ) CO2 = 340.+170.+170.*cos(2*3.14* (year-13.)/30.) ! sinus wave
  if( log_exp .eq. 26 ) CO2 = 2*340. ! step function
  if( log_exp .eq. 26 .and. year >= 1980 ) CO2 = 340. ! step function
  if( log_exp .eq. 27  ) CO2      = 340.
  if( log_exp .eq. 27  ) dS0      = 27.0 ! [W/m2]
  if( log_exp .eq. 27  ) sw_solar = (1365+dS0)/1365*sw_solar_ctrl ! change S0 by dS0
  if( log_exp .eq. 28  ) CO2      = 340.
  if( log_exp .eq. 28  ) dS0      = 1.0 ! amplitude solar cycle [W/m2]
  if( log_exp .eq. 28  ) sw_solar = (1365+dS0*sin(2*3.14* (year)/11.))/1365*sw_solar_ctrl ! 11yrs cycle

  if( log_exp .eq. 30  ) CO2      = 200. ! paleo 231 kyr BP
  if( log_exp .eq. 30  ) sw_solar = sw_solar_scnr !
  if( log_exp .eq. 31  ) CO2      = 340.
  if( log_exp .eq. 31  ) sw_solar = sw_solar_scnr
  if( log_exp .eq. 32  ) CO2      = 200.

! partial scenarios
! 2xCO2 NH
  if( log_exp .eq. 40 )  CO2      = 2*340.
  if( log_exp .eq. 40 )  co2_part(:,1:24)     = 0.5
! 2xCO2 SH
  if( log_exp .eq. 41 )  CO2      = 2*340.
  if( log_exp .eq. 41 )  co2_part(:,25:48)    = 0.5
! 2xCO2 tropics
  if( log_exp .eq. 42 )  CO2      = 2*340.
  if( log_exp .eq. 42 )  co2_part(:,1:15)     = 0.5
  if( log_exp .eq. 42 )  co2_part(:,33:48)    = 0.5
  if( log_exp .eq. 42 )  co2_part(4:4:96,33)  = 1.0
  if( log_exp .eq. 42 )  co2_part(4:4:96,15)  = 1.0
! 2xCO2 extra tropics
  if( log_exp .eq. 43 )  CO2      = 2*340.
  if( log_exp .eq. 43 )  co2_part(:,16:32)    = 0.5
  if( log_exp .eq. 43 )  co2_part(4:4:96,32)  = 1.0
  if( log_exp .eq. 43 )  co2_part(4:4:96,16)  = 1.0
! 2xCO2 ocean
  if( log_exp .eq. 44 )  CO2      = 2*340.
  if( log_exp .eq. 44 )  where (z_topo  > 0.) co2_part = 0.5
  if( log_exp .eq. 44 )  where (icmn_ctrl1  >= 0.5 ) co2_part = 0.5
! 2xCO2 land/ice
  if( log_exp .eq. 45 )  CO2      = 2*340.
  if( log_exp .eq. 45 )  where (z_topo  <= 0.) co2_part = 0.5
  if( log_exp .eq. 45 )  where (icmn_ctrl1  >= 0.5 ) co2_part = 1.0
! 2xCO2 winter
  if( log_exp .eq. 46 )  CO2      = 340.
  if( log_exp .eq. 46 .and. mod(it,nstep_yr) .le. 181)  CO2 = 2*340.
  if( log_exp .eq. 46 .and. mod(it,nstep_yr) .ge. 547)  CO2 = 2*340.
! 2xCO2 summer
  if( log_exp .eq. 47 )  CO2      = 2*340.
  if( log_exp .eq. 47 .and. mod(it,nstep_yr) .le. 181)  CO2 = 340.
  if( log_exp .eq. 47 .and. mod(it,nstep_yr) .ge. 547)  CO2 = 340.

! IPCC A1B scenario
  if( log_exp .eq. 95 ) then
     CO2_1950=310.;  CO2_2000=370.;  CO2_2050=520.
     if (year <= 2000.) CO2=CO2_1950 + 60./50.*(year-1950.)
     if (year  > 2000. .and. year <= 2050.) CO2=CO2_2000 + 150./50.*(year-2000.)
     if (year  > 2050. .and. year <= 2100.) CO2=CO2_2050 + 180./50.*(year-2050.)
  end if

! IPCC RCP scenarios
  if( log_exp .ge. 96 .and. log_exp .le. 99 ) then ! IPCC RCP scenarios
      if(mod(it,nstep_yr) .eq. 1) read (26,*) t, CO2
  end if

! own CO2 scenario
  if( log_exp .eq. 100 ) then
      if(mod(it,nstep_yr) .eq. 1) read (26,*) t, CO2
  end if

! Forced Climate Change run
  if( log_exp .eq. 230 ) Tsurf = Tclim(:,:,ityr) ! Keep temp on external boundary condition

! Forced ENSO run
  if( log_exp .eq. 240 .or. log_exp .eq. 241 ) Tsurf = Tclim(:,:,ityr)  ! Keep temp on external boundary condition


end subroutine forcing


!+++++++++++++++++++++++++++++++++++++++
subroutine ice_accumulation(ice_dH, Ts1, Ta1, dq_rain, wz_vapor) 
!+++++++++++++++++++++++++++++++++++++++
! ice sheet : accumulation process
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics,  ONLY: r_qviwv, Tl_ice2, ice_svlm
  implicit none
  real, dimension(xdim,ydim) :: Ts1                     ! surface temperature [K] 
  real, dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real, dimension(xdim,ydim) :: ice_dH                  ! ice thickness tendency [kg/m2/s]
  real, dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]
  real, dimension(xdim,ydim) :: wz_vapor                ! surface pressure change coefficient [1]
  ice_dH = 0.
  dq_rain = -1./86400.
  ! ice sheet accumulation from snowfall, kg/m2/s   
!  where((Ts1 <= Tl_ice2) .and. (Ta1 <= Tl_ice2)) ice_dH = - ice_svlm*dq_rain*r_qviwv*wz_vapor; 
  where((Ts1 <= Tl_ice2) .and. (Ta1 <= Tl_ice2)) ice_dH = - ice_svlm*dq_rain; 
end subroutine ice_accumulation

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_sheet(it, ionum, irec, mon, ice_Ts0, ice_H0, ice_Ts1, ice_H1, ice_dH, Fn_surf, dT_ocean, z_topo)
!+++++++++++++++++++++++++++++++++++++++
  USE mo_numerics, ONLY: xdim, ydim, dt
  USE mo_physics,  ONLY: cap_ice, rho_ice, cp_ice, d_ice_max, cap_surf
  implicit none
  real, dimension(xdim,ydim) :: ice_Ts0                 ! ice surface temperature (forward) [K]
  real, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature (current) [K]
  real, dimension(xdim,ydim) :: Fn_surf                 ! surface net heat flux [J/m2]
  real, dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real, dimension(xdim,ydim) :: U_gtmlt                 ! ocean temperature change [K]

  real, dimension(xdim,ydim) :: ice_H0                  ! ice thickness (forward) [m]
  real, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (current) [m]
  real, dimension(xdim,ydim) :: ice_dH                  ! ice thickness tendency [m]
  real, dimension(xdim,ydim) :: melt                    ! ice melting accumulation [m]
  real, dimension(xdim,ydim) :: z_topo                  ! surface topography [m]
  real, dimension(xdim,ydim) :: ice_cover               ! ice cover type [-1~0:partial ice sheet;0 land;0~1 partial sea ice]

  real                       :: Tmin_limit              ! no very low Tsurf/Tatmoss;  numerical stability

  integer                    :: ionum                   ! file unit of GREB
  integer                    :: ice_iounit              ! file unit of ice sheet data
  integer                    :: it, irec, mon           ! work variable for count

  Tmin_limit = 40
  cap_ice    = cap_surf
  ! snow heat capacity
  where((ice_H1 <  d_ice_max).and.(ice_H1 > 0.)) cap_ice = ice_H1*cp_ice*rho_ice ! heat capacity of snow [J/K/m^2]
  where((ice_H1 >= d_ice_max).and.(ice_H1 > 0.)) cap_ice = d_ice_max*cp_ice*rho_ice ! heat capacity limitation of snow [J/K/m^2]
  
  ! snowfall accumulation
  where(z_topo >= 0) ice_H0 = ice_H1 + dt*ice_dH 
  
  ! ice surface temperature
  ice_Ts0  = ice_Ts1 + dT_ocean + dt*Fn_surf / cap_ice 
  where(ice_Ts0 .le. Tmin_limit ) ice_Ts0 = Tmin_limit ! no very low Tsurf;  numerical stability
  ! fusion correction
  call ice_fusion( ice_Ts0, ice_H0, melt, U_gtmlt)
  ! thickness after fusion
  
 !  where(ice_H0 > 0.) cap_surf = cap_land                       ! ice sheet heat capacity
  ! ice sheet output
  ice_iounit = 100 + ionum
  call ice_output(it, ice_iounit, irec, mon, ice_Ts0, ice_H0, U_gtmlt)
  
end subroutine ice_sheet

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_fusion( ice_Ts0, ice_H0, melt,U_gtmlt) 
!+++++++++++++++++++++++++++++++++++++++
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics,  ONLY: ci_latent, Tl_ice2, cap_ice, cap_land, rho_ice
  implicit none
  real, dimension(xdim,ydim) :: ice_Ts0                 ! ice surface temperature [K]
  real, dimension(xdim,ydim) :: ice_H0                  ! ice thickness [m]
  real, dimension(xdim,ydim) :: melt                    ! ice melting accumulation [m]
  real, dimension(xdim,ydim) :: U_gtmlt                 ! internal energy greater than melting point [J/m2]
  real, dimension(xdim,ydim) :: Lm_max                  ! potential energy of snow fusion [J/m2] 

  melt    = 0.
  
  where( ice_H0 > 0.)
      U_gtmlt = (ice_Ts0 - Tl_ice2) * cap_ice
      Lm_max  = ci_latent * rho_ice * ice_H0
  end where
  ! surface snow totally melts away
  where((ice_Ts0 >= Tl_ice2) .and. (U_gtmlt >  Lm_max) .and. (ice_H0 >0.))
        melt    = ice_H0
        ice_Ts0 = Tl_ice2 + (U_gtmlt - Lm_max) / cap_land
  end where
  ! surface snow partially melts
  where((ice_Ts0 >= Tl_ice2) .and. (U_gtmlt <= Lm_max) .and. (ice_H0 >0.))
        melt    = U_gtmlt / (rho_ice*ci_latent)
        ice_Ts0 = Tl_ice2
  end where
  ice_H0 = ice_H0 - melt

end subroutine ice_fusion

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_output(it, ice_iounit, irec, mon, ice_Ts0, ice_H0, U_gtmlt )
!+++++++++++++++++++++++++++++++++++++++
! ice sheet : output file
  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, nstep_yr, time_scnr &
&                          , time_ctrl, ireal, dt
  USE mo_physics,      ONLY: jday, log_exp, r_qviwv, wz_vapor, cap_ice
  use mo_diagnostics,  ONLY: ice_Tsmm, ice_Tsmn_ctrl, ice_Hmn_ctrl, ice_mask_ctrl
  implicit none
  real, external              :: gmean
  real, dimension(xdim,ydim)  :: ice_H0, ice_Ts0, ice_mask, melt, ice_melt, U_gtmlt 
  integer, parameter          :: nvar = 3              ! number of output variable
  integer                     :: it,irec,mon,iyrec     ! work variable for count, controled by subroutine output
  integer                     :: ndm                   ! total time for mean calculation
  integer                     :: ice_iounit            ! written file uint

  ! diagnostics: monthly means
  ice_mask = 0.; ice_melt = 0.
  where(ice_H0 > 0.) ice_mask = 1.
  ice_Tsmm = ice_Tsmm+ice_Ts0
  ice_melt =  U_gtmlt
  
! control output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. ice_iounit == 201 ) then
     ndm=jday_mon(mon)*ndt_days
     if (it/float(ndt_days)  > 365*(time_ctrl-1)) then
         if (log_exp .eq. 1 .or. log_exp .eq. 310 ) then
         write(ice_iounit,rec=nvar*irec+1)  ice_Tsmm/ndm
         write(ice_iounit,rec=nvar*irec+2)  ice_H0
         write(ice_iounit,rec=nvar*irec+3)  ice_mask
         irec=irec+1;
         else
         ice_Tsmn_ctrl(:,:,mon)  = ice_Tsmm/ndm
         ice_Hmn_ctrl(:,:,mon)   = ice_H0
         ice_mask_ctrl(:,:,mon)  = ice_melt
         end if
     end if
     ice_Tsmm=0.
     mon=mon+1; if (mon==13) mon=1
  end if

! scenario output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. ice_iounit == 202 ) then

     ndm=jday_mon(mon)*ndt_days
     write(ice_iounit,rec=nvar*irec+1)  ice_Tsmm/ndm
     write(ice_iounit,rec=nvar*irec+2)  ice_H0
     write(ice_iounit,rec=nvar*irec+3)  ice_melt
     irec = irec + 1
     iyrec=floor(float((irec-1)/12))
     ice_Tsmm=0.; ice_melt=0.
     mon=mon+1; if (mon==13) mon=1
  end if

end subroutine ice_output
!TB
!+++++++++++++++++++++++++++++++++++++++
function gmean(data)
!+++++++++++++++++++++++++++++++++++++++

use mo_numerics,		ONLY: xdim, ydim, dlat

! declare variables
real, dimension(xdim,ydim) 	:: data, w
real, dimension(ydim)		:: lat

do i=1,ydim
lat(i) = -90-(dlat*0.5)+i*dlat
end do
do i=1,xdim
w(i,:) = cos(2.*3.14159*lat/360.)
end do

gmean = sum(data*w)/sum(w)

end function
