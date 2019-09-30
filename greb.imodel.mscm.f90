
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
  parameter( rho_ice   = 910.)           ! density of incompressible polycrystalline ice [kg/m^3]
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
end module

! ice sheet model for GREB
! start 26 Sep 2019
program ice_sheet_model
  implicit none
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndt_days  = 24*3600           ! number of timesteps per day
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: ndt_yr    = ndays_yr*ndt_days        ! number of timesteps per day
  real               :: dt        = ndt_yr/2           ! time step [s]
  integer            :: ireal     = 4
  integer            :: j, i, irec, nstep_end
  real,parameter     :: pi        = 3.1416 
  real,parameter     :: dx        = 2*pi/(xdim)
  real, dimension(xdim,ydim) :: iceH_ini, ice_H1, ice_H0, ice_vx, diceH_crcl
  real, dimension(xdim,ydim) :: term_l1, term_l2, term_l3, term_r

  nstep_end = 96*10
 
  open(301,file='ice_scheme_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  ice_vx = dx/ndt_yr
  do j = 1,xdim
      do i = 1,ydim
          iceH_ini(j,i) = exp(-dx*(j-48)**2)
      !iceH_ini(j) = sin(j*2*pi/xdim)
      end do
  end do
  ice_H1  = iceH_ini

  do irec = 1,nstep_end
      do i = 1,ydim
      j = 1
      term_l2 = 1
      diceH_crcl(j,i) = ice_vx(j+1,i)*ice_H1(j+1,i) - ice_vx(xdim,i)*ice_H1(xdim,i)
      term_l1(j,i)    = dt/(4*dx)*ice_vx(xdim,i)
      term_l3(j,i)    = -dt/(4*dx)*ice_vx(2,i)
      j = xdim
      diceH_crcl(j,i) = ice_vx(1,i)*ice_H1(1,i) - ice_vx(j-1,i)*ice_H1(j-1,i)
      term_l1(j,i)    = dt/(4*dx)*ice_vx(xdim-1,i)
      term_l3(j,i)    = -dt/(4*dx)*ice_vx(1,i)
      do j = 2, xdim-1
          diceH_crcl(j,i) = ice_vx(j+1,i)*ice_H1(j+1,i) - ice_vx(j-1,i)*ice_H1(j-1,i) 
          term_l1(j,i)    = dt/(4*dx)*ice_vx(j-1,i)
          term_l3(j,i)    = -dt/(4*dx)*ice_vx(j+1,i)
      end do
      end do
      diceH_crcl = diceH_crcl/(4*dx)
!      dt = -sum(diceH_crcl*ice_H2)*dx**2/(sum(diceH_crcl*diceH_crcl)*dx**2)*2.

      term_r  = ice_H1 + dt*diceH_crcl
      !call Guass_Seidel(term_l1,term_l2,term_l3,ice_H0,term_r,ice_H1)
      call advection_finite_volume_pwl_slimit(ice_H1,diceH_crcl,dt,ice_vx)
      ice_H0 = ice_H1+diceH_crcl

  write(301,rec=irec) ice_H0
  
  ice_H1 = ice_H0

  end do
end

subroutine Guass_Seidel(coff_0,coff_1,coff_2,x0,f,x1)
  implicit none
  integer, parameter    :: xdim = 96, ydim = 48          ! field dimensions
  integer               :: iter, k,j, n, itmax
  real                  :: omega, epsiron
  real,dimension(xdim,ydim)  :: coff_0,coff_1,coff_2,x0,x1,f   ! input coefficient
  n       = xdim
  itmax   = 1000
  epsiron = 1e-10
  omega   = 1.
  do iter = 1,itmax
  do j = 1,ydim
      x0(1,j) = (1-omega)*x1(1,j)+omega*(f(1,j)-coff_0(1,j)*x1(n,j)-coff_2(1,j)*x1(2,j))/coff_1(1,j)
      do k=2,n-1
          x0(k,j) = (1-omega)*x1(k,j)+omega*(f(k,j)-coff_0(k,j)*x0(k-1,j)-coff_2(k,j)*x0(k+1,j))/coff_1(k,j)
      end do
      x0(n,j) = (1-omega)*x1(n,j)+omega*(f(n,j)-coff_0(n,j)*x1(n-1,j)-coff_2(n,j)*x0(1,j))/coff_1(n,j)

      if(sum((x0-x1)**2)<epsiron) then
          exit
      end if
  end do

      x1 = x0

  end do
end subroutine


subroutine Tridiagonal_matrix_algorithm(coff_0,coff_1,coff_2,x,f)
! According to "The Theory of Splines and Their Applications" (1976) P15 by Ahlberg et al.
! solving tridiagonal matrix equation set with periodic condition
  implicit none
  integer, parameter    :: xdim = 96, ydim = 48          ! field dimensions
  integer               :: k, n
  real,dimension(xdim)  :: coff_0, coff_1, coff_2, x, f  ! input coefficient
  real,dimension(xdim)  :: p, q, u, s                    ! variables for solving matrix
  real,dimension(xdim)  :: t, v                          ! variables for substitute
  n = xdim
  p(1)=coff_1(1);q(1)=-coff_2(1)/p(1);u(1)=f(1)/p(1);s(1)=-coff_0(1)/p(1);   ! initialization
  ! x(k) = q(k)*x(k+1)+s(k)*x(n)+u(k)
  do k= 2,n-1
      p(k) = coff_0(k)*q(k-1)+coff_1(k)    ! coefficient for diagonal element
      q(k) = -coff_2(k)/p(k)               ! coefficient for diagonal element next step
      u(k) = (f(k)+coff_0(k)*u(k-1))/p(k)  ! coefficient for offset
      s(k) = -coff_0(k)*s(k-1)/p(k)        ! coefficient for xn 
  end do
  ! x(k) = t(k)*x(n)+v(k)
  t(n)=1;v(n)=0 ! initialization
  do k= n-1,1,-1
      t(k) = q(k)*t(k+1)+s(k)              ! coefficient for xn
      v(k) = q(k)*v(k+1)+u(k)              ! coefficient for offset
  end do
  ! xn
  x(n) = (f(n)-v(n-1)*coff_0(n)-v(1)*coff_2(n))/(coff_1(n)+coff_0(n)*t(n-1)+coff_2(n)*t(1))
  ! xk
  do k= 1,n-1
      x(k) = t(k)*x(n) + v(k)
  end do
  
end subroutine


subroutine advection_finite_volume_pwl_slimit(T1, dX_advec,dt, wz)
!+++++++++++++++++++++++++++++++++++++++
    USE mo_numerics, ONLY: xdim, ydim, dlon, dlat, dt_crcl
    USE mo_physics,  ONLY: pi, z_topo, uclim, vclim, ityr, z_vapor, log_exp
    IMPLICIT NONE
    real, dimension(xdim,ydim), intent(in)  :: T1, wz
    real, dimension(xdim,ydim), intent(out) :: dX_advec
    real, dimension(xdim,ydim)              :: Twz,adv_flux, adv_no
    integer :: i
    integer, dimension(ydim):: ilat = (/(i,i=1,ydim)/)
    real, dimension(ydim)   :: lat, dx
    real                    :: sa, sb, sc, sd
    real                    :: uA, vB, uC, vD
    real                    :: fA, fB, fC, fD
    real    :: deg, dy, dt
    integer :: j, k, km1, kp1, kp2, kp3, jm1, jp1, jp2, jp3

    deg      = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
    dy       = dlat*deg
    lat      = dlat*ilat-dlat/2.-90.
    !dx       = dlon*deg*cos(2.*pi/360.*lat)
    dx       = 2*pi/xdim
    adv_no    = 0.
    adv_flux = wz
    dX_advec = 0.

    !< Get advection variable by multiplying variable with wz
    Twz = T1

    !< Whole domain except boundaries and corners
    do k=2, ydim-3
        do j=2, xdim-3
                km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=j-1; jp1=j+1; jp2=j+2; jp3=j+3

                !< Slopes
                sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
                sb = (Twz(j,kp1)-Twz(j,k)) / dy
                sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
                sd = (Twz(j,k)-Twz(j,km1)) / dy

                !< Mean winds on boundaries of cell
                uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
                vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
                uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
                vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

                !< Flux over cell coundaries
                !< fA
                if (uA > 0.) then
                    sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
                    fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
                else if (uA < 0.) then
                    sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
                    fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
                else
                    fA = 0.
                end if

                !< fB
                if (vB > 0.) then
                    sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
                    fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
                else if (vB < 0.) then
                    sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
                    fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
                else
                    fB = 0.
                end if

                !< fC
                if (uC > 0.) then
                    sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
                    fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
                else if (uC < 0.) then
                    sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
                    fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
                else
                    fC = 0.
                end if

                !< fD
                if (vD > 0.) then
                    sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
                    fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
                else if (vD < 0.) then
                    sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
                    fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
                else
                    fD = 0.
                end if

                dX_advec(j,k) = dt * ( fA - fB - fC + fD )

           end do
    end do

    !< Boundaries
    !< Left boundary
    j=1
    do k=2, ydim-3
        km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=xdim; jp1=j+1; jp2=j+2; jp3=j+3

        !< Slopes
        sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
        sb = (Twz(j,kp1)-Twz(j,k)) / dy
        sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
        sd = (Twz(j,k)-Twz(j,km1)) / dy

        !< Mean winds on boundaries of cell
        uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
        vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
        uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
        vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

        !< Flux over cell coundaries
        !< fA
        if (uA > 0.) then
            sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
            fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else if (uA < 0.) then
            sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else
            fA = 0.
        end if

        !< fB
        if (vB > 0.) then
            sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else if (vB < 0.) then
            sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
            fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else
            fB = 0.
        end if

        !< fC
        if (uC > 0.) then
            sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else if (uC < 0.) then
            sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
            fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else
            fC = 0.
        end if

        !< fD
        if (vD > 0.) then
            sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
            fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else if (vD < 0.) then
            sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else
            fD = 0.
        end if

        dX_advec(j,k) = dt * ( fA - fB - fC + fD )

    end do

    !< Right boundary
    j=xdim
    do k=2, ydim-3
        km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=j-1; jp1=1; jp2=2; jp3=3

        !< Slopes
        sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
        sb = (Twz(j,kp1)-Twz(j,k)) / dy
        sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
        sd = (Twz(j,k)-Twz(j,km1)) / dy

        !< Mean winds on boundaries of cell
        uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
        vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
        uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
        vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

        !< Flux over cell coundaries
        !< fA
        if (uA > 0.) then
            sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
            fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else if (uA < 0.) then
            sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else
        fA = 0.
        end if

        !< fB
        if (vB > 0.) then
            sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else if (vB < 0.) then
            sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
            fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else
        fB = 0.
        end if

        !< fC
        if (uC > 0.) then
            sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else if (uC < 0.) then
            sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
            fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else
            fC = 0.
        end if

        !< fD
        if (vD > 0.) then
            sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
            fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else if (vD < 0.) then
            sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else
            fD = 0.
        end if

        dX_advec(j,k) = dt * ( fA - fB - fC + fD )

    end do

    j=xdim-1
    do k=2, ydim-3
        km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=j-1; jp1=j+1; jp2=1; jp3=2

        !< Slopes
        sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
        sb = (Twz(j,kp1)-Twz(j,k)) / dy
        sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
        sd = (Twz(j,k)-Twz(j,km1)) / dy

        !< Mean winds on boundaries of cell
        uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
        vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
        uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
        vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

        !< Flux over cell coundaries
        !< fA
        if (uA > 0.) then
            sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
            fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else if (uA < 0.) then
            sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else
            fA = 0.
        end if

        !< fB
        if (vB > 0.) then
            sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else if (vB < 0.) then
            sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
            fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else
            fB = 0.
        end if

        !< fC
        if (uC > 0.) then
            sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else if (uC < 0.) then
            sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
            fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else
            fC = 0.
        end if

        !< fD
        if (vD > 0.) then
            sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
            fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else if (vD < 0.) then
            sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else
            fD = 0.
        end if

        dX_advec(j,k) = dt * ( fA - fB - fC + fD )

    end do

    j=xdim-2
    do k=2, ydim-3
        km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=j-1; jp1=j+1; jp2=j+2; jp3=1

        !< Slopes
        sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
        sb = (Twz(j,kp1)-Twz(j,k)) / dy
        sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
        sd = (Twz(j,k)-Twz(j,km1)) / dy

        !< Mean winds on boundaries of cell
        uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
        vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
        uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
        vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

        !< Flux over cell coundaries
        !< fA
        if (uA > 0.) then
            sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
            fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else if (uA < 0.) then
            sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else
            fA = 0.
        end if

        !< fB
        if (vB > 0.) then
            sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else if (vB < 0.) then
            sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
            fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else
            fB = 0.
        end if

        !< fC
        if (uC > 0.) then
            sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else if (uC < 0.) then
            sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
            fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else
            fC = 0.
        end if

        !< fD
        if (vD > 0.) then
            sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
            fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else if (vD < 0.) then
            sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else
            fD = 0.
        end if

        dX_advec(j,k) = dt * ( fA - fB - fC + fD )

    end do

    j=xdim-3
    do k=2, ydim-3
        km1=k-1; kp1=k+1; kp2=k+2; kp3=k+3; jm1=j-1; jp1=j+1; jp2=j+2; jp3=xdim

        !< Slopes
        sa = (Twz(j,k)-Twz(jm1,k)) / dx(k)
        sb = (Twz(j,kp1)-Twz(j,k)) / dy
        sc = (Twz(jp1,k)-Twz(j,k)) / dx(k)
        sd = (Twz(j,k)-Twz(j,km1)) / dy

        !< Mean winds on boundaries of cell
        uA = ( adv_flux(jm1,k)+adv_flux(j,k) ) / 2.
        vB = ( adv_no(j,kp1)+adv_no(j,k) ) / 2.
        uC = ( adv_flux(jp1,k)+adv_flux(j,k) ) / 2.
        vD = ( adv_no(j,km1)+adv_no(j,k) ) / 2.

        !< Flux over cell coundaries
        !< fA
        if (uA > 0.) then
            sa = minmod( (Twz(j,k)-Twz(jm1,k)) / dx(k), (Twz(jp1,k)-Twz(j,k)) / dx(k) )
            fA = (uA * Twz(jm1,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else if (uA < 0.) then
            sa = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fA = (uA * Twz(j,k) + 0.5*uA*sa*(dx(k) - uA*dt)) / dx(k)
        else
            fA = 0.
        end if

        !< fB
        if (vB > 0.) then
            sb = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fB = (vB * Twz(j,k) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else if (vB < 0.) then
            sb = minmod( (Twz(j,kp2)-Twz(j,kp1)) / dy, (Twz(j,kp3)-Twz(j,kp2)) / dy )
            fB = (vB * Twz(j,kp1) + 0.5*vB*sb*(dy - vB*dt)) / dy
        else
            fB = 0.
        end if

        !< fC
        if (uC > 0.) then
            sc = minmod( (Twz(jp1,k)-Twz(j,k)) / dx(k), (Twz(jp2,k)-Twz(jp1,k)) / dx(k) )
            fC = (uC * Twz(j,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else if (uC < 0.) then
            sc = minmod( (Twz(jp2,k)-Twz(jp1,k)) / dx(k), (Twz(jp3,k)-Twz(jp2,k)) / dx(k) )
            fC = (uC * Twz(jp1,k) + 0.5*uC*sc*(dx(k) - uC*dt)) / dx(k)
        else
            fC = 0.
        end if

        !< fD
        if (vD > 0.) then
            sd = minmod( (Twz(j,k)-Twz(j,km1)) / dy, (Twz(j,kp1)-Twz(j,k)) / dy )
            fD = (vD * Twz(j,km1) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else if (vD < 0.) then
            sd = minmod( (Twz(j,kp1)-Twz(j,k)) / dy, (Twz(j,kp2)-Twz(j,kp1)) / dy )
            fD = (vD * Twz(j,k) + 0.5*vD*sd*(dy - vD*dt)) / dy
        else
            fD = 0.
        end if

        dX_advec(j,k) = dt * ( fA - fB - fC + fD )

    end do

    contains
        real function minmod(a, b)
            real, intent(in)  :: a, b
            !real, intent(out) :: minmod

            if( abs(a) < abs(b) .AND. a*b > 0. ) then
                minmod = a
            else if ( abs(a) > abs(b) .AND. a*b > 0. ) then
                minmod = b
            else
                minmod = 0.
            end if

        end function minmod

end subroutine advection_finite_volume_pwl_slimit

