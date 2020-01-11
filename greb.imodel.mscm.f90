! ice sheet model for GREB
! start 26 Sep 2019
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_2d.cpp

subroutine ice_sheet(it, ionum, irec, mon, ice_H1, ice_T1, ice_Ts1, Ta1, dT_ocean, Fn_surf, &
&                    dq_rain, wz_vapor, z_surf, ice_H0, ice_T0, ice_Ts0)
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
&                        dlon, dlat, jday_mon, icestep, dt_ice
  USE mo_physics, ONLY: pi, cp_ice, rho_ice, d_ice_max, d_ice_mov,&
&                       Tl_ice2, jday, z_topo, garma_air
  implicit none
  character(len=37)  :: dirname
  real,parameter     :: gs_layer = 1/sqrt(3.)
  real,dimension(4)  :: zeta = (/-1.,-gs_layer,gs_layer,1./)
  integer                    :: ice_iunit              ! file unit of ice sheet data
  integer                    :: ionum                   ! file unit of GREB
  integer                    :: i, j, k, it, irec, mon  ! work variable for count
  real, dimension(xdim,ydim) :: Fn_surf                 ! surface net heat flux [J/m2]
  real, dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real, dimension(xdim,ydim) :: Ta1, ice_Ta1            ! air temperature [K]
  real, dimension(xdim,ydim) :: ice_snf                  ! snowfall accumulation rate [kg/m2/s]
  real, dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]
  real, dimension(xdim,ydim) :: wz_vapor                ! surface pressure change coefficient [1]
  real, dimension(xdim,ydim) :: z_surf                  ! surface height [m]

  real, dimension(xdim,ydim)     :: ice_H1, ice_H0, ice_zs, ice_mask
  real, dimension(xdim,ydim)     :: ice_Ts1, ice_Ts0
  real, dimension(xdim,ydim)     :: term_mass, dice_Ts, dT_hcorrect, term_hadv
  real, dimension(xdim,ydim+1)   :: crantx, cranty, fu, fv, ice_vx, ice_vy
  real, dimension(xdim,ydim,4)   :: ice_T1, ice_T0, ice_Tcoef
  real, dimension(xdim,ydim,4)   :: term_sig, term_dif, term_tadv, dice_T
  real, dimension(ydim)          :: lat, dxlat
  real                           :: deg, dx, dy, dyy
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)

  ice_zs    = z_topo+ice_H1

  term_hadv = 0.; term_mass = 0.; dice_Ts = 0.
  term_tadv = 0.; term_dif = 0.; term_sig = 0.; dice_T = 0.
  if (  (jday == 1 .or. jday == icestep + 1)                &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) ) then
      ! ice temperature profile estimate and vertical diffusion
      ice_T1(:,:,4) = ice_Ts1
      do j = 1,xdim
         do i = 1,ydim
            call ice_regression_coef(ice_T1(j,i,:), ice_Tcoef(j,i,:))
            call ice_temperature_diffusion(ice_T1(j,i,:), ice_H1(j,i), zeta, term_dif(j,i,:))
         end do
      end do
      ! mean velocity 
      ice_vx = 0.; ice_vy = 0.
      call ice_sheet_velocity_and_stress(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, 0., term_sig, 'mean')
      ! ice thickness advection
      call ice_advection(ice_vx,ice_vy,ice_H1,term_hadv, dyy, dxlat, lat)
      
      do k = 1,4
          ! velocity and stress term at 4 layers
          ice_vx = 0.; ice_vy = 0.
          call ice_sheet_velocity_and_stress(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta(k), term_sig(:,:,k), 'levs')
          term_sig(:,:,k) = dt_ice*term_sig(:,:,k)/(cp_ice*rho_ice)
          ! ice temperature advection
          call ice_advection(ice_vx,ice_vy,ice_T1(:,:,k),term_tadv(:,:,k), dyy, dxlat, lat)
      end do
  end if
  ! mass balance
  call ice_mass_balance(ice_H1, ice_Ts1, Ta1, dq_rain, wz_vapor, Fn_surf, dT_ocean, dice_Ts, term_mass)

  ! ice thickness equation
  ice_H0  = ice_H1 + term_hadv + term_mass
  where(z_topo == -0.1) ice_H0 = 0; ! may be changed
  ! ice surface temperature equation
  dT_hcorrect = (ice_H0-ice_H1)*garma_air
  ice_Ts0 = ice_Ts1 + dT_ocean + dice_Ts + dT_hcorrect
  ! ice temperature equation
  dice_T = term_tadv + term_dif + term_sig
  where(dice_T > Tl_ice2-ice_T1) dice_T = 0.9*(Tl_ice2-ice_T1) ! numeric stability
  ice_T0 = ice_T1 + dice_T
  do k = 1,4
      where(ice_H1 == 0 .and. ice_H0 > 0) ice_T0(:,:,k) = ice_Ts0 
  end do
  
  !print*,ice_H0(38,41), ice_vy(38,41),cranty(38,42),ice_vy(38,42)*86400
  ice_H1 = ice_H0
  ice_T1 = ice_T0
  z_surf = z_topo + ice_H0
  
  ! ice sheet output
  ice_iunit = 100 + ionum
  call ice_output(it, ice_iunit, irec, mon, ice_Ts0, ice_H0, z_surf, term_mass)

end subroutine

subroutine ice_mass_balance(ice_H1, ice_Ts1, Ta1, dq_rain, wz_vapor, Fn_surf, dT_ocean, dice_Ts, term_mass)
  USE mo_numerics, ONLY: xdim, ydim, dt
  USE mo_physics,  ONLY: z_topo, cap_ice, rho_ice, cp_ice, d_ice_max, cap_surf
  implicit none
  real, dimension(xdim,ydim) :: ice_Ts0                 ! ice surface temperature (forward) [K]
  real, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature (current) [K]
  real, dimension(xdim,ydim) :: Fn_surf                 ! surface net heat flux [J/m2]
  real, dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real, dimension(xdim,ydim) :: U_gtmlt                 ! ocean temperature change [K]

  real, dimension(xdim,ydim) :: ice_H0                  ! ice thickness (forward) [m]
  real, dimension(xdim,ydim) :: ice_H1                  ! ice thickness (current) [m]
  real, dimension(xdim,ydim) :: ice_snf                 ! snowfall accumulation rate [m]
  real, dimension(xdim,ydim) :: ice_melt                ! ice melting accumulation [m]
  real, dimension(xdim,ydim) :: ice_fus                 ! ice melting latent heat [W/m2]
  real, dimension(xdim,ydim) :: ice_cover               ! ice cover type [-1~0:partial ice sheet;0 land;0~1 partial sea ice]

  real, dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real, dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]
  real, dimension(xdim,ydim) :: wz_vapor                ! surface pressure change coefficient [1]

  real                       :: Tmin_limit              ! no very low Tsurf/Tatmoss;  numerical stability
  real, dimension(xdim,ydim) :: dice_Ts                 ! ice temperature total tendency [K]
  real, dimension(xdim,ydim) :: term_lat               !  ice temperature tendency due to fusion latent heat [K] 
  real, dimension(xdim,ydim) :: term_mass               ! ice thickness tendency due to mass balance [m] 

  dice_Ts = 0.; term_mass = 0.
  ! ice accumulation, thickness increase by snowfall
  call ice_accumulation(ice_snf, ice_Ts1, Ta1, dq_rain, wz_vapor) 
  ! ice fusion
  call ice_fusion( ice_fus, ice_melt, ice_Ts1, ice_H1, Fn_surf, dT_ocean) 
  
  ! mass balance 
  where(z_topo >= 0) term_mass = dt* (ice_snf + ice_melt)
  
  ! ice surface temperature
  !Tmin_limit = 40
  dice_Ts  = dt*(Fn_surf + ice_fus ) / cap_ice 
  ! thickness after fusion
  !  where(ice_H0 > 0.) cap_surf = cap_land                       ! ice sheet heat capacity

end subroutine

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_accumulation(ice_snf, ice_Ts1, Ta1, dq_rain, wz_vapor) 
!+++++++++++++++++++++++++++++++++++++++
! ice sheet : accumulation process
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics,  ONLY: r_qviwv, Tl_ice2, ice_svlm
  implicit none
  real, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature [K] 
  real, dimension(xdim,ydim) :: Ta1                     ! air temperature [K]
  real, dimension(xdim,ydim) :: ice_snf                  ! snowfall accumulation rate [kg/m2/s]
  real, dimension(xdim,ydim) :: dq_rain                 ! precipitation estimate (always negative) [m]
  real, dimension(xdim,ydim) :: wz_vapor                ! surface pressure change coefficient [1]
  ice_snf = 0.
  ! ice sheet accumulation from snowfall, kg/m2/s   
  where((ice_Ts1 <= Tl_ice2) .and. (Ta1 <= Tl_ice2)) ice_snf = - ice_svlm*dq_rain*r_qviwv*wz_vapor; 
end subroutine ice_accumulation

!+++++++++++++++++++++++++++++++++++++++
subroutine ice_fusion( ice_fus, ice_melt, ice_Ts1, ice_H1, Fn_surf, dT_ocean) 
!+++++++++++++++++++++++++++++++++++++++
  USE mo_numerics, ONLY: xdim, ydim, dt
  USE mo_physics,  ONLY: ci_latent, cp_ice, Tl_ice2, cap_ice, cap_surf, cap_land, rho_ice, d_ice_max
  implicit none
  real, dimension(xdim,ydim) :: ice_Tse                 ! ice surface temperature estimate [K]
  real, dimension(xdim,ydim) :: ice_Ts1                 ! ice surface temperature [K]
  real, dimension(xdim,ydim) :: ice_H1                  ! ice thickness [m]
  real, dimension(xdim,ydim) :: dT_ocean                ! ocean temperature change [K]
  real, dimension(xdim,ydim) :: ice_melt                ! ice melting rate [m/s]
  real, dimension(xdim,ydim) :: ice_fus                 ! ice melting latent heat [W/m2]
  real, dimension(xdim,ydim) :: U_gtmlt                 ! internal energy greater than melting point [J/m2]
  real, dimension(xdim,ydim) :: Lm_max                  ! potential energy of snow fusion [J/m2] 
  real, dimension(xdim,ydim) :: Fn_surf                 ! surface heat flux without fusion [W/m2] 

  cap_ice    = cap_surf
  ice_melt  = 0.; ice_fus   = 0.
 
  ! snow heat capacity
!  where((ice_H1 <  d_ice_max).and.(ice_H1 > 0.)) cap_ice = ice_H1*cp_ice*rho_ice ! heat capacity of snow [J/K/m^2]
  where(ice_H1 >= d_ice_max) cap_ice = d_ice_max*cp_ice*rho_ice ! heat capacity limitation of snow [J/K/m^2]
  
  ice_Tse  = ice_Ts1 + dT_ocean + dt*Fn_surf / cap_ice 
  where( ice_H1 > 0.)
      U_gtmlt = (ice_Tse - Tl_ice2) * cap_ice
      Lm_max  = ci_latent * rho_ice * ice_H1
  end where
  ! surface snow totally melts away
  where((ice_Ts1 >= Tl_ice2) .and. (U_gtmlt >  Lm_max) .and. (ice_H1 >0.))
        ice_melt    = - ice_H1 / dt
        ice_fus     = - (U_gtmlt - Lm_max) / dt
  end where
  ! surface snow partially melts
  where((ice_Ts1 >= Tl_ice2) .and. (U_gtmlt <= Lm_max) .and. (ice_H1 >0.))
        ice_melt    = - U_gtmlt / (dt*rho_ice*ci_latent)
        ice_fus     = - U_gtmlt / dt
  end where

end subroutine ice_fusion

subroutine ice_sheet_velocity_and_stress(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta, sigma_e, opt)
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics, ONLY: pi, rho_ice, grav, d_ice_mov 
  implicit none
  integer, external  :: cycle_ind
  character(len=4)   :: opt
  integer            :: i, j, k1, k2, k3, ind_zonal_real
  real               :: iceZ_jm1, iceZ_im1, iceZ_c, iceH_c,iceH_jm1, iceH_im1,dZ
  real               :: sighrz, sighrzjm, sighrzim
  real               :: T_layer, T_Gauss, T_Gaussjm, T_Gaussim, coef, coef_jm, coef_im
  real               :: cons_int, cons_intjm, cons_intim, term_poly, term_polyjm, term_polyim
  real               :: dZdx, dZdy
  real               :: deg, dyy, zeta, zh, zr
  real,dimension(4)  :: coef_poly  = (/1,3,3,1/)
  real,dimension(3)  :: gs_nodem   = (/0.8920, -0.4678, 0.1598/)
  real,dimension(3)  :: gs_weightm = (/3.4961,  2.5689, 0.3350/)
  real,dimension(3,4):: gs_node    = reshape((/-0.7746, 0.0000, 0.7746,&
                                               -0.8228,-0.1811, 0.5753,&
                                               -0.8540,-0.3060, 0.4100,&
                                               -0.8758,-0.3976, 0.2735/), (/3,4/))
  real,dimension(3,4):: gs_weight  = reshape((/ 0.5556, 0.8889, 0.5556,&
                                                0.8037, 0.9170, 0.2793,&
                                                1.2571, 1.1700, 0.2396,&
                                                2.063,  1.6736, 0.2637/),(/3,4/))
  real, dimension(xdim,ydim)     :: iceZ_ini, ice_H1, ice_H0, ice_Hx, ice_Hy, ice_zs
  real, dimension(xdim,ydim)     :: sigma_e
  real, dimension(xdim,ydim,4)   :: ice_Tcoef
  real, dimension(xdim,ydim+1)   :: ice_vx, ice_vy
  real, dimension(ydim)          :: lat, dxlat, ccx
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)
  

  sigma_e = 0.
  zh = (zeta+1)/2
  zr = 1-zh
  
  do i = 1,ydim+1
      do j = 1,xdim
          call meridion_shift(ice_H1, j, i, iceH_c)
          call meridion_shift(ice_H1, j, i-1, iceH_im1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_H1, ind_zonal_real, i, iceH_jm1)
          if(iceH_c < d_ice_mov .and. iceH_im1 < d_ice_mov .and. iceH_jm1 < d_ice_mov) cycle

          ! effective stress horizontal part, vertical part is left for intergal
          sighrz      = sigma_horizon(ice_zs, ice_H1, j  , i, dxlat, dyy)
          sighrzjm    = sigma_horizon(ice_zs, ice_H1, j-1, i, dxlat, dyy)
          sighrzim    = sigma_horizon(ice_zs, ice_H1, j, i-1, dxlat, dyy)
          ! vertical intergal part
          cons_int = 0.; cons_intjm = 0.; cons_intim = 0.
          do k1 = 1,4
              do k2 = 1,3
                  ! temperature at Gauss node
                  T_Gauss = 0.; T_Gaussjm = 0.; T_Gaussim = 0.; T_layer = 0.; term_poly = 0.
                  do k3 = 1,4
                      call meridion_shift(ice_Tcoef(:,:,k3), j, i, coef  )
                      T_layer   = T_layer   + coef*zeta**(k3-1)
                      if(opt == 'mean') then
                          T_Gauss   = T_Gauss   + coef*gs_nodem(k2)**(k3-1)
                      else
                          T_Gauss   = T_Gauss   + coef*gs_node(k2,k1)**(k3-1)
                      end if
                      ind_zonal_real = cycle_ind(j-1,xdim) 
                      call meridion_shift(ice_Tcoef(:,:,k3), ind_zonal_real, i, coef_jm)
                      if(opt == 'mean') then
                          T_Gaussjm = T_Gaussjm + coef_jm*gs_nodem(k2)**(k3-1)
                      else
                          T_Gaussjm = T_Gaussjm + coef_jm*gs_node(k2,k1)**(k3-1)
                      end if
                      call meridion_shift(ice_Tcoef(:,:,k3), j, i-1, coef_im)
                      if(opt == 'mean') then
                          T_Gaussim = T_Gaussim + coef_im*gs_nodem(k2)**(k3-1)
                      else
                          T_Gaussim = T_Gaussim + coef_im*gs_node(k2,k1)**(k3-1)
                      end if
                  end do
                  if(T_layer   > 273.15) T_layer   = 273.15
                  if(T_Gauss   > 273.15) T_Gauss   = 273.15
                  if(T_Gaussjm > 273.15) T_Gaussjm = 273.15
                  if(T_Gaussim > 273.15) T_Gaussim = 273.15
                  
                  term_poly = 0.; term_polyjm = 0.; term_polyim = 0.
                  ! polynomial coefficient
                  if(opt == 'mean') then
                      if(iceH_c   /= 0) term_poly   = (iceH_c  /2)**4*gs_weightm(k2)/(4.*2.)  
                      if(iceH_jm1 /= 0) term_polyjm = (iceH_jm1/2)**4*gs_weightm(k2)/(4.*2.)  
                      if(iceH_im1 /= 0) term_polyim = (iceH_im1/2)**4*gs_weightm(k2)/(4.*2.)  
                  elseif(zh > 0.) then
                      if(iceH_c   /= 0) term_poly   = coef_poly(k1)*(zh*iceH_c /2.)**4*(2.*zr/zh)**(4-k1)*gs_weight(k2,k1) 
                      if(iceH_jm1 /= 0) term_polyjm = coef_poly(k1)*(zh*iceH_jm1/2.)**4*(2.*zr/zh)**(4-k1)*gs_weight(k2,k1) 
                      if(iceH_im1 /= 0) term_polyim = coef_poly(k1)*(zh*iceH_im1/2.)**4*(2.*zr/zh)**(4-k1)*gs_weight(k2,k1) 
                  end if
                  cons_int    = cons_int   + consititutive_equation(T_Gauss  )*term_poly  *sighrz**2
                  cons_intjm  = cons_intjm + consititutive_equation(T_Gaussjm)*term_polyjm*sighrzjm**2
                  cons_intim  = cons_intim + consititutive_equation(T_Gaussim)*term_polyim*sighrzim**2
              end do
          end do

          ! ice surface gradient
          call meridion_shift(ice_zs, j, i, iceZ_c)
          call meridion_shift(ice_zs, j, i-1, iceZ_im1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_zs, ind_zonal_real, i, iceZ_jm1)
          if (i>ydim) then
              dZdx = 0.
          else
              dZ = (iceZ_c - iceZ_jm1)
              dZdx = dZ/dxlat(i)
          end if

          dZ = (iceZ_c - iceZ_im1)
          dZdy = dZ/dyy

          ice_vx(j,i)   = -2.*rho_ice*grav*dZdx*(cons_int+cons_intjm)/2.
          ice_vy(j,i)   = -2.*rho_ice*grav*dZdy*(cons_int+cons_intim)/2.
          sigma_e(j,i) = 2*(sighrz*zr*iceH_c)**4*consititutive_equation(T_layer)
      end do
  end do

  contains 
       real function consititutive_equation(T_Gauss)
       USE mo_physics, ONLY: A_fact, actene, R_gas, ice_Tcons, enh_fact
       implicit none
       real               :: T_Gauss
       T_Gauss = 0.
       if(T_Gauss > ice_Tcons .and. T_Gauss /= 0) then
           consititutive_equation =  A_fact(1)*exp(-actene(1)/R_gas/T_Gauss)*enh_fact
       else
           consititutive_equation =  A_fact(2)*exp(-actene(2)/R_gas/T_Gauss)*enh_fact
       end if
       end function 

       real function sigma_horizon(ice_zs, ice_H1, ind_zonal, ind_merid, dxlat, dy)
       USE mo_numerics, ONLY: xdim, ydim
       USE mo_physics, ONLY: rho_ice, grav 
       implicit none
       integer, external  :: cycle_ind
       integer            :: ind_zonal, ind_merid, ind_zonal_real
       real               :: iceZ_jp1, iceZ_jm1, iceZ_ip1, iceZ_im1
       real               :: iceH_jp1, iceH_jm1, iceH_ip1, iceH_im1
       real               :: dZdx_c, dZdy_c, dZ
       real               :: dx, dy
       real, dimension(ydim)      :: dxlat
       real, dimension(xdim,ydim) :: ice_zs, ice_H1
       ! zonal
       ind_zonal_real = cycle_ind(ind_zonal-1,xdim)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid, iceZ_jm1)
       call meridion_shift(ice_H1, ind_zonal_real, ind_merid, iceH_jm1)
       ind_zonal_real = cycle_ind(ind_zonal+1,xdim)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid, iceZ_jp1)
       call meridion_shift(ice_H1, ind_zonal_real, ind_merid, iceH_jp1)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid-1, iceZ_im1)
       call meridion_shift(ice_H1, ind_zonal_real, ind_merid-1, iceH_im1)
       call meridion_shift(ice_zs, ind_zonal_real, ind_merid+1, iceZ_ip1)
       call meridion_shift(ice_H1, ind_zonal_real, ind_merid+1, iceH_ip1)
       if   (ind_merid.gt.ydim) then
            dx = dxlat(2*ydim+1-ind_merid)
       else if(ind_merid.lt.1 ) then
            dx = dxlat(1-ind_merid)
       else 
            dx = dxlat(ind_merid)
       end if
       dZ     = (iceZ_jp1 - iceZ_jm1)
       dZdx_c = dZ/(2.*dx)
       dZ     = (iceZ_ip1 - iceZ_im1)
       dZdy_c = dZ/(2.*dy)
       sigma_horizon  = rho_ice*grav*sqrt(dZdx_c**2+dZdy_c**2)
       end function 

end subroutine

subroutine ice_advection(vx, vy, var1, term_adv, dyy, dxlat, lat)
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
                         dlon, dlat, dt_ice
  USE mo_physics, ONLY: pi
  implicit none
  integer, external  :: cycle_ind
  integer            :: i, j

  real, dimension(xdim,ydim)     :: var1, varx, vary, term_adv
  real, dimension(xdim,ydim+1)   :: crantx, cranty, fu, fv, vx, vy
  real, dimension(ydim)          :: lat, dxlat
  real                           :: dyy

  fu     = 0.; fv     = 0.

  do j = 1,xdim
        do i = 1,ydim
            crantx(j,i) = vx(j,i)*dt_ice/dxlat(i)
            cranty(j,i) = vy(j,i)*dt_ice/dyy
            varx(j,i) = var1(j,i)
            vary(j,i) = var1(j,i)*cos(lat(i)/180*pi)
        end do
        i = ydim+1
        cranty(j,i) = vy(j,i)*dt_ice/dyy
  end do

  call flux_operator(varx, vary, crantx, cranty, fu, fv)

  do j = 1,xdim
      do i = 1,ydim
          term_adv(j,i) = - (fu(cycle_ind(j+1,xdim),i)-fu(j,i)) - (fv(j,i+1)-fv(j,i))/cos(lat(i)/180*pi)
      end do   
  end do

end subroutine

subroutine ice_regression_coef(T1, Tcoef)
   implicit none
   integer             :: k1, k2
   real                :: sum_var
   real,dimension(4,4) :: rg_coef   = reshape((/-0.250, 0.750, 0.750, -0.250,&
                                             0.250,	-1.299,	1.299, -0.250,&
                                             0.750,	-0.750,	-0.750,	0.750,&
                                             -0.750,	1.299,	-1.299,	0.750/),(/4,4/))
   real, dimension(4)  :: T1, Tcoef
   do k1 = 1,4
       sum_var = 0
       do k2 = 1,4
            sum_var = sum_var + rg_coef(k2,k1)*T1(k2)
       end do
       Tcoef(k1) = sum_var
   end do
end subroutine

subroutine ice_temperature_diffusion(T1_ini, H1_ini, zeta, dif)
  USE mo_numerics, ONLY: xdim, ydim, dt_ice
  USE mo_physics, ONLY: d_ice_mov, ice_kappa, cp_ice, rho_ice
  implicit none
  integer, parameter   :: tstp_dif = 91 
  integer, parameter   :: kdim     = 4 
  
  integer              :: i, j, k, nday
  real                 :: dTdz2, T0_diff, T1_diff 
  real,dimension(kdim)    :: Tcoef, dT_diff, zeta, dzeta 

  real                 :: ice_H1, H1_ini
  real, dimension(kdim)   :: ice_T1, ice_T0, T1_ini
  real, dimension(kdim)   :: dif
  ice_H1 = H1_ini; ice_T1 = T1_ini
  dzeta(kdim) = zeta(kdim) - zeta(kdim-1)
  do k = 1,kdim-1
      dzeta(k) = zeta(k+1) - zeta(k)
  end do

  do nday = 1,tstp_dif
      call ice_regression_coef(ice_T1, Tcoef)
      if(ice_H1 < d_ice_mov) then 
          dT_diff = 0
      else
          do k = 2,kdim-1
              dTdz2 = 2*((ice_T1(k+1)-ice_T1(k))/dzeta(k)-(ice_T1(k)-ice_T1(k-1))/dzeta(k-1))/(dzeta(k+1)+dzeta(k)) 
              dT_diff(k) = dt_ice/tstp_dif*dTdz2*ice_kappa/(cp_ice*rho_ice)
          end do
          dT_diff(1) = dT_diff(2)
      end if
      ice_T0 = ice_T1 + dT_diff
      ice_T1 = ice_T0
  end do
  dif = ice_T0 - T1_ini

end subroutine

subroutine flux_operator(ice_Hx, ice_Hy, crantx, cranty, fu, fv)
  USE mo_numerics, ONLY: xdim, ydim
  implicit none
  integer, external             :: cycle_ind
  integer                       :: i, j, k 
  real                          :: iceH_indp2, iceH_indp1, iceH_ind, iceH_indm1, iceH_indm2 
  real                          :: xl, xr, dx1, dx2, dx3, fl, df, f6 
  real                          :: yl, yr, dy1, dy2, dy3
  integer                       :: crant_intx, crant_inty, adv_ind 
  real                          :: crant_frax, crant_fray 
  real, dimension(xdim,ydim)    :: ice_Hx, ice_Hy
  real, dimension(xdim,ydim+1)  :: ice_vx, ice_vy, crantx, cranty, fu, fv
  
  do i = 1,ydim
      do j = 1,xdim
      if(crantx(j,i)==0. .and. cranty(j,i)==0.) cycle
          ! x direction
          crant_intx = int(crantx(j,i))
          crant_frax = crantx(j,i) - float(crant_intx)
          
          ! integer flux
          if(crant_intx .ge. 1) then
              do k = 1,crant_intx
                 fu(j,i) = fu(j,i) + ice_Hx(cycle_ind(j-k,xdim),i) 
              end do
          else if(crant_intx .le. -1) then
              do k = 1,-crant_intx
                 fu(j,i) = fu(j,i) - ice_Hx(cycle_ind(j-1+k,xdim),i)
              end do
          endif
          ! fraction flux
          if(crantx(j,i) > 0) then
              adv_ind = cycle_ind(j-1-crant_intx, xdim)
          else
              adv_ind = cycle_ind(j-crant_intx, xdim)
          end if
          if(crant_frax > 0) then
              xl = 1-crant_frax; xr = 1
          else
              xl = 0; xr = -crant_frax
          end if
          iceH_indm2 = ice_Hx(cycle_ind(adv_ind-2,xdim), i); iceH_indm1 = ice_Hx(cycle_ind(adv_ind-1,xdim), i); 
          iceH_ind   = ice_Hx(adv_ind, i);
          iceH_indp2 = ice_Hx(cycle_ind(adv_ind+2,xdim), i); iceH_indp1 = ice_Hx(cycle_ind(adv_ind+1,xdim), i);
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dx1 = xr-xl; dx2 = xr*xr-xl*xl; dx3 = xr*xr*xr-xl*xl*xl;
          fu(j,i) = fu(j,i) + sign(fl*dx1+0.5*df*dx2+f6*(0.5*dx2-dx3/3.0), crant_frax);
          ! y direction
          crant_inty = int(cranty(j,i))
          crant_fray = cranty(j,i) - float(crant_inty)

          ! integer flux
          if(crant_inty .ge. 1) then
              do k = 1,crant_inty
                 call meridion_shift(ice_Hy, j, i-k, iceH_ind)
                 fv(j,i) = fv(j,i) + iceH_ind
              end do
          else if(crant_inty .le. -1) then
              do k = 1,-crant_inty
                 call meridion_shift(ice_Hy, j, i-1+k, iceH_ind)
                 fv(j,i) = fv(j,i) - iceH_ind
              end do
          endif
          ! fraction flux
          if(cranty(j,i) > 0) then
              adv_ind = i-1-crant_inty
          else
              adv_ind = i-crant_inty
          end if
          if(crant_fray > 0) then
              yl = 1-crant_fray; yr = 1
          else
              yl = 0; yr = -crant_fray
          end if
         
          if ((adv_ind>ydim-2).or.(adv_ind<2)) then
              call meridion_shift(ice_Hy, j, adv_ind-2, iceH_indm2)
              call meridion_shift(ice_Hy, j, adv_ind-1, iceH_indm1)
              call meridion_shift(ice_Hy, j, adv_ind  , iceH_ind  )
              call meridion_shift(ice_Hy, j, adv_ind+1, iceH_indp1)
              call meridion_shift(ice_Hy, j, adv_ind+2, iceH_indp2)
          else
              iceH_indm2 = ice_Hy(j,adv_ind-2); iceH_indm1 = ice_Hy(j,adv_ind-1); iceH_ind = ice_Hy(j,adv_ind);
              iceH_indp2 = ice_Hy(j,adv_ind+2); iceH_indp1 = ice_Hy(j,adv_ind+1);
          end if
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dy1 = yr-yl; dy2 = yr*yr-yl*yl; dy3 = yr*yr*yr-yl*yl*yl;
          fv(j,i) = fv(j,i) + sign(fl*dy1+0.5*df*dy2+f6*(0.5*dy2-dy3/3.0), crant_fray);         
      end do
  end do
  i = ydim + 1
  do j = 1,xdim
          ! y direction
          crant_inty = int(cranty(j,i))
          crant_fray = cranty(j,i) - float(crant_inty)

          ! integer flux
          if(crant_inty .ge. 1) then
              do k = 1,crant_inty
                 call meridion_shift(ice_Hy, j, i-k, iceH_ind)
                 fv(j,i) = fv(j,i) + iceH_ind
              end do
          else if(crant_inty .le. -1) then
              do k = 1,-crant_inty
                 call meridion_shift(ice_Hy, j, i-1+k, iceH_ind)
                 fv(j,i) = fv(j,i) - iceH_ind
              end do
          endif
          ! fraction fluy
          if(cranty(j,i) > 0) then
              adv_ind = i-1-crant_inty
          else
              adv_ind = i-crant_inty
          end if
          if(crant_fray > 0) then
              yl = 1-crant_fray; yr = 1
          else
              yl = 0; yr = -crant_fray
          end if
         
          if ((adv_ind>ydim-2).or.(adv_ind<2)) then
              call meridion_shift(ice_Hy, j, adv_ind-2, iceH_indm2)
              call meridion_shift(ice_Hy, j, adv_ind-1, iceH_indm1)
              call meridion_shift(ice_Hy, j, adv_ind  , iceH_ind  )
              call meridion_shift(ice_Hy, j, adv_ind+1, iceH_indp1)
              call meridion_shift(ice_Hy, j, adv_ind+2, iceH_indp2)
          else
              iceH_indm2 = ice_Hy(j,adv_ind-2); iceH_indm1 = ice_Hy(j,adv_ind-1); iceH_ind = ice_Hy(j,adv_ind);
              iceH_indp2 = ice_Hy(j,adv_ind+2); iceH_indp1 = ice_Hy(j,adv_ind+1);
          end if
          call ppm(iceH_indm2, iceH_indm1, iceH_ind, iceH_indp1, iceH_indp2, fl, df, f6)
          dy1 = yr-yl; dy2 = yr*yr-yl*yl; dy3 = yr*yr*yr-yl*yl*yl;
          fv(j,i) = fv(j,i) + sign(fl*dy1+0.5*df*dy2+f6*(0.5*dy2-dy3/3.0), crant_fray);         
  end do
end subroutine

subroutine meridion_shift(ice_H1, ind_zonal, ind_merid, ice_nounb)
 USE mo_numerics, ONLY: xdim, ydim
 implicit none
 integer, external          :: cycle_ind
 integer                    :: ind_zonal, ind_merid
 real                       :: ice_nounb
 real, dimension(xdim,ydim) :: ice_H1
 if   (ind_merid.gt.ydim) then
      ice_nounb = ice_H1(cycle_ind(ind_zonal+xdim/2,xdim),2*ydim+1-ind_merid)
 else if(ind_merid.lt.1 ) then
      ice_nounb = ice_H1(cycle_ind(ind_zonal+xdim/2,xdim), 1-ind_merid)
 else 
      ice_nounb = ice_H1(ind_zonal,ind_merid)
 end if
end subroutine

subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)
  implicit none
  real :: fm2, fm1, f, fp1, fp2 
  real :: dfl, df, dfr, fl, fr, f6
  dfl     = mismatch(fm2, fm1, f)
  df      = mismatch(fm1, f  , fp1)
  dfr     = mismatch(f  , fp1, fp2)
  fl      = 0.5*(fm1+f)+(dfl-df)/6.
  fr      = 0.5*(fp1+f)+(df-dfr)/6.

  fl = f-sign(min(abs(df), abs(fl-f)), df);
  fr = f+sign(min(abs(df), abs(fr-f)), df);
  f6 = 6*f-3*(fl+fr);
  df = fr-fl;

  contains
      real function mismatch(fm1, f, fp1)  ! half the slope calculation in (Colella & Woodward, 1984)
      real fm1, f, fp1, df, dfMin, dfMax
      if ((fp1-f)*(f-fm1) .gt. 0) then
          df       = (fp1 - fm1)*0.5
          dfMin    = f - min(fm1, f , fp1)
          dfMax    = max(fm1, f, fp1) - f
          mismatch = sign(min(abs(df), 2*dfMin, 2*dfMax), df) ! monontonic
          !mismatch = sign(min(abs(df), 2*f), df) ! positive defined
      else
          mismatch = 0
      end if
      end function mismatch
end

integer function cycle_ind(x,xdim)
  implicit none
  integer :: x, xdim
  if(x < 1) then
    cycle_ind = xdim + x
   else if(x > xdim) then
    cycle_ind = x - xdim
   else
    cycle_ind = x
   endif
end function

subroutine ice_output(it, ice_iunit, irec, mon, ice_Ts0, ice_H0, z_surf, term_mass)
!+++++++++++++++++++++++++++++++++++++++
! ice sheet : output file
  USE mo_numerics,     ONLY: xdim, ydim, jday_mon, ndt_days, nstep_yr, time_scnr &
&                          , time_ctrl, ireal, dt
  USE mo_physics,      ONLY: jday, log_exp, r_qviwv, wz_vapor, cap_ice
  use mo_diagnostics,  ONLY: ice_Tsmm, ice_Tsmn_ctrl, ice_Hmn_ctrl, ice_mask_ctrl
  implicit none
  real, external              :: gmean
  real, dimension(xdim,ydim)  :: ice_H0, ice_Ts0, ice_mask, melt, term_mass, U_gtmlt, term_massmm
  real, dimension(xdim,ydim)  :: z_surf 
  integer, parameter          :: nvar = 3              ! number of output variable
  integer                     :: it,irec,mon,iyrec     ! work variable for count, controled by subroutine output
  integer                     :: ndm                   ! total time for mean calculation
  integer                     :: ice_iunit            ! written file uint

  ! diagnostics: monthly means
  ice_mask = 0.; term_massmm = 0.
  where(ice_H0 > 0.) ice_mask = 1.
  ice_Tsmm = ice_Tsmm+ice_Ts0
  term_massmm = term_massmm + term_mass
  
! control output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. ice_iunit == 201 ) then
     ndm=jday_mon(mon)*ndt_days
     if (it/float(ndt_days)  > 365*(time_ctrl-1)) then
         if (log_exp .eq. 1 .or. log_exp .eq. 310 ) then
         write(ice_iunit,rec=nvar*irec+1)  ice_Tsmm/ndm
         write(ice_iunit,rec=nvar*irec+2)  ice_H0
         write(ice_iunit,rec=nvar*irec+3)  z_surf
         else
         ice_Tsmn_ctrl(:,:,mon)  = ice_Tsmm/ndm
         ice_Hmn_ctrl(:,:,mon)   = ice_H0
         ice_mask_ctrl(:,:,mon)  = z_surf
         end if
     end if
     ice_Tsmm=0.
  end if

! scenario output
  if (       jday == sum(jday_mon(1:mon))                   &
&      .and. it/float(ndt_days) == nint(it/float(ndt_days)) &
&      .and. ice_iunit == 202 ) then

     ndm=jday_mon(mon)*ndt_days
     write(ice_iunit,rec=nvar*irec+1)  ice_Tsmm/ndm
     write(ice_iunit,rec=nvar*irec+2)  ice_H0
     write(ice_iunit,rec=nvar*irec+3)  z_surf

     write(203,rec=               12*iyrec+mon) gmean(ice_Tsmm/ndm - ice_Tsmn_ctrl(:,:,mon))
     write(203,rec=1*12*time_scnr+12*iyrec+mon) gmean(ice_H0 -ice_Hmn_ctrl(:,:,mon))
     write(203,rec=2*12*time_scnr+12*iyrec+mon) gmean(z_surf -ice_mask_ctrl(:,:,mon))
     ice_Tsmm=0.; term_massmm=0.
  end if

end subroutine ice_output
