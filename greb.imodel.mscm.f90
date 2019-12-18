! ice sheet model for GREB
! start 26 Sep 2019
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_2d.cpp

include "greb.model.mscm.f90"


program ice_sheet_model
  USE mo_numerics, ONLY: xdim, ydim, ndt_days, ndays_yr, ireal, &
                         dlon, dlat
  USE mo_physics, ONLY: pi, cp_ice, rho_ice, ice_kappa, z_topo, d_ice_max, Tl_ice2
  implicit none
  real,external      :: gmean
  integer, external  :: cycle_ind
  character(len=120) :: filename
  character(len=37)  :: dirname
  integer, parameter :: ndt_yr    = ndays_yr*ndt_days        ! number of timesteps per day
  real,parameter     :: gs_layer = 1./sqrt(3.)
  real,dimension(4)  :: zeta = (/-1.,-gs_layer,gs_layer,1./)
  real               :: dt        = 365*86400/2           ! time step [s]
  integer            :: i, j, k,  nrec, irec, nstep_end,year 
  real               :: iceH_indp2, iceH_indp1, iceH_ind, iceH_indm1, iceH_indm2
  real               :: dTdz2, T0_diff, T1_diff 

  real, dimension(xdim,ydim)     :: ice_H1, ice_H0, ice_Hx, ice_Hy, ice_zs, ice_mask
  real, dimension(xdim,ydim)     :: ice_Hm1, ice_Hmn 
  real, dimension(xdim,ydim)     :: term_mass
  real, dimension(xdim,ydim)     :: ice_Tx, ice_Ty, ice_Ts1, ice_Ts0, ice_Tsmn
  real, dimension(xdim,ydim+1,4) :: ice_vx3, ice_vy3, vy3, vx3
  real, dimension(xdim,ydim+1)   :: crantx, cranty, fu, fv, ice_vx, ice_vy, ice_vmx, ice_vmy
  real, dimension(xdim,ydim,4)   :: ice_T1, ice_T0, ice_Tcoef
  real, dimension(xdim,ydim,4)   :: term_sig, term_dif, term_adv, dice_T
  real, dimension(ydim)      :: lat, dxlat
  real    :: deg, dx, dy, dyy
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/)

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)

  nstep_end = 19999
  year = 1

  
 ! initialization 
  open(19,file='../input/global.topography.bin',      		ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  read(19,rec=1)  z_topo
 
  dirname = "/Volumes/YMI/research_data/GREB_10kyr"
  filename = dirname // '/scenario_ice_sheet.exp-310.solar.231K.bin'
  open(301,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  filename = dirname // '/output/ice_thickness.bin'
  open(302,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  filename = dirname //  '/output/ice_temperature.bin'
  open(303,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim*4)
  filename = dirname // '/output/ice_velocity_x.bin'
  open(304,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  filename = dirname // '/output/ice_velocity_y.bin'
  open(305,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  
  nrec = 0
  do k=1,6
      read(301,rec=nrec+1) ice_Tsmn
      read(301,rec=nrec+2) ice_Hmn
      read(301,rec=nrec+3) ice_mask
      nrec = nrec + 3
  end do
  ice_Ts1   = ice_Tsmn
  ice_H1    = ice_Hmn
  ice_Hm1   = ice_Hmn

  do k = 1,4
      ice_T1(:,:,k) = ice_Ts1
  end do
  
  do irec = 1,nstep_end
  
  do k=1,6
      read(301,rec=nrec+1) ice_Tsmn
      read(301,rec=nrec+2) ice_Hmn
      read(301,rec=nrec+3) ice_mask
      nrec = nrec + 3
  end do

  ice_zs    = z_topo+ice_H1
  ice_Ts1   = ice_Tsmn
  term_mass = ice_Hmn - ice_Hm1
  ice_Hm1   = ice_Hmn
  ice_T1(:,:,4) = ice_Ts1
  do j = 1,xdim
     do i = 1,ydim
        call ice_regression_coef(ice_T1(j,i,:), ice_Tcoef(j,i,:))
        call ice_temperature_diffusion(ice_T1(j,i,:), ice_H1(j,i), zeta, term_dif(j,i,:))
     end do
  end do
  ! thickness eqation (outer operator)
  ! velocity 
  ice_vx = 0.
  ice_vy = 0.
  call ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, 0., term_sig, 'mean',irec)
  
  fu     = 0.
  fv     = 0.

  do j = 1,xdim
        do i = 1,ydim
            crantx(j,i) = ice_vx(j,i)*dt/dxlat(i)
            cranty(j,i) = ice_vy(j,i)*dt/dyy
            ice_Hx(j,i) = ice_H1(j,i)
            ice_Hy(j,i) = ice_H1(j,i)*cos(lat(i)/180*pi)
        end do
        i = ydim+1
        cranty(j,i) = ice_vy(j,i)*dt/dyy
  end do

  call flux_operator(ice_Hx, ice_Hy, crantx, cranty, fu, fv)

  do j = 1,xdim
      do i = 1,ydim
          ice_H0(j,i) = ice_H1(j,i) - (fu(cycle_ind(j+1,xdim),i)-fu(j,i)) - (fv(j,i+1)-fv(j,i))/cos(lat(i)/180*pi) &
                      + term_mass(j,i)
      end do   
  end do
  ! temperature eqation (outer operator)
  do k = 1,4
      ! velocity 
      ice_vx = 0.
      ice_vy = 0.
      call ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta(k), term_sig(:,:,k), 'levs',irec)
      ice_vx3(:,:,k) = ice_vx
      ice_vy3(:,:,k) = ice_vy
      
      fu     = 0.
      fv     = 0.
    
      do j = 1,xdim
            do i = 1,ydim
                crantx(j,i) = ice_vx(j,i)*dt/dxlat(i)
                cranty(j,i) = ice_vy(j,i)*dt/dyy
                ice_Tx(j,i) = ice_T1(j,i,k)
                ice_Ty(j,i) = ice_T1(j,i,k)*cos(lat(i)/180*pi)
            end do
            i = ydim+1
            cranty(j,i) = ice_vy(j,i)*dt/dyy
      end do
    
      call flux_operator(ice_Tx, ice_Ty, crantx, cranty, fu, fv)
    
      do j = 1,xdim
          do i = 1,ydim
              term_adv(j,i,k) =  - (fu(cycle_ind(j+1,xdim),i)-fu(j,i)) - (fv(j,i+1)-fv(j,i))/cos(lat(i)/180*pi)
              term_sig(j,i,k) = dt*term_sig(j,i,k)/(cp_ice*rho_ice)
          end do   
      end do

  end do
  dice_T = term_adv + term_dif + term_sig
  where(dice_T > Tl_ice2-ice_T1) dice_T = 0.9*(Tl_ice2-ice_T1) ! numeric stability
  ice_T0 = ice_T1 + dice_T
  
  !print*,ice_H0(38,41), ice_vy(38,41),cranty(38,42),ice_vy(38,42)*86400
   if(mod(irec,2) == 1) then
   print *, "YEAR ", "ice surface temperature ", "ice thickness [m] ", "Greenland ", "Antarctic ", "Tibetan" !TB
   print *, year, gmean(ice_T0(:,:,4))-273.15,gmean(ice_H1) !TB
   year = year + 1
   end if

  write(302,rec=irec) ice_H0
  write(303,rec=irec) ice_T0
  write(304,rec=irec) ice_vx3
  write(305,rec=irec) ice_vy3
  
  ice_H1 = ice_H0
  ice_T1 = ice_T0

  end do
end

subroutine ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta, term_sig, opt, irec)
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics, ONLY: pi, rho_ice, grav, d_ice_max 
  implicit none
  integer, external  :: cycle_ind
  character(len=4)   :: opt
  integer            :: i, j, k1, k2, k3, ind_zonal_real,irec
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
  real, dimension(xdim,ydim)     :: term_sig
  real, dimension(xdim,ydim,4)   :: ice_Tcoef
  real, dimension(xdim,ydim+1)   :: ice_vx, ice_vy
  real, dimension(ydim)          :: lat, dxlat, ccx
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)

  term_sig = 0.
  zh = (zeta+1)/2
  zr = 1-zh
  
  do i = 1,ydim+1
      do j = 1,xdim
          call meridion_shift(ice_H1, j, i, iceH_c)
          call meridion_shift(ice_H1, j, i-1, iceH_im1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_H1, ind_zonal_real, i, iceH_jm1)
          if(iceH_c < d_ice_max .and. iceH_im1 < d_ice_max .and. iceH_jm1 < d_ice_max) cycle

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
 !             if ( dZ.gt. iceH_c  ) dZ =  iceH_c   ! gradient larger than thickness at j+1
 !             if ( dZ.lt.-iceH_jm1) dZ = -iceH_jm1 ! gradient larger than thickness at j-1
              dZdx = dZ/dxlat(i)
          end if

          dZ = (iceZ_c - iceZ_im1)
 !         if ( dZ.gt. iceH_c  ) dZ =  iceH_c   ! gradient larger than thickness at j+1
 !         if ( dZ.lt.-iceH_im1) dZ = -iceH_im1 ! gradient larger than thickness at j-1
          dZdy = dZ/dyy

          ice_vx(j,i)   = -2.*rho_ice*grav*dZdx*(cons_int+cons_intjm)/2.
          ice_vy(j,i)   = -2.*rho_ice*grav*dZdy*(cons_int+cons_intim)/2.
          term_sig(j,i) = 2*(sighrz*zr*iceH_c)**4*consititutive_equation(T_layer)
      end do
  end do

  contains 
       real function consititutive_equation(T_Gauss)
       USE mo_physics, ONLY: A_fact, actene, R_gas, ice_Tcons, enh_fact
       implicit none
       real               :: T_Gauss
       if(T_Gauss > ice_Tcons) then
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
 !!       if ( dZ.gt. iceH_jp1) dZ =  iceH_jp1 ! gradient larger than thickness at j+1
 !       if ( dZ.lt.-iceH_jm1) dZ = -iceH_jm1 ! gradient larger than thickness at j-1
       dZdx_c = dZ/(2.*dx)
       dZ     = (iceZ_ip1 - iceZ_im1)
 !       if ( dZ.gt. iceH_ip1) dZ =  iceH_ip1 ! gradient larger than thickness at i+1
 !       if ( dZ.lt.-iceH_im1) dZ = -iceH_im1 ! gradient larger than thickness at i-1
       dZdy_c = dZ/(2.*dy)
       sigma_horizon  = rho_ice*grav*sqrt(dZdx_c**2+dZdy_c**2)
       end function 

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
  USE mo_numerics, ONLY: xdim, ydim
  USE mo_physics, ONLY: d_ice_max, ice_kappa, cp_ice, rho_ice
  implicit none
  integer, parameter   :: tstp_dif = 500 
  integer, parameter   :: kdim     = 4 
  real                 :: dt        = 365*86400/2           ! time step [s]
  
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
      if(ice_H1 < d_ice_max) then 
          dT_diff = 0
      else
          do k = 2,kdim-1
              dTdz2 = 2*((ice_T1(k+1)-ice_T1(k))/dzeta(k)-(ice_T1(k)-ice_T1(k-1))/dzeta(k-1))/(dzeta(k+1)+dzeta(k)) 
              dT_diff(k) = dt/tstp_dif*dTdz2*ice_kappa/(cp_ice*rho_ice)
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
