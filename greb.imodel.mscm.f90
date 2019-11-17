! ice sheet model for GREB
! start 26 Sep 2019
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_2d.cpp
program ice_sheet_model
  implicit none
  character(len=100) :: filename
  integer, external  :: cycle_ind
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndt_days  = 24*3600           ! number of timesteps per day
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: ndt_yr    = ndays_yr*ndt_days        ! number of timesteps per day
  real,parameter     :: gs_layer = 1./sqrt(3.)
  real,dimension(4)  :: zeta = (/-1.,-gs_layer,gs_layer,1./)
  real               :: dt        = ndt_yr/2           ! time step [s]
  integer            :: ireal     = 4
  integer            :: i, j, k, irec, nstep_end 
  real               :: iceH_indp2, iceH_indp1, iceH_ind, iceH_indm1, iceH_indm2 
  real,parameter     :: pi        = 3.1416 
  real, dimension(xdim,ydim)     :: iceH_ini, ice_H1, ice_H0, ice_Hx, ice_Hy, ice_zs, z_topo
  real, dimension(xdim,ydim+1,4) :: ice_vx3, ice_vy3, vy3, vx3
  real, dimension(xdim,ydim+1)   :: crantx, cranty, fu, fv, ice_vx, ice_vy, ice_vmx, ice_vmy
  real, dimension(xdim,ydim,4)   :: ice_T1, ice_T0, ice_Tcoef
  real, dimension(ydim)      :: lat, dxlat
  real    :: deg, dx, dy, dyy
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/)

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  !deg = pi/180.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)

  nstep_end = 10*96
  
  open(301,file='ice_scheme_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(302,file='ice_scheme_test_zs.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(306,file='ice_scheme_test_thk.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  open(303,file='ice_scheme_test_tmp.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim*4)
  open(304,file='ice_scheme_test_vx.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  open(305,file='ice_scheme_test_vy.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  filename = '/Users/zxie0012/Documents/ice_sheet_model/model/my_ice_model/sico_out/zs.bin'
  open(11,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  read(11,rec=1) ice_zs
  filename = '/Users/zxie0012/Documents/ice_sheet_model/model/my_ice_model/sico_out/thk.bin'
  open(12,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  read(12,rec=1) ice_H1
  filename = '/Users/zxie0012/Documents/ice_sheet_model/model/my_ice_model/sico_out/temp.bin'
  open(15,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim*4)
  read(15,rec=1) ice_T1
  filename = '/Users/zxie0012/Documents/ice_sheet_model/model/my_ice_model/sico_out/vx.bin'
  open(13,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  read(13,rec=1) vx3
  filename = '/Users/zxie0012/Documents/ice_sheet_model/model/my_ice_model/sico_out/vy.bin'
  open(14,file=filename,ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*(ydim+1)*4)
  read(14,rec=1) vy3
  z_topo = ice_zs - ice_H1
  call ice_regression_coef(ice_T1, ice_Tcoef)
  do k = 2,4
      call ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta(k), 'levs')
      ice_vx3(:,:,k) = ice_vx
      ice_vy3(:,:,k) = ice_vy
  end do
  write(302,rec=1) ice_zs
  write(303,rec=1) ice_T1
  write(304,rec=1) ice_vx3
  write(305,rec=1) ice_vy3
  write(306,rec=1) ice_H1

  do irec = 1,nstep_end
  
  call ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta(k), 'mean')

  ! outer operator
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
          ice_H0(j,i) = ice_H1(j,i) - (fu(cycle_ind(j+1,xdim),i)-fu(j,i)) - (fv(j,i+1)-fv(j,i))/cos(lat(i)/180*pi)
      end do   
  end do
  
  write(301,rec=irec) ice_H0
  
  ice_H1 = ice_H0
  ice_zs = z_topo+ice_H0
  end do
end

subroutine ice_sheet_velocity(ice_zs,ice_H1,ice_Tcoef, dyy, dxlat, ice_vx, ice_vy, zeta, opt)
  implicit none
  integer, external  :: cycle_ind
  character(len=4)   :: opt
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  real,parameter     :: g         = 9.8
  integer            :: i, j, k1, k2, k3, ind_zonal_real
  real               :: iceZ_jm1, iceZ_im1, iceZ_c, iceH_c,iceH_jm1, iceH_im1,dZ
  real               :: sigma, sigmajm, sigmaim
  real               :: T_Gauss, T_Gaussjm, T_Gaussim, coef, coef_jm, coef_im
  real               :: cons_int, cons_intjm, cons_intim, term_poly 
  real               ::  dZdx, dZdy
  real               :: deg, dyy, zeta, zh, zr
  real,parameter     :: pi         = 3.1416 
  real,parameter     :: rho_ice    = 910.
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
  real, dimension(xdim,ydim,4)   :: ice_Tcoef
  real, dimension(xdim,ydim+1)   :: ice_vx, ice_vy
  real, dimension(ydim)          :: lat, dxlat, ccx
  integer, dimension(ydim)       :: ilat = (/(i,i=1,ydim)/)

  ice_vx = 0.
  ice_vy = 0.
  zh = (zeta+1)/2
  zr = 1-zh
  
  do i = 1,ydim+1
      do j = 1,xdim
          call meridion_shift(ice_H1, j, i, iceH_c)
          call meridion_shift(ice_H1, j, i-1, iceH_im1)
          ind_zonal_real = cycle_ind(j-1,xdim) 
          call meridion_shift(ice_H1, ind_zonal_real, i, iceH_jm1)
          if(iceH_c .eq. 0) cycle

          ! effective stress horizontal part, vertical part is left for intergal
          sigma      = sigma_horizon(ice_zs, ice_H1, j  , i, dxlat, dyy)
          sigmajm    = sigma_horizon(ice_zs, ice_H1, j-1, i, dxlat, dyy)
          sigmaim    = sigma_horizon(ice_zs, ice_H1, j, i-1, dxlat, dyy)
          ! if(i.eq.33.and.j.eq.71) print*,sigma
          ! vertical intergal part
          cons_int = 0.; cons_intjm = 0.; cons_intim = 0.
          do k1 = 1,4
              do k2 = 1,3
                  ! temperature at Gauss node
                  T_Gauss = 0.; T_Gaussjm = 0.; T_Gaussim = 0.
                  do k3 = 1,4
                      call meridion_shift(ice_Tcoef(:,:,k3), j, i, coef  )
                      if(opt .eq. 'mean') then
                          T_Gauss   = T_Gauss   + coef*gs_nodem(k2)**(k3-1)
                      else
                          T_Gauss   = T_Gauss   + coef*gs_node(k2,k1)**(k3-1)
                      end if
                      if(T_Gauss .eq. 0) print*,i,j, iceZ_c
                      ind_zonal_real = cycle_ind(j-1,xdim) 
                      call meridion_shift(ice_Tcoef(:,:,k3), ind_zonal_real, i, coef_jm)
                      if(opt .eq. 'mean') then
                          T_Gaussjm = T_Gaussjm + coef_jm*gs_nodem(k2)**(k3-1)
                      else
                          T_Gaussjm = T_Gaussjm + coef_jm*gs_node(k2,k1)**(k3-1)
                      end if
                      call meridion_shift(ice_Tcoef(:,:,k3), j, i-1, coef_im)
                      if(opt .eq. 'mean') then
                          T_Gaussim = T_Gaussim + coef_im*gs_nodem(k2)**(k3-1)
                      else
                          T_Gaussim = T_Gaussim + coef_im*gs_node(k2,k1)**(k3-1)
                      end if
                  end do
                  if(T_Gauss   > 273.15) T_Gauss   = 273.15
                  if(T_Gaussjm > 273.15) T_Gaussjm = 273.15
                  if(T_Gaussim > 273.15) T_Gaussim = 273.15
                  
                  ! constitutive equation
                  if(opt .eq. 'mean') then
                      term_poly  = (iceH_c/2)**5*gs_weightm(k2)/4./iceH_c  
                  else
                      term_poly  = coef_poly(k1)*(zh*iceH_c/2.)**4*(2.*zr/zh)**(4-k1)*gs_weight(k2,k1) 
                  end if
                  cons_int    = cons_int   + constitutive_equation(T_Gauss  )*term_poly*sigma**2
                  cons_intjm  = cons_intjm + constitutive_equation(T_Gaussjm)*term_poly*sigmajm**2
                  cons_intim  = cons_intim + constitutive_equation(T_Gaussim)*term_poly*sigmaim**2
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
              !if ( dZ.gt. iceH_c  ) dZ =  iceH_c   ! gradient larger than thickness at j+1
              !if ( dZ.lt.-iceH_jm1) dZ = -iceH_jm1 ! gradient larger than thickness at j-1
              dZdx = dZ/dxlat(i)
          end if

          dZ = (iceZ_c - iceZ_im1)
          !if ( dZ.gt. iceH_c  ) dZ =  iceH_c   ! gradient larger than thickness at j+1
          !if ( dZ.lt.-iceH_im1) dZ = -iceH_im1 ! gradient larger than thickness at j-1
          dZdy = dZ/dyy

          ice_vx(j,i)  = -2.*rho_ice*g*dZdx*(cons_int+cons_intjm)/2.
          ice_vy(j,i)  = -2.*rho_ice*g*dZdy*(cons_int+cons_intim)/2.
          !if(i.eq.33.and.j.eq.71) print*,ice_vx(j,i)*365*86400
      end do
  end do

  contains 
       real function constitutive_equation(T_Gauss)
       implicit none
       real,parameter     :: R_air     = 8.314
       real               :: T_Gauss
       if(T_Gauss > 263.15) then
           constitutive_equation =  1.96e3*exp(-1.39e5/R_air/T_Gauss)*3
       else
           constitutive_equation =  3.985e-13*exp(-6e4/R_air/T_Gauss)*3
       end if
       end function 

       real function sigma_horizon(ice_zs, ice_H1, ind_zonal, ind_merid, dxlat, dy)
       implicit none
       integer, external  :: cycle_ind
       integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
       integer            :: ind_zonal, ind_merid, ind_zonal_real
       real,parameter     :: g         = 9.8
       real,parameter     :: rho_ice   = 910.
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
       ! meridianal
       ind_zonal_real = cycle_ind(ind_zonal,xdim)
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
       ! if ( dZ.gt. iceH_jp1) dZ =  iceH_jp1 ! gradient larger than thickness at j+1
       ! if ( dZ.lt.-iceH_jm1) dZ = -iceH_jm1 ! gradient larger than thickness at j-1
       dZdx_c = dZ/(2.*dx)
       dZ     = (iceZ_ip1 - iceZ_im1)
       ! if ( dZ.gt. iceH_ip1) dZ =  iceH_ip1 ! gradient larger than thickness at i+1
       ! if ( dZ.lt.-iceH_im1) dZ = -iceH_im1 ! gradient larger than thickness at i-1
       dZdy_c = dZ/(2.*dy)
       sigma_horizon  = rho_ice*g*sqrt(dZdx_c**2+dZdy_c**2)
       end function 

end subroutine

subroutine ice_regression_coef(ice_T1, ice_Tcoef)
   implicit none
   integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
   integer            :: i, j, k1, k2
   real               :: sum_var
   real,parameter     :: pi        = 3.1416 
   real,dimension(4,4):: rg_coef   = reshape((/-0.250, 0.750, 0.750, -0.250,&
                                             0.250,	-1.299,	1.299, -0.250,&
                                             0.750,	-0.750,	-0.750,	0.750,&
                                             -0.750,	1.299,	-1.299,	0.750/),(/4,4/))
   real, dimension(xdim,ydim,4)   :: ice_T1, ice_Tcoef
   do i = 1,ydim
      do j = 1,xdim
         do k1 = 1,4
             sum_var = 0
             do k2 = 1,4
                  sum_var = sum_var + rg_coef(k2,k1)*ice_T1(j,i,k2)
             end do
             ice_Tcoef(j,i,k1) = sum_var
         end do
      end do
   end do
end subroutine

subroutine flux_operator(ice_Hx, ice_Hy, crantx, cranty, fu, fv)
  implicit none
  integer, external             :: cycle_ind
  integer, parameter            :: xdim = 96, ydim = 48          ! field dimensions
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
 implicit none
 integer, external          :: cycle_ind
 integer, parameter         :: xdim = 96, ydim = 48          ! field dimensions
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
