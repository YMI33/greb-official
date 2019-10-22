! ice sheet model for GREB
! start 26 Sep 2019
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_2d.cpp
program ice_sheet_model
  implicit none
  integer, external  :: cycle_ind
  integer, parameter :: xdim = 96, ydim = 48          ! field dimensions
  integer, parameter :: ndt_days  = 24*3600           ! number of timesteps per day
  integer, parameter :: ndays_yr  = 365               ! number of days per year
  integer, parameter :: ndt_yr    = ndays_yr*ndt_days        ! number of timesteps per day
  real               :: dt        = ndt_yr/2           ! time step [s]
  integer            :: ireal     = 4
  integer            :: i, j, k, irec, nstep_end 
  real               :: iceH_indp2, iceH_indp1, iceH_ind, iceH_indm1, iceH_indm2 
  real,parameter     :: pi        = 3.1416 
  real, dimension(xdim,ydim)     :: iceH_ini, ice_H1, ice_H0, ice_Hx, ice_Hy
  real, dimension(xdim,ydim+1)   :: ice_vx, ice_vy, crantx, cranty, fu, fv
  real, dimension(ydim)      :: lat, dxlat, ccx
  real    :: deg, dx, dy, dyy
  real, parameter    :: dlon      = 360./xdim         ! linear increment in lon
  real, parameter    :: dlat      = 180./ydim         ! linear increment in lat
  integer, dimension(ydim)   :: ilat = (/(i,i=1,ydim)/)

  deg = 2.*pi*6.371e6/360.;   ! length of 1deg latitude [m]
  deg = pi/180.;   ! length of 1deg latitude [m]
  dx = dlon; dy=dlat; dyy=dy*deg
  lat = dlat*ilat-dlat/2.-90.;  dxlat=dx*deg*cos(2.*pi/360.*lat)

  nstep_end = 10*96
 
  open(301,file='ice_scheme_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  do i = 1,ydim
      do j = 1,xdim
          ice_vx(j,i) = dxlat(i)/ndt_yr
          ice_vx(j,i) = 0 
      end do
  end do
!  ice_vy   = 0.
  do i = 1,ydim+1
      do j = 1,xdim/4
         ice_vy(j,i)   = -dyy/ndt_yr
      end do
      do j = xdim/4+1,xdim*3/4
         ice_vy(j,i)   = dyy/ndt_yr
      end do
      do j = 3*xdim/4+1,xdim
         ice_vy(j,i)   = -dyy/ndt_yr
      end do
  end do
  iceH_ini = 0. 
  do j = 1,xdim
      do i = 1,ydim
          iceH_ini(j,i) = exp(-((i-24)/10.)**2-((j-48)/10.)**2)
          !iceH_ini(j,i) = sin(j*2*pi/96)+1
      end do
  end do
  ice_H1 = iceH_ini

  do irec = 1,nstep_end
  
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
  end do
end

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
