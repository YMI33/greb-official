
! ice sheet model for GREB
! start 26 Sep 2019
! reference code : https://github.com/dongli/IAP-CGFD/blob/master/advection/ffsl/main_1d.cpp
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
  integer            :: indm2, indm1, indp1, indp2 
  real               :: xl, xr, dx1, dx2, dx3, fl, df, f6 
  integer            :: crant_int, adv_ind 
  real               :: crant, crant_fra 
  real,parameter     :: pi        = 3.1416 
  real,parameter     :: dx        = 2*pi/(xdim)
  real, dimension(xdim,ydim)   :: iceH_ini, ice_H2, ice_H1, ice_H0, ice_vx, diceH_crcl
  real, dimension(xdim,ydim)   :: fu

  nstep_end = 96*10
 
  open(301,file='ice_scheme_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  ice_vx = 4*dx/ndt_yr
  do j = 1,xdim
  do i = 1,ydim
      iceH_ini(j,i) = exp(-dx*(j-48)**2)
  end do
  end do
  ice_H1 = iceH_ini
  ice_H2 = ice_H1

  do irec = 1,nstep_end
  fu     = 0.
  
  do i = 1,ydim
      do j = 1,xdim
          crant  = ice_vx(j,i)*dt/dx
          crant_int = int(crant)
          crant_fra = crant - float(crant_int)
          ! integer flux
          if(crant_int .ge. 1) then
              do k = 1,crant_int
                 fu(j,i) = fu(j,i) + ice_H1(cycle_ind(j+1-k,xdim),i) 
              end do
          else if(crant_int .le. -1) then
              do k = 1,-crant_int
                 fu(j,i) = fu(j,i) - ice_H1(cycle_ind(j+k,xdim),i)
              end do
          endif
          ! fraction flux
          if(crant > 0) then
              adv_ind = cycle_ind(j-crant_int, xdim)
          else
              adv_ind = cycle_ind(j+1-crant_int, xdim)
          end if
          if(crant_fra > 0) then
              xl = 1-crant_fra; xr = 1
          else
              xl = 0; xr = -crant_fra
          end if
          
          indm1=cycle_ind(adv_ind-1,xdim);indm2=cycle_ind(adv_ind-2,xdim);
          indp1=cycle_ind(adv_ind+1,xdim);indp2=cycle_ind(adv_ind+2,xdim);
          call ppm(ice_H1(indm2,i),ice_H1(indm1,i),ice_H1(adv_ind,i),ice_H1(indp1,i),ice_H1(indp2,i),fl,df,f6)
          dx1 = xr-xl; dx2 = xr*xr-xl*xl; dx3 = xr*xr*xr-xl*xl*xl;
          fu(j,i) = fu(j,i) + sign(fl*dx1+0.5*df*dx2+f6*(0.5*dx2-dx3/3.0), crant_fra);
      end do
  end do
  do i = 1,ydim
      do j = 1,xdim
          ice_H0(j,i) = ice_H1(j,i) - (fu(j,i)-fu(cycle_ind(j-1,xdim),i))
      end do   
  end do

  write(301,rec=irec) ice_H0
  
  ice_H1 = ice_H0
!  if(irec<2) then 
!     ice_H1 = ice_H0
!  else
!     ice_H2 = ice_H1
!     ice_H1 = ice_H0
!  end if

  end do

end

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
      df       = (fp1 - fm1)*0.5
      dfMin    = f - min(fm1, f , fp1)
      dfMax    = max(fm1, f, fp1) - f
      mismatch = sign(min(abs(df)*0.5, dfMin, dfMax), df)

      end function mismatch

end

integer function cycle_ind(x,xdim)
  implicit none
  integer :: x, xdim, xdiff
  if(x < 1) then
      cycle_ind = xdim + x
  else if(x > xdim) then
      cycle_ind = x - xdim
  else
      cycle_ind = x
  endif
end function





