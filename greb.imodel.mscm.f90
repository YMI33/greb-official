
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
  integer            :: j, irec, nstep_end,i 
  real,parameter     :: pi        = 3.1416 
  real,parameter     :: dx        = 2*pi/(xdim)
  real, dimension(xdim,ydim) :: iceH_ini, ice_H2, ice_H1, ice_H0, ice_vx, diceH_crcl

  nstep_end = 96*10
 
  open(301,file='ice_scheme_test.bin',ACCESS='DIRECT',FORM='UNFORMATTED', RECL=ireal*xdim*ydim)
  ice_vx = dx/ndt_yr
  do j = 1,xdim
  do i = 1,ydim
      iceH_ini(j,i) = exp(-dx*(j-48)**2)
  end do
  end do
  ice_H1 = iceH_ini
  ice_H2 = ice_H1

  do irec = 1,nstep_end
  do i = 1,ydim
      j = 1
      diceH_crcl(j,i) = ice_vx(j+1,i)*ice_H1(j+1,i) - ice_vx(xdim,i)*ice_H1(xdim,i)
      j = xdim
      diceH_crcl(j,i) = ice_vx(1,i)*ice_H1(1,i) - ice_vx(j-1,i)*ice_H1(j-1,i)
      do j = 2, xdim-1
          diceH_crcl(j,i) = ice_vx(j+1,i)*ice_H1(j+1,i) - ice_vx(j-1,i)*ice_H1(j-1,i) 
      end do
!      dt = -sum(diceH_crcl*ice_H2)*dx**2/(sum(diceH_crcl*diceH_crcl)*dx**2)*2.
  end do
      dt = ndt_yr/2
      diceH_crcl = diceH_crcl/(2*dx)
      ice_H0 = ice_H2 - dt*diceH_crcl
  write(301,rec=irec) ice_H0
  
  if(irec<2) then
       ice_H1 = ice_H0
  else
       ice_H2 = ice_H1
       ice_H1 = ice_H0
  end if

  end do

end


