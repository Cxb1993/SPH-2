program grid2d

  use types
  use vars, only : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  use grid, only : grid_create,grid_destroy
  use kernels, only : kernel_select,kernel_support
  implicit none
  
  integer(i4b), parameter :: ndim = 2
  integer(i4b), parameter :: n = 20
  real(dp), dimension(ndim), parameter :: xmin=(/ 0.0_dp, 0.0_dp /), xmax=(/ 1.0_dp, 1.0_dp /)
  real(dp), dimension(ndim), parameter :: grav=(/ 0.0_dp, 0.0_dp /) 
  integer(i4b), dimension(ndim), parameter :: periodic = (/ 0, 0 /)
  real(dp), parameter :: dr = 0.49_dp/(n-1)
  real(dp), parameter :: h = 1.3_dp*dr, rho0 = 1.0_dp, gam = 7.0_dp
  real(dp), parameter :: cs = 10.0_dp, p0 = rho0*cs**2/gam, mu0 = 0.1_dp
  character(len=80) :: kernelname = 'm4',dir = '/mnt/disk2/tmattner/sph/',job = 'spindown2d'
  character(len=120) :: filnam
  integer(i4b) :: nsteps=100
 
  integer(i4b) :: i,j,k,np,ntheta
  real(dp) :: dtheta,theta,dt,r,x0,x1

! calculate number of particles

  np = 1
  do i=2,n
     np = np + int(2.0_dp*pi_d*(i-1))
  end do
  dt = min(0.1_dp*dr/cs,0.1_dp*dr**2*rho0/mu0)

  open(10,file='input')
  write(10,'(a)')trim(dir)
  write(10,'(a)')trim(job)
  write(10,*)ndim
  write(10,*)np
  write(10,*)xmin
  write(10,*)xmax
  write(10,*)periodic
  write(10,*)h
  write(10,'(a)')trim(kernelname)
  write(10,*)rho0
  write(10,*)gam
  write(10,*)p0
  write(10,*)grav
  write(10,*)dt
  write(10,*)nsteps
  write(10,*)0
  write(10,*)0.0_dp
  close(10)

  write(filnam,'(3a,i3.3)')'mkdir ',trim(dir),trim(job)
  call system(trim(filnam))
 
  call vars_create(ndim,np)

  u = 0.0_dp
  bc = 0
  k = 0

  do i=1,n
     if (i == 1) then
        k = k + 1
        x(1,k) = 0.5_dp
        x(2,k) = 0.5_dp
        u(1,k) = 0.0_dp
        u(2,k) = 0.0_dp
     else
        r = dr*(i-1)
        ntheta = int(2.0_dp*pi_d*(i-1))
        dtheta = 2.0_dp*pi_d/ntheta
        do j=1,ntheta
           k = k + 1
           theta = j*dtheta
           x(1,k) = r*cos(theta) + 0.5_dp
           x(2,k) = r*sin(theta) + 0.5_dp
           if (i > n-3) then
              u(1,k) = 0.0_dp
              u(2,k) = 0.0_dp
              bc(k) = 1
           else
              u(1,k) = -r*sin(theta)
              u(2,k) = r*cos(theta)
           end if
        end do
     end if
  end do

 
  call vars_params(h,gam,rho0,p0,grav)
  call kernel_select(kernelname)
  call grid_create(xmin,xmax,kernel_support(h),periodic)
  m = 1.0_dp
  call vars_density
  r = sum(rho)/np ! average density for m = 1
  m = rho0/r      ! average density = rho0
  rho = rho/r     ! density equation is linear
  call vars_pressure
  mu = mu0

  write(filnam,'(4a,i3.3)')trim(dir),trim(job),'/',trim(job),0
  call vars_write(trim(filnam))
  call grid_destroy
  call vars_destroy
  
  
end program grid2d
