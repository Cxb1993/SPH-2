program taylorgreen3d

  use types
  use vars, only : vars_create, vars_destroy, vars_write, vars_params, vars_density, vars_pressure
  use vars, only : x, u, m, rho, p, mu, bc
  use grid, only : grid_create, grid_destroy
  use kernels, only : kernel_select, kernel_support
  implicit none
  
  integer(i4b), parameter                  :: ndim = 3
  integer(i4b), dimension(ndim), parameter :: n = (/ 64, 64, 64 /)
  real(dp), dimension(ndim), parameter     :: xmin = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  real(dp), dimension(ndim), parameter     :: xmax = (/ 2.0_dp*pi, 2.0_dp*pi, 2.0_dp*pi /)
  real(dp), dimension(ndim), parameter     :: grav = (/ 0.0_dp, 0.0_dp, 0.0_dp /) 
  integer(i4b), dimension(ndim), parameter :: periodic = (/ 1, 1, 1 /)
  real(dp), parameter                      :: h = 1.2_dp*2.0_dp*pi/64.0_dp
  real(dp), parameter                      :: rho0 = 1.0_dp
  real(dp), parameter                      :: gam = 7.0_dp
  real(dp), parameter                      :: mu0 = 0.01d0
  real(dp), parameter                      :: cs = 10.0_dp, p0 = rho0*cs**2/gam
  character(len=80)                        :: kernelname = 'm4'
  character(len=80)                        :: dir = '/mnt/disk2/tmattner/SPH/'
  character(len=80)                        :: job = 'taylorgreen3d'
  integer(i4b)                             :: nsteps=100
 
   character(len=120) :: filnam
  integer(i4b)        :: i, j, k, kk, np
  real(dp)            :: dx(ndim), dt, r, x0, x1
  
  np = n(1)*n(2)*n(3)
  dx = (xmax - xmin)/n
  dt = min(0.1_dp*minval(dx)/cs, 0.1_dp*minval(dx)**2*rho0/mu0)

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

  write(filnam, '(3a,i3.3)') 'mkdir ', trim(dir), trim(job)
  call system(trim(filnam))
 
  call vars_create(ndim,np)

  u = 0.0_dp
  bc = 0
  do k=1,n(3)
     do j=1,n(2)
        do i=1,n(1)
           kk = (k-1)*n(1)*n(2) + (j-1)*n(1) + i
           x(1,kk) = dx(1)*(i-1) + 0.5_dp*dx(1)  
           x(2,kk) = dx(2)*(j-1) + 0.5_dp*dx(2)
           x(3,kk) = dx(3)*(k-1) + 0.5_dp*dx(3)
           u(1,kk) = sin(x(1,kk))*cos(x(2,kk))
           u(2,kk) = -cos(x(1,kk))*sin(x(2,kk))
           u(3,kk) = 0.0_dp
        end do
     end do
  end do

  call vars_params(h, gam, rho0, p0, grav)
  call kernel_select(kernelname)
  call grid_create(xmin, xmax, kernel_support(h), periodic)
  m = 1.0_dp
  call vars_density
  r = sum(rho)/np ! average density for m = 1
  m = rho0/r      ! average density = rho0
  rho = rho/r     ! density equation is linear
  call vars_pressure
  mu = mu0

  write(filnam, '(4a,i3.3)') trim(dir), trim(job), '/', trim(job), 0
  call vars_write(trim(filnam))
  call grid_destroy
  call vars_destroy
  
end program taylorgreen3d
