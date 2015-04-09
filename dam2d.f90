program grid2d

  use types
  use vars, only : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  use grid, only : grid_create,grid_destroy
  use kernels, only : kernel_select,kernel_support
  implicit none
  
  integer(i4b), parameter :: ndim = 2
  integer(i4b), dimension(ndim), parameter :: n = (/ 20, 20 /)
  real(dp), dimension(ndim), parameter :: xmin=(/ 0.0_dp, 0.0_dp /), xmax=(/ 1.0_dp, 1.0_dp /)
  real(dp), dimension(ndim), parameter :: grav=(/ 0.0_dp, -1.0_dp /) 
  integer(i4b), dimension(ndim), parameter :: periodic = (/ 1, 0 /)
  real(dp) :: h,cs,p0
  real(dp), parameter :: rho0 = 1.0_dp, gam = 7.0_dp, mu0 = 1.d-8
!  real(dp) :: h = 0.5_dp*0.13_dp, rho0 = 1.0_dp, gam = 7.0_dp
!  real(dp), parameter :: cs = 10.0_dp, p0 = 300.0_dp, mu0 = 1.d-3
  character(len=80) :: kernelname = 'm4',dir = '/home/tmattner/sph/',job = 'dam'
  character(len=120) :: filnam
  integer(i4b) :: nsteps=100
 
  integer(i4b) :: i,j,k,np
  real(dp) :: dx(ndim),dt,r
  
  np = n(1)*n(2)
  dx = (xmax - xmin)/n
!  dt = min(0.1_dp*minval(dx)/cs,0.1_dp*minval(dx)**2*rho0/mu0)

  p0 = 200.0_dp/7.0_dp*rho0*abs(grav(2))*(xmax(2) - xmin(2))
  cs = sqrt(gam*p0/rho0)
  dt = min(0.5_dp*minval(dx)/cs,0.5_dp*minval(dx)**2*rho0/mu0)
  h = 1.3_dp*maxval(dx)

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
  do j=1,n(2)
     do i=1,n(1)
        k = (j-1)*n(2) + i
        x(1,k) = dx(1)*(i-1) + 0.5_dp*dx(1)
        x(2,k) = dx(2)*(j-1) + 0.5_dp*dx(2)
        if (j < 5) then
           bc(k) = 1
 !       else
 !          x(1,k) = 0.5*x(1,k) + 0.25_dp
        end if
     end do
  end do 

  call vars_params(h,gam,rho0,p0,grav)
  call kernel_select(kernelname)
  call grid_create(xmin,xmax,kernel_support(h),periodic)
  m = rho0*dx(1)*dx(2)
  rho = rho0
!  do i=1,np
!     if (bc(i) == 0) then
!        rho(i) = rho0*(1.0_dp + (gam - 1.0_dp)/gam*rho0*grav(2)*(x(2,i) - xmax(2))/p0)**(1.0_dp/(gam - 1.0_dp))
!     end if
!  end do
  call vars_pressure
  mu = mu0

  
  write(filnam,'(4a,i3.3)')trim(dir),trim(job),'/',trim(job),0
  call vars_write(trim(filnam))
  call grid_destroy
  call vars_destroy
  
  
end program grid2d
