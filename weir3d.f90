program weir3d

  use types
  use vars, only : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  use grid, only : grid_create,grid_destroy
  use kernels, only : kernel_select,kernel_support
  implicit none
  
  integer(i4b), parameter :: ndim = 3
  real(dp), dimension(ndim), parameter :: grav=(/ 0.0_dp, -10.0_dp, 0.0_dp /) 
  integer(i4b), dimension(ndim), parameter :: periodic = (/ 0, 0, 1 /)
  real(dp), parameter :: dx = 0.1_dp
  real(dp), dimension(ndim) :: xmin, xmax
  real(dp), parameter :: rho0 = 1.0_dp, gam = 7.0_dp
  real(dp), parameter :: mu0 = 1.d-8
  character(len=80) :: kernelname = 'm4',dir = '/mnt/disk2/tmattner/sph/',job = 'weir3d'
  character(len=120) :: filnam
  integer(i4b) :: nsteps=100
 
  integer(i4b) :: i,j,k,np,cnt
  real(dp) :: h,dt,p0,cs

  cnt = 0

  np = 134100
  call vars_create(ndim,np)

  do i=1,4
     do j=1,100
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1)
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,50
     do j=1,4
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*4
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do
     
  do i=1,4
     do j=1,25
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*54
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,50
     do j=1,4
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*58
            x(2,cnt) = dx*(j-1) + dx*21
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,4
     do j=1,25
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*108
            x(2,cnt) = dx*(j-1) + dx*21
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,50
     do j=1,4
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*112
            x(2,cnt) = dx*(j-1) + dx*42
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,4
     do j=1,58
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*162
            x(2,cnt) = dx*(j-1) + dx*42
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        end do
     end do
  end do

  do i=1,50
     do j=1,25
        do k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*112
            x(2,cnt) = dx*(j-1) + dx*46
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 0
        end do
     end do
  end do

  np = cnt
  print *,np

  h = 1.3_dp*dx
  do i=1,2
     xmin(i) = minval(x(i,:)) - 1.d-8
     xmax(i) = maxval(x(i,:)) + 1.d-8
  end do
  xmin(3) = minval(x(3,:)) - 0.5_dp*dx ! for periodic bcs
  xmax(3) = maxval(x(3,:)) + 0.5_dp*dx ! for periodic bcs
  
  p0 = 200.0_dp/7.0_dp*rho0*abs(grav(2))*25.0_dp*dx
  cs = sqrt(gam*p0/rho0)
  dt = min(0.5_dp*dx/cs,0.5_dp*dx**2*rho0/mu0)

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

  u = 0.0_dp

  call vars_params(h,gam,rho0,p0,grav)
  call kernel_select(kernelname)
  call grid_create(xmin,xmax,kernel_support(h),periodic)
  
  m = rho0*dx**3
  rho = rho0
!  call vars_density
!  do i=1,np
!     if (bc(i) == 0) then
!        rho(i) = rho0*(1.0_dp + (gam - 1.0_dp)/gam*rho0*grav(2)*(x(2,i) - dx*70)/p0)**(1.0_dp/(gam - 1.0_dp))
!     end if
!  end do
  call vars_pressure
  mu = mu0

  
  write(filnam,'(4a,i3.3)')trim(dir),trim(job),'/',trim(job),0
  call vars_write(trim(filnam))
  call grid_destroy
  call vars_destroy
  
  
end program weir3d
