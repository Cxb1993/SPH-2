PROGRAM grid2d

  USE types
  USE vars, ONLY : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  USE grid, ONLY : grid_create,grid_destroy
  USE kernels, ONLY : kernel_select,kernel_support
  IMPLICIT NONE
  
  INTEGER(I4B), PARAMETER :: ndim = 2
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: n = (/ 10, 10 /)
  REAL(DP), DIMENSION(ndim),  PARAMETER :: xmin=(/ 0.0_dp, 0.0_dp /), xmax=(/ 1.0_dp, 1.0_dp /)
  REAL(DP), DIMENSION(ndim), PARAMETER :: grav=(/ 0.0_dp, 0.0_dp /) 
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: periodic = (/ 1, 1 /)
  REAL(DP), PARAMETER :: h = 0.12_dp, rho0 = 1.0_dp, gam = 7.0_dp 
  REAL(DP), PARAMETER :: cs = 10.0_dp, p0 = rho0*cs**2/gam, mu0 = 1.D-6
  CHARACTER(LEN=80) :: kernelname = 'm4',dir = '/home/tmattner/sph/',job = 'test'
  CHARACTER(LEN=120) :: filnam
  INTEGER(I4B) :: nsteps=100
 
  INTEGER(I4B) :: i,j,k,np
  REAL(DP) :: dx(ndim),dt,r

  np = n(1)*n(2)
  dx = (xmax - xmin)/n
  dt = MIN(0.1_dp*MINVAL(dx)/cs,0.1_dp*MINVAL(dx)**2*rho0/mu0)

  OPEN(10,FILE='input')
  WRITE(10,'(A)')TRIM(dir)
  WRITE(10,'(A)')TRIM(job)
  WRITE(10,*)ndim
  WRITE(10,*)np
  WRITE(10,*)xmin
  WRITE(10,*)xmax
  WRITE(10,*)periodic
  WRITE(10,*)h
  WRITE(10,'(A)')TRIM(kernelname)
  WRITE(10,*)rho0
  WRITE(10,*)gam
  WRITE(10,*)p0
  WRITE(10,*)grav
  WRITE(10,*)dt
  WRITE(10,*)nsteps
  WRITE(10,*)0
  WRITE(10,*)0.0_dp
  CLOSE(10)

  WRITE(filnam,'(3A,I3.3)')'mkdir ',TRIM(dir),TRIM(job)
  CALL SYSTEM(TRIM(filnam))
  

  CALL vars_create(ndim,np)

  DO i=1,n(1)
     DO j=1,n(2)
        k = (i-1)*n(2) + j
        x(1,k) = dx(1)*(i-1) + 0.5_dp*dx(1)
        x(2,k) = dx(2)*(j-1) + 0.5_dp*dx(2)
     END DO
  END DO   

  CALL vars_params(h,gam,rho0,p0,grav)
  CALL kernel_select(kernelname)
  CALL grid_create(xmin,xmax,kernel_support(h),periodic)
  m = 1.0_dp
  CALL vars_density
  r = SUM(rho)/np ! average density for m = 1
  m = rho0/r      ! average density = rho0
  rho = rho/r     ! density equation is linear

  CALL RANDOM_NUMBER(u)
  rho = rho + 0.0001_dp*u(1,:)
  u = 0.0_dp
  
  
  CALL vars_pressure
  mu = mu0
  bc = 0
  
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),0
  CALL vars_write(TRIM(filnam))
  CALL grid_destroy
  CALL vars_destroy
  

END PROGRAM grid2d
