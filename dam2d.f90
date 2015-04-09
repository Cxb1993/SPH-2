PROGRAM grid2d

  USE types
  USE vars, ONLY : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  USE grid, ONLY : grid_create,grid_destroy
  USE kernels, ONLY : kernel_select,kernel_support
  IMPLICIT NONE
  
  INTEGER(I4B), PARAMETER :: ndim = 2
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: n = (/ 20, 20 /)
  REAL(DP), DIMENSION(ndim), PARAMETER :: xmin=(/ 0.0_dp, 0.0_dp /), xmax=(/ 1.0_dp, 1.0_dp /)
  REAL(DP), DIMENSION(ndim), PARAMETER :: grav=(/ 0.0_dp, -1.0_dp /) 
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: periodic = (/ 1, 0 /)
  REAL(DP) :: h,cs,p0
  REAL(DP), PARAMETER :: rho0 = 1.0_dp, gam = 7.0_dp, mu0 = 1.D-8
!  REAL(DP) :: h = 0.5_dp*0.13_dp, rho0 = 1.0_dp, gam = 7.0_dp
!  REAL(DP), PARAMETER :: cs = 10.0_dp, p0 = 300.0_dp, mu0 = 1.D-3
  CHARACTER(LEN=80) :: kernelname = 'm4',dir = '/home/tmattner/sph/',job = 'dam'
  CHARACTER(LEN=120) :: filnam
  INTEGER(I4B) :: nsteps=100
 
  INTEGER(I4B) :: i,j,k,np
  REAL(DP) :: dx(ndim),dt,r
  
  np = n(1)*n(2)
  dx = (xmax - xmin)/n
!  dt = MIN(0.1_dp*MINVAL(dx)/cs,0.1_dp*MINVAL(dx)**2*rho0/mu0)

  p0 = 200.0_dp/7.0_dp*rho0*ABS(grav(2))*(xmax(2) - xmin(2))
  cs = SQRT(gam*p0/rho0)
  dt = MIN(0.5_dp*MINVAL(dx)/cs,0.5_dp*MINVAL(dx)**2*rho0/mu0)
  h = 1.3_dp*MAXVAL(dx)

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

  u = 0.0_dp
  bc = 0
  DO j=1,n(2)
     DO i=1,n(1)
        k = (j-1)*n(2) + i
        x(1,k) = dx(1)*(i-1) + 0.5_dp*dx(1)
        x(2,k) = dx(2)*(j-1) + 0.5_dp*dx(2)
        IF (j < 5) THEN
           bc(k) = 1
 !       ELSE
 !          x(1,k) = 0.5*x(1,k) + 0.25_dp
        END IF
     END DO
  END DO 

  CALL vars_params(h,gam,rho0,p0,grav)
  CALL kernel_select(kernelname)
  CALL grid_create(xmin,xmax,kernel_support(h),periodic)
  m = rho0*dx(1)*dx(2)
  rho = rho0
!  DO i=1,np
!     IF (bc(i) == 0) THEN
!        rho(i) = rho0*(1.0_dp + (gam - 1.0_dp)/gam*rho0*grav(2)*(x(2,i) - xmax(2))/p0)**(1.0_dp/(gam - 1.0_dp))
!     END IF
!  END DO
  CALL vars_pressure
  mu = mu0

  
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),0
  CALL vars_write(TRIM(filnam))
  CALL grid_destroy
  CALL vars_destroy
  
  
END PROGRAM grid2d
