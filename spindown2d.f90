PROGRAM grid2d

  USE types
  USE vars, ONLY : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  USE grid, ONLY : grid_create,grid_destroy
  USE kernels, ONLY : kernel_select,kernel_support
  IMPLICIT NONE
  
  INTEGER(I4B), PARAMETER :: ndim = 2
  INTEGER(I4B), PARAMETER :: n = 20
  REAL(DP), DIMENSION(ndim), PARAMETER :: xmin=(/ 0.0_dp, 0.0_dp /), xmax=(/ 1.0_dp, 1.0_dp /)
  REAL(DP), DIMENSION(ndim), PARAMETER :: grav=(/ 0.0_dp, 0.0_dp /) 
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: periodic = (/ 0, 0 /)
  REAL(DP), PARAMETER :: dr = 0.49_dp/(n-1)
  REAL(DP), PARAMETER :: h = 1.3_dp*dr, rho0 = 1.0_dp, gam = 7.0_dp
  REAL(DP), PARAMETER :: cs = 10.0_dp, p0 = rho0*cs**2/gam, mu0 = 0.1_dp
  CHARACTER(LEN=80) :: kernelname = 'm4',dir = '/mnt/disk2/tmattner/SPH/',job = 'spindown2d'
  CHARACTER(LEN=120) :: filnam
  INTEGER(I4B) :: nsteps=100
 
  INTEGER(I4B) :: i,j,k,np,ntheta
  REAL(DP) :: dtheta,theta,dt,r,x0,x1

! Calculate number of particles

  np = 1
  DO i=2,n
     np = np + INT(2.0_dp*pi_d*(i-1))
  END DO
  dt = MIN(0.1_dp*dr/cs,0.1_dp*dr**2*rho0/mu0)

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
  k = 0

  DO i=1,n
     IF (i == 1) THEN
        k = k + 1
        x(1,k) = 0.5_dp
        x(2,k) = 0.5_dp
        u(1,k) = 0.0_dp
        u(2,k) = 0.0_dp
     ELSE
        r = dr*(i-1)
        ntheta = INT(2.0_dp*pi_d*(i-1))
        dtheta = 2.0_dp*pi_d/ntheta
        DO j=1,ntheta
           k = k + 1
           theta = j*dtheta
           x(1,k) = r*COS(theta) + 0.5_dp
           x(2,k) = r*SIN(theta) + 0.5_dp
           IF (i > n-3) THEN
              u(1,k) = 0.0_dp
              u(2,k) = 0.0_dp
              bc(k) = 1
           ELSE
              u(1,k) = -r*SIN(theta)
              u(2,k) = r*COS(theta)
           END IF
        END DO
     END IF
  END DO

 
  CALL vars_params(h,gam,rho0,p0,grav)
  CALL kernel_select(kernelname)
  CALL grid_create(xmin,xmax,kernel_support(h),periodic)
  m = 1.0_dp
  CALL vars_density
  r = SUM(rho)/np ! average density for m = 1
  m = rho0/r      ! average density = rho0
  rho = rho/r     ! density equation is linear
  CALL vars_pressure
  mu = mu0

  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),0
  CALL vars_write(TRIM(filnam))
  CALL grid_destroy
  CALL vars_destroy
  
  
END PROGRAM grid2d
