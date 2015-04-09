PROGRAM grid2d

  USE types
  USE vars, ONLY : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  USE grid, ONLY : grid_create,grid_destroy
  USE kernels, ONLY : kernel_select,kernel_support
  IMPLICIT NONE
  
  INTEGER(I4B), PARAMETER :: ndim = 3
  REAL(DP), DIMENSION(ndim), PARAMETER :: grav=(/ 0.0_dp, 0.0_dp, -980.0_dp /) 
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: periodic = (/ 0, 0, 0 /)
  REAL(DP), DIMENSION(ndim) :: xmin, xmax
  REAL(DP), PARAMETER :: rho0 = 1.D-3, gam = 7.0_dp
  REAL(DP), PARAMETER :: mu0 = 1.D-5
  CHARACTER(LEN=80) :: kernelname = 'm4',dir = '/home/tmattner/sph/',job = 'dam3d'
  CHARACTER(LEN=120) :: filnam
  INTEGER(I4B) :: nsteps=100
 
  INTEGER(I4B) :: i,j,k,np
  REAL(DP) :: h,dx,dt,p0,cs

  OPEN(11,FILE='dam_breaking.init')
  READ(11,*)
  READ(11,*)
  READ(11,*)np
  READ(11,*)dx
  CALL vars_create(ndim,np)
  DO i=1,np
     READ(11,'(I6,3(D16.8),1(I1))')j,x(1,i),x(3,i),x(2,i),bc(i)
     bc(i) = bc(i) - 1
     IF (bc(i) > 1) bc(i) = 1
  END DO
  CLOSE(11)
  
  h = 1.3_dp*dx
  DO i=1,ndim
     xmin(i) = MINVAL(x(i,:)) - 1.D-8
     xmax(i) = MAXVAL(x(i,:)) + 1.D-8
  END DO
  
  p0 = 200.0_dp/7.0_dp*rho0*ABS(grav(3))*55.0_dp
  cs = SQRT(gam*p0/rho0)
  dt = MIN(0.5_dp*dx/cs,0.5_dp*dx**2*rho0/mu0)

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

  u = 0.0_dp

  CALL vars_params(h,gam,rho0,p0,grav)
  CALL kernel_select(kernelname)
  CALL grid_create(xmin,xmax,kernel_support(h),periodic)
  
  m = rho0*dx**3
  rho = rho0
!  DO i=1,np
!     IF (bc(i) == 0) THEN
!        rho(i) = rho0*(1.0_dp + (gam - 1.0_dp)/gam*rho0*grav(3)*(x(3,i) - 55.0_dp)/p0)**(1.0_dp/(gam - 1.0_dp))
!     END IF
!  END DO
  CALL vars_pressure
  mu = mu0

  
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),0
  CALL vars_write(TRIM(filnam))
  CALL grid_destroy
  CALL vars_destroy
  
  
END PROGRAM grid2d
