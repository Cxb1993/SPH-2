PROGRAM weir3d

  USE types
  USE vars, ONLY : vars_create,vars_destroy,vars_write,vars_params,vars_density,vars_pressure,x,u,m,rho,p,mu,bc
  USE grid, ONLY : grid_create,grid_destroy
  USE kernels, ONLY : kernel_select,kernel_support
  IMPLICIT NONE
  
  INTEGER(I4B), PARAMETER :: ndim = 3
  REAL(DP), DIMENSION(ndim), PARAMETER :: grav=(/ 0.0_dp, -10.0_dp, 0.0_dp /) 
  INTEGER(I4B), DIMENSION(ndim), PARAMETER :: periodic = (/ 0, 0, 1 /)
  REAL(DP), PARAMETER :: dx = 0.1_dp
  REAL(DP), DIMENSION(ndim) :: xmin, xmax
  REAL(DP), PARAMETER :: rho0 = 1.0_dp, gam = 7.0_dp
  REAL(DP), PARAMETER :: mu0 = 1.D-8
  CHARACTER(LEN=80) :: kernelname = 'm4',dir = '/mnt/disk2/tmattner/SPH/',job = 'weir3d'
  CHARACTER(LEN=120) :: filnam
  INTEGER(I4B) :: nsteps=100
 
  INTEGER(I4B) :: i,j,k,np,cnt
  REAL(DP) :: h,dt,p0,cs

  cnt = 0

  np = 134100
  CALL vars_create(ndim,np)

  DO i=1,4
     DO j=1,100
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1)
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,50
     DO j=1,4
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*4
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO
     
  DO i=1,4
     DO j=1,25
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*54
            x(2,cnt) = dx*(j-1)
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,50
     DO j=1,4
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*58
            x(2,cnt) = dx*(j-1) + dx*21
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,4
     DO j=1,25
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*108
            x(2,cnt) = dx*(j-1) + dx*21
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,50
     DO j=1,4
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*112
            x(2,cnt) = dx*(j-1) + dx*42
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,4
     DO j=1,58
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*162
            x(2,cnt) = dx*(j-1) + dx*42
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 1
        END DO
     END DO
  END DO

  DO i=1,50
     DO j=1,25
        DO k=1,50
            cnt = cnt + 1
            x(1,cnt) = dx*(i-1) + dx*112
            x(2,cnt) = dx*(j-1) + dx*46
            x(3,cnt) = dx*(k-1)
            bc(cnt) = 0
        END DO
     END DO
  END DO

  np = cnt
  PRINT *,np

  h = 1.3_dp*dx
  DO i=1,2
     xmin(i) = MINVAL(x(i,:)) - 1.D-8
     xmax(i) = MAXVAL(x(i,:)) + 1.D-8
  END DO
  xmin(3) = MINVAL(x(3,:)) - 0.5_dp*dx ! for periodic bcs
  xmax(3) = MAXVAL(x(3,:)) + 0.5_dp*dx ! for periodic bcs
  
  p0 = 200.0_dp/7.0_dp*rho0*ABS(grav(2))*25.0_dp*dx
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
!  CALL vars_density
!  DO i=1,np
!     IF (bc(i) == 0) THEN
!        rho(i) = rho0*(1.0_dp + (gam - 1.0_dp)/gam*rho0*grav(2)*(x(2,i) - dx*70)/p0)**(1.0_dp/(gam - 1.0_dp))
!     END IF
!  END DO
  CALL vars_pressure
  mu = mu0

  
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),0
  CALL vars_write(TRIM(filnam))
  CALL grid_destroy
  CALL vars_destroy
  
  
END PROGRAM weir3d
