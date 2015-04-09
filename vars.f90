MODULE vars

  USE types
  IMPLICIT NONE

  INTERFACE vars_step
    MODULE PROCEDURE vars_step_rk3
  END INTERFACE
  
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: x,u,du
  REAL(DP), DIMENSION(:), ALLOCATABLE :: m,rho,drho,p,mu
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: bc
  INTEGER(I4B) :: ndim,np
  REAL(DP) :: h,gam,rho0,p0,simtime
  REAL(DP), DIMENSION(:), ALLOCATABLE :: grav
  
CONTAINS

  SUBROUTINE vars_create(nd,n)
  
    INTEGER(I4B), INTENT(IN) :: nd,n
    
    ndim = nd
    np = n
    ALLOCATE(x(ndim,np),u(ndim,np),du(ndim,np),m(np),rho(np),drho(np),p(np),mu(np),bc(np),grav(ndim))
    
  END SUBROUTINE vars_create
  
  SUBROUTINE vars_destroy()
  
    DEALLOCATE(x,u,du,m,rho,drho,p,mu,bc,grav)
    
  END SUBROUTINE vars_destroy
  
  SUBROUTINE vars_params(hdum,gamdum,rho0dum,p0dum,gravdum)
  
    REAL(DP), INTENT(IN) :: hdum,gamdum,rho0dum,p0dum
    REAL(DP), DIMENSION(ndim), INTENT(IN) :: gravdum
    
    h = hdum
    gam = gamdum
    rho0 = rho0dum
    p0 = p0dum
    grav = gravdum
    
  END SUBROUTINE vars_params
  
  SUBROUTINE vars_read(filename)
  
    CHARACTER(*), INTENT(IN) :: filename
    
    INTEGER(I4B) :: i
    
    OPEN(1,FILE=TRIM(filename))
    DO i=1,np
       READ(1,*)x(:,i),u(:,i),m(i),rho(i),p(i),mu(i),bc(i)
    END DO
    CLOSE(1)
    
  END SUBROUTINE vars_read
  
  SUBROUTINE vars_write(filename)
  
    CHARACTER(*), INTENT(IN) :: filename
    
    INTEGER(I4B) :: i
    
    OPEN(1,FILE=TRIM(filename))
    DO i=1,np
       WRITE(1,*)x(:,i),u(:,i),m(i),rho(i),p(i),mu(i),bc(i)
    END DO
    CLOSE(1)
    
  END SUBROUTINE vars_write
  
  SUBROUTINE vars_rates(a,b)
  
    USE kernels, ONLY : kernel_support,kernel_grad!,w
    USE grid, ONLY : xdiff
  
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    
    REAL(DP), DIMENSION(3), PARAMETER :: visc = (/ 12.0_dp, 16.0_dp, 20.0_dp /) ! Don't know for 1d
    REAL(DP), PARAMETER :: eps = 1.D-8
    INTEGER(I4B) :: na,nb,i,j
    REAL(DP), DIMENSION(ndim) :: rab,uab,wab,piab
    REAL(DP) :: supp,twosupp,supp2,rab2,pab!,fab,w0
    
    na = SIZE(a,1)
    nb = SIZE(b,1)

!    rab = 0.0_dp
!    rab(1) = h
!    w0 = w(rab,h)
    
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
    !$OMP PARALLEL SHARED(DU, DRHO) PRIVATE(I, J, RAB, RAB2, UAB, WAB, PIAB, PAB)
    !$OMP DO 
    DO i=1,na
       SELECT CASE(bc(a(i)))
       CASE (1)
          du(:,a(i)) = 0.0_dp
          drho(a(i)) = 0.0_dp
          DO j=1,nb
             rab = x(:,a(i)) - x(:,b(j))
             WHERE (rab > twosupp)
                rab = rab - xdiff
             ELSEWHERE (rab < -twosupp)
                rab = rab + xdiff
             END WHERE   
             rab2 = SUM(rab**2)
             IF (rab2 < supp2) THEN
                uab = u(:,a(i)) - u(:,b(j)) 
                wab = kernel_grad(rab,h)
                drho(a(i)) = drho(a(i)) + m(b(j))/rho(b(j))*DOT_PRODUCT(uab,wab)
             END IF
          END DO
          drho(a(i)) = drho(a(i))*rho(a(i))
       CASE DEFAULT
          du(:,a(i)) = 0.0_dp
          drho(a(i)) = 0.0_dp
          DO j=1,nb
             rab = x(:,a(i)) - x(:,b(j))
             WHERE (rab > twosupp)
                rab = rab - xdiff
             ELSEWHERE (rab < -twosupp)
                rab = rab + xdiff
             END WHERE   
             rab2 = SUM(rab**2)
             IF (rab2 < supp2) THEN
                uab = u(:,a(i)) - u(:,b(j)) 
                wab = kernel_grad(rab,h)
                piab = -visc(ndim)*mu(a(i))*mu(b(j))/rho(a(i))/rho(b(j))/(mu(a(i)) + mu(b(j))) &
                       /(rab2 + eps*h**2)*DOT_PRODUCT(uab,rab)
                pab = p(b(j))/rho(b(j))**2 + p(a(i))/rho(a(i))**2
!                fab = 0.01_dp*ABS(pab)*(w(rab,h)/w0)**4 ! anticlumping
                drho(a(i)) = drho(a(i)) + m(b(j))/rho(b(j))*DOT_PRODUCT(uab,wab)
!                du(:,a(i)) = du(:,a(i)) - m(b(j))*(pab + piab + fab)*wab
                du(:,a(i)) = du(:,a(i)) - m(b(j))*(pab + piab)*wab
             END IF
          END DO
          drho(a(i)) = drho(a(i))*rho(a(i))
          du(:,a(i)) = du(:,a(i)) + grav
       END SELECT
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  
  END SUBROUTINE vars_rates

  SUBROUTINE vars_rho(a,b)
  
    USE kernels, ONLY : kernel_support,kernel_eval
    USE grid, ONLY : xdiff
  
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    
    INTEGER(I4B) :: na,nb,i,j
    REAL(DP), DIMENSION(ndim) :: rab
    REAL(DP) :: supp,twosupp,supp2,rab2
    
    na = SIZE(a,1)
    nb = SIZE(b,1)
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
    DO i=1,na
       rho(a(i)) = 0.0_dp
       DO j=1,nb
          rab = x(:,a(i)) - x(:,b(j))
! apply periodic correction if any component > 2*support
          WHERE (rab > twosupp)
             rab = rab - xdiff
          ELSEWHERE (rab < -twosupp)
             rab = rab + xdiff
          END WHERE   
! some physics
          rab2 = SUM(rab**2)
          IF (rab2 < supp2) THEN
             rho(a(i)) = rho(a(i)) + m(b(j))*kernel_eval(rab,h)
          END IF
       END DO
    END DO
  
  END SUBROUTINE vars_rho
  
  SUBROUTINE vars_pressure()
  
    p = p0*((rho/rho0)**gam - 1.0_dp)
  
  END SUBROUTINE vars_pressure

  SUBROUTINE vars_density()
  
    USE grid, ONLY : grid_clear,grid_insert,grid_loop,grid_periodic_map

!    INTEGER(I4B), DIMENSION(np) :: ii
!    INTEGER(I4B) :: i
!    FORALL (i=1:np) ii(i) = i
  
    CALL grid_clear
    CALL grid_insert(x)
    CALL grid_loop(vars_rho)

!    CALL vars_rho(ii,ii)
  END SUBROUTINE vars_density
  
  SUBROUTINE vars_step_euler(dt)
  
    USE grid, ONLY : grid_clear,grid_insert,grid_loop,grid_periodic_map
  
    REAL(DP), INTENT(IN) :: dt
    
    INTEGER(I4B) :: i

!    INTEGER(I4B), DIMENSION(np) :: ii
!    INTEGER(I4B) :: i
!    FORALL (i=1:np) ii(i) = i
    
    CALL vars_pressure
    CALL grid_clear
    CALL grid_insert(x)
    CALL grid_loop(vars_rates)
!    CALL vars_rates(ii,ii)
    
    rho = rho + dt*drho
    x = x + dt*u
    u = u + dt*du

    PRINT '(3(E16.10,X))',(SUM(m*u(i,:)),i=1,ndim)
    
    CALL grid_periodic_map(x)
    
  END SUBROUTINE vars_step_euler

 SUBROUTINE vars_step_rk3(dt)
  
    USE grid, ONLY : grid_clear,grid_insert,grid_loop,grid_periodic_map
  
    REAL(DP), INTENT(IN) :: dt

    INTEGER(I4B), PARAMETER :: nstep = 3
    REAL(DP), DIMENSION(nstep) :: w = (/ 1.0_dp/6.0_dp, 3.0_dp/10.0_dp,  8.0_dp/15.0_dp /)
    REAL(DP), DIMENSION(nstep) :: a = (/ 0.0_dp, -5.0_dp/9.0_dp,  -153.0_dp/128.0_dp  /)
    REAL(DP), DIMENSION(nstep) :: b = (/ 1.0_dp/3.0_dp, 15.0_dp/16.0_dp, 8.0_dp/15.0_dp /)
    REAL(DP), DIMENSION(nstep) :: c = (/ 1.0_dp/3.0_dp, 5.0_dp/12.0_dp,  1.0_dp/4.0_dp  /)

    REAL(DP), DIMENSION(ndim,np) :: us,xs,xx
    REAL(DP), DIMENSION(np) :: rhos
    INTEGER(I4B) :: i
    
    xx = x
    rhos = 0.0_dp
    us = 0.0_dp
    xs = 0.0_dp
    DO i=1,nstep
       CALL vars_pressure
       CALL grid_clear
       x = xx ! takes into account boundary crossings
       CALL grid_periodic_map(x)
       CALL grid_insert(x)
       CALL grid_loop(vars_rates)

       rhos = a(i)*rhos + dt*drho
       rho  = b(i)*rhos + rho
       xs = a(i)*xs + dt*u
       xx = b(i)*xs + xx
       us = a(i)*us + dt*du
       u  = b(i)*us + u
       simtime = simtime + c(i)*dt
    END DO
    x = xx
    CALL grid_periodic_map(x)
    PRINT '(3(E16.10,X))',(SUM(m*u(i,:)),i=1,ndim)
    
  END SUBROUTINE vars_step_rk3
  


  
  
END MODULE vars
