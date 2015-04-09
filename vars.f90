module vars

  use types
  implicit none

  interface vars_step
    module procedure vars_step_rk3
  end interface
  
  real(dp), dimension(:,:), allocatable   :: x, u, du
  real(dp), dimension(:), allocatable     :: m, rho, drho, p, mu
  integer(i4b), dimension(:), allocatable :: bc
  integer(i4b)                            :: ndim, np
  real(dp)                                :: h, gam, rho0, p0, simtime
  real(dp), dimension(:), allocatable     :: grav
  
contains

  subroutine vars_create(nd, n)
  
    integer(i4b), intent(in) :: nd, n
    
    ndim = nd
    np = n
    allocate(x(ndim,np), u(ndim,np), du(ndim,np))
    allocate(m(np), rho(np), drho(np), p(np))
    allocate(mu(np), bc(np), grav(ndim))
    
  end subroutine vars_create
  
  subroutine vars_destroy()
  
    deallocate(x, u, du, m, rho, drho, p, mu, bc, grav)
    
  end subroutine vars_destroy
  
  subroutine vars_params(hdum, gamdum, rho0dum, p0dum, gravdum)
  
    real(dp), intent(in)                  :: hdum, gamdum, rho0dum, p0dum
    real(dp), dimension(ndim), intent(in) :: gravdum
    
    h = hdum
    gam = gamdum
    rho0 = rho0dum
    p0 = p0dum
    grav = gravdum
    
  end subroutine vars_params
  
  subroutine vars_read(filename)
  
    character(*), intent(in) :: filename
    
    integer(i4b) :: i
    
    open(1, file=trim(filename))
    do i=1,np
       read(1,*) x(:,i), u(:,i), m(i), rho(i), p(i), mu(i), bc(i)
    end do
    close(1)
    
  end subroutine vars_read
  
  subroutine vars_write(filename)
  
    character(*), intent(in) :: filename
    
    integer(i4b) :: i
    
    open(1, file=trim(filename))
    do i=1,np
       write(1,*) x(:,i), u(:,i), m(i), rho(i), p(i), mu(i), bc(i)
    end do
    close(1)
    
  end subroutine vars_write
  
  subroutine vars_rates(a, b)
  
    use kernels, only : kernel_support, kernel_grad!,w
    use grid, only : xdiff
  
    integer(i4b), dimension(:), intent(in) :: a, b
    
    real(dp), dimension(3), parameter :: visc = (/ 12.0_dp, 16.0_dp, 20.0_dp /) ! don't know for 1d
    real(dp), parameter               :: eps = 1.d-8
    integer(i4b)                      :: na, nb, i, j
    real(dp), dimension(ndim)         :: rab, uab, wab, piab
    real(dp)                          :: supp, twosupp, supp2, rab2, pab!,fab,w0
    
    na = size(a,1)
    nb = size(b,1)

!    rab = 0.0_dp
!    rab(1) = h
!    w0 = w(rab,h)
    
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
 !   !$omp parallel default(private) shared(du, drho)
 !   !$omp do
    do i=1,na
       select case(bc(a(i)))
       case (1)
          du(:,a(i)) = 0.0_dp
          drho(a(i)) = 0.0_dp
          do j=1,nb
             rab = x(:,a(i)) - x(:,b(j))
             where (rab > twosupp)
                rab = rab - xdiff
             elsewhere (rab < -twosupp)
                rab = rab + xdiff
             end where   
             rab2 = sum(rab**2)
             if (rab2 < supp2) then
                uab = u(:,a(i)) - u(:,b(j)) 
                wab = kernel_grad(rab, h)
                drho(a(i)) = drho(a(i)) + m(b(j))/rho(b(j))*dot_product(uab,wab)
             end if
          end do
          drho(a(i)) = drho(a(i))*rho(a(i))
       case default
          du(:,a(i)) = 0.0_dp
          drho(a(i)) = 0.0_dp
          do j=1,nb
             rab = x(:,a(i)) - x(:,b(j))
             where (rab > twosupp)
                rab = rab - xdiff
             elsewhere (rab < -twosupp)
                rab = rab + xdiff
             end where   
             rab2 = sum(rab**2)
             if (rab2 < supp2) then
                uab = u(:,a(i)) - u(:,b(j)) 
                wab = kernel_grad(rab, h)
                piab = -visc(ndim)*mu(a(i))*mu(b(j))/rho(a(i))/rho(b(j))/(mu(a(i)) + mu(b(j))) &
                       /(rab2 + eps*h**2)*dot_product(uab,rab)
                pab = p(b(j))/rho(b(j))**2 + p(a(i))/rho(a(i))**2
!                fab = 0.01_dp*abs(pab)*(w(rab,h)/w0)**4 ! anticlumping
                drho(a(i)) = drho(a(i)) + m(b(j))/rho(b(j))*dot_product(uab,wab)
!                du(:,a(i)) = du(:,a(i)) - m(b(j))*(pab + piab + fab)*wab
                du(:,a(i)) = du(:,a(i)) - m(b(j))*(pab + piab)*wab
             end if
          end do
          drho(a(i)) = drho(a(i))*rho(a(i))
          du(:,a(i)) = du(:,a(i)) + grav
       end select
    end do
 !   !$omp end do
 !   !$omp end parallel
  
  end subroutine vars_rates

  subroutine vars_rho(a, b)
  
    use kernels, only : kernel_support, kernel_eval
    use grid, only : xdiff
  
    integer(i4b), dimension(:), intent(in) :: a, b
    
    integer(i4b)              :: na, nb, i, j
    real(dp), dimension(ndim) :: rab
    real(dp)                  :: supp, twosupp, supp2, rab2
    
    na = size(a,1)
    nb = size(b,1)
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
    do i=1,na
       rho(a(i)) = 0.0_dp
       do j=1,nb
          rab = x(:,a(i)) - x(:,b(j))
! apply periodic correction if any component > 2*support
          where (rab > twosupp)
             rab = rab - xdiff
          elsewhere (rab < -twosupp)
             rab = rab + xdiff
          end where   
! some physics
          rab2 = sum(rab**2)
          if (rab2 < supp2) then
             rho(a(i)) = rho(a(i)) + m(b(j))*kernel_eval(rab, h)
          end if
       end do
    end do
  
  end subroutine vars_rho
  
  subroutine vars_pressure()
  
    p = p0*((rho/rho0)**gam - 1.0_dp)
  
  end subroutine vars_pressure

  subroutine vars_density()
  
    use grid, only : grid_clear, grid_insert, grid_loop, grid_periodic_map

!    integer(i4b), dimension(np) :: ii
!    integer(i4b) :: i
!    forall (i=1:np) ii(i) = i
  
    call grid_clear
    call grid_insert(x)
    call grid_loop(vars_rho)

!    call vars_rho(ii,ii)
  end subroutine vars_density
  
  subroutine vars_step_euler(dt)
  
    use grid, only : grid_clear, grid_insert, grid_loop, grid_periodic_map
  
    real(dp), intent(in) :: dt
    
    integer(i4b) :: i

!    integer(i4b), dimension(np) :: ii
!    integer(i4b) :: i
!    forall (i=1:np) ii(i) = i
    
    call vars_pressure
    call grid_clear
    call grid_insert(x)
    call grid_loop(vars_rates)
!    call vars_rates(ii,ii)
    
    rho = rho + dt*drho
    x = x + dt*u
    u = u + dt*du

    print '(3(e16.10,x))', (sum(m*u(i,:)), i=1,ndim)
    
    call grid_periodic_map(x)
    
  end subroutine vars_step_euler

 subroutine vars_step_rk3(dt)
  
    use grid, only : grid_clear, grid_insert, grid_loop, grid_periodic_map
  
    real(dp), intent(in) :: dt

    integer(i4b), parameter    :: nstep = 3
    real(dp), dimension(nstep) :: w = (/ 1.0_dp/6.0_dp, 3.0_dp/10.0_dp,  8.0_dp/15.0_dp /)
    real(dp), dimension(nstep) :: a = (/ 0.0_dp, -5.0_dp/9.0_dp,  -153.0_dp/128.0_dp  /)
    real(dp), dimension(nstep) :: b = (/ 1.0_dp/3.0_dp, 15.0_dp/16.0_dp, 8.0_dp/15.0_dp /)
    real(dp), dimension(nstep) :: c = (/ 1.0_dp/3.0_dp, 5.0_dp/12.0_dp,  1.0_dp/4.0_dp  /)

    real(dp), dimension(ndim,np) :: us, xs, xx
    real(dp), dimension(np)      :: rhos
    integer(i4b)                 :: i
    
    xx = x
    rhos = 0.0_dp
    us = 0.0_dp
    xs = 0.0_dp
    do i=1,nstep
       call vars_pressure
       call grid_clear
       x = xx ! takes into account boundary crossings
       call grid_periodic_map(x)
       call grid_insert(x)
       call grid_loop(vars_rates)

       rhos = a(i)*rhos + dt*drho
       rho  = b(i)*rhos + rho
       xs = a(i)*xs + dt*u
       xx = b(i)*xs + xx
       us = a(i)*us + dt*du
       u  = b(i)*us + u
       simtime = simtime + c(i)*dt
    end do
    x = xx
    call grid_periodic_map(x)
    print '(3(e16.10,x))', (sum(m*u(i,:)), i=1,ndim)
    
  end subroutine vars_step_rk3
  
end module vars
