module kernels

  use types
  implicit none
  
  interface kernel_eval
    module procedure w
  end interface

  real(dp), dimension(3), parameter :: a0 = (/ 1.0_dp/sqrtpi_d, 1.0_dp/pi_d, 1.0_dp/pi_d/sqrtpi_d /)
  real(dp), dimension(3), parameter :: a4 = (/ 1.0_dp, 15.0_dp/7.0_dp/pi_d, 1.5_dp/pi_d /)
  real(dp), dimension(3), parameter :: a5 = (/ 1.0_dp/24.0_dp, 96.0_dp/1199.0_dp/pi_d, 0.05_dp/pi_d /) 
  real(dp), dimension(3), parameter :: a6 = (/ 1.0_dp/120.0_dp, 7.0_dp/478.0_dp/pi_d, 1.0_dp/120.0_dp/pi_d /)
  real(dp), parameter :: supp0 = 6.0_dp, supp4 = 2.0_dp, supp5 = 2.5_dp, supp6 = 3.0_dp

  character(len=80)      :: kernelname
  real(dp), dimension(3) :: alpha
  real(dp)               :: support

contains

  subroutine kernel_select(kernel)
   
    character(*) :: kernel

    kernelname = kernel
    select case(trim(kernelname))
    case('gaussian')
       alpha = a0
       support = supp0
    case('m4')
       alpha = a4
       support = supp4
    case('m5')
       alpha = a5
       support = supp5
    case('m6')
       alpha = a6
       support = supp6
    case default
       alpha = a4
       support = supp4
    end select

  end subroutine kernel_select
  
  function kernel_support(h)
  
    real(dp), intent(in) :: h
    real(dp) :: kernel_support
  
    kernel_support = support*h
    
  end function kernel_support

  function f(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: f

    select case(trim(kernelname))
    case('gaussian')
       f = m0(x,h)
    case('m4')
       f = m4(x,h)
    case('m5')
       f = m5(x,h)
    case('m6')
       f = m6(x,h)
    case default
       f = m4(x,h)
    end select  

  end function f

 function df(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: df

    select case(trim(kernelname))
    case('gaussian')
       df = dm0(x,h)
    case('m4')
       df = dm4(x,h)
    case('m5')
       df = dm5(x,h)
    case('m6')
       df = dm6(x,h)
    case default
       df = dm4(x,h)
    end select  

  end function df

 function d2f(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: d2f

    select case(trim(kernelname))
    case('gaussian')
       d2f = d2m0(x,h)
    case('m4')
       d2f = d2m4(x,h)
    case('m5')
       d2f = d2m5(x,h)
    case('m6')
       d2f = d2m6(x,h)
    case default
       d2f = d2m4(x,h)
    end select  

  end function d2f

  function w(x, h)

    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in)               :: h
    real(dp)                           :: w

    integer(i4b) :: ndim
    real(dp)     :: r

    ndim = size(x,1)
    r = sqrt(sum(x**2))
    w = f(r,h)*alpha(ndim)/h**ndim

  end function w

  function dw(i, x, h)

    integer(i4b), intent(in)           :: i
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in)               :: h
    real(dp)                           :: dw

    integer(i4b) :: ndim
    real(dp)     :: r

    ndim = size(x,1)
    r = sqrt(sum(x**2))
    if (r /= 0.0_dp) then
       dw = x(i)/r*df(r,h)*alpha(ndim)/h**ndim
    else
       dw = 0.0_dp
    end if

  end function dw
  
  function kernel_grad(x, h) result(grad)
  
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in)               :: h
    real(dp), dimension(size(x,1))     :: grad

    integer(i4b) :: ndim
    real(dp)     :: r

    ndim = size(x,1)
    r = sqrt(sum(x**2))
    if (r /= 0.0_dp) then
       grad = 1.0_dp/r*df(r,h)*alpha(ndim)/h**ndim*x
    else
       grad = 0.0_dp
    end if

  end function kernel_grad

  function d2w(i, j, x, h)

    integer(i4b), intent(in)           :: i, j
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in)               :: h
    real(dp)                           :: d2w

    real(dp), parameter :: tol = 1.d-14
    integer(i4b)        :: ndim
    real(dp)            :: r, tmp, a

    ndim = size(x,1)
    r = sqrt(sum(x**2))
    a = alpha(ndim)/h**ndim
    if (r >= tol) then
       tmp = a*df(r,h)/r
       if (i == j) then
          d2w = tmp
       else
          d2w = 0.0_dp
       end if
       d2w = d2w + x(i)*x(j)/r**2*(a*d2f(r,h) - tmp)
    else
       if (i == j) then
          d2w = a*d2f(r,h)
       else
          d2w = 0.0_dp
       end if
    end if

  end function d2w

  function m4(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: m4

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       m4 = ((2.0_dp - q)**3 - 4.0_dp*(1.0_dp - q)**3)/6.0_dp
    case(1)
       m4 = ((2.0_dp - q)**3)/6.0_dp
    case default
       m4 = 0.0_dp
    end select

  end function m4

  function dm4(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: dm4

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       dm4 = (-3.0_dp*(2.0_dp - q)**2 + 12.0_dp*(1.0_dp - q)**2)/6.0_dp
    case(1)
       dm4 = (-3.0_dp*(2.0_dp - q)**2)/6.0_dp
    case default
       dm4 = 0.0_dp
    end select
    dm4 = dm4*sign(1.0_dp/h, x)

  end function dm4

  function d2m4(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: d2m4

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       d2m4 = (6.0_dp*(2.0_dp - q) - 24.0_dp*(1.0_dp - q))/6.0_dp
    case(1)
       d2m4 = (6.0_dp*(2.0_dp - q))/6.0_dp
    case default
       d2m4 = 0.0_dp
    end select
    d2m4 = d2m4/h**2

  end function d2m4

  function m5(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: m5

    real(dp) :: q

    q = abs(x)/h
    select case(int(2.0_dp*q))
    case(0)
       m5 = ((2.5_dp - q)**4 - 5.0_dp*(1.5_dp - q)**4 &
            + 10.0_dp*(0.5_dp - q)**4)
    case(1:2)
       m5 = ((2.5_dp - q)**4 - 5.0_dp*(1.5_dp - q)**4)
    case(3:4)
       m5 = ((2.5_dp - q)**4)
    case default
       m5 = 0.0_dp
    end select

  end function m5

  function dm5(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: dm5

    real(dp) :: q

    q = abs(x)/h
    select case(int(2.0_dp*q))
    case(0)
       dm5 = (-4.0_dp*(2.5_dp - q)**3 + 20.0_dp*(1.5_dp - q)**3 &
            - 40.0_dp*(0.5_dp - q)**3)
    case(1:2)
       dm5 = (-4.0_dp*(2.5_dp - q)**3 + 20.0_dp*(1.5_dp - q)**3)
    case(3:4)
       dm5 = (-4.0_dp*(2.5_dp - q)**3)
    case default
       dm5 = 0.0_dp
    end select
    dm5 = dm5*sign(1.0_dp/h, x)

  end function dm5

  function d2m5(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: d2m5

    real(dp) :: q

    q = abs(x)/h
    select case(int(2.0_dp*q))
    case(0)
       d2m5 = (12.0_dp*(2.5_dp - q)**2 - 60.0_dp*(1.5_dp - q)**2 &
            + 120.0_dp*(0.5_dp - q)**2)
    case(1:2)
       d2m5 = (12.0_dp*(2.5_dp - q)**2 - 60.0_dp*(1.5_dp - q)**2)
    case(3:4)
       d2m5 = (12.0_dp*(2.5_dp - q)**2)
    case default
       d2m5 = 0.0_dp
    end select
    d2m5 = d2m5/h**2

  end function d2m5

  function m6(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: m6

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       m6 = (3.0_dp - q)**5 - 6.0_dp*(2.0_dp - q)**5 + 15.0_dp*(1.0_dp - q)**5
    case(1)
       m6 = (3.0_dp - q)**5 - 6.0_dp*(2.0_dp - q)**5
    case(2)
       m6 = (3.0_dp - q)**5
    case default
       m6 = 0.0_dp
    end select

  end function m6

  function dm6(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: dm6

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       dm6 = -5.0_dp*(3.0_dp - q)**4 + 30.0_dp*(2.0_dp - q)**4 - 75.0_dp*(1.0_dp - q)**4
    case(1)
       dm6 = -5.0_dp*(3.0_dp - q)**4 + 30.0_dp*(2.0_dp - q)**4
    case(2)
       dm6 = -5.0_dp*(3.0_dp - q)**4
    case default
       dm6 = 0.0_dp
    end select
    dm6 = dm6*sign(1.0_dp/h, x)

  end function dm6

  function d2m6(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: d2m6

    real(dp) :: q

    q = abs(x)/h
    select case(int(q))
    case(0)
       d2m6 = 20.0_dp*(3.0_dp - q)**3 - 120.0_dp*(2.0_dp - q)**3 + 300.0_dp*(1.0_dp - q)**3
    case(1)
       d2m6 = 20.0_dp*(3.0_dp - q)**3 - 120.0_dp*(2.0_dp - q)**3
    case(2)
       d2m6 = 20.0_dp*(3.0_dp - q)**3
    case default
       d2m6 = 0.0_dp
    end select
    d2m6 = d2m6/h**2

  end function d2m6

  function m0(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: m0

    real(dp) :: q

    q = abs(x)/h
    m0 = exp(-q**2)

  end function m0

  function dm0(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: dm0

    real(dp) :: q

    q = abs(x)/h
    dm0 = -2.0_dp*q*exp(-q**2)
    dm0 = dm0*sign(1.0_dp/h, x)

  end function dm0

  function d2m0(x, h)

    real(dp), intent(in) :: x, h
    real(dp)             :: d2m0

    real(dp) :: q

    q = abs(x)/h
    d2m0 = (4.0_dp*q**2 - 2.0_dp)*exp(-q**2)/h**2

  end function d2m0
    
end module kernels
