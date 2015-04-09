MODULE kernels

  USE types
  IMPLICIT NONE
  
  INTERFACE kernel_eval
    MODULE PROCEDURE w
  END INTERFACE

  REAL(DP), DIMENSION(3), PARAMETER :: a0 = (/ 1.0_dp/sqrtpi_d, 1.0_dp/pi_d, 1.0_dp/pi_d/sqrtpi_d /)
  REAL(DP), DIMENSION(3), PARAMETER :: a4 = (/ 1.0_dp, 15.0_dp/7.0_dp/pi_d, 1.5_dp/pi_d /)
  REAL(DP), DIMENSION(3), PARAMETER :: a5 = (/ 1.0_dp/24.0_dp, 96.0_dp/1199.0_dp/pi_d, 0.05_dp/pi_d /) 
  REAL(DP), DIMENSION(3), PARAMETER :: a6 = (/ 1.0_dp/120.0_dp, 7.0_dp/478.0_dp/pi_d, 1.0_dp/120.0_dp/pi_d /)
  REAL(DP), PARAMETER :: supp0 = 6.0_dp, supp4 = 2.0_dp, supp5 = 2.5_dp, supp6 = 3.0_dp

  CHARACTER(LEN=80) :: kernelname
  REAL(DP), DIMENSION(3) :: alpha
  REAL(DP) :: support

CONTAINS

  SUBROUTINE kernel_select(kernel)
   
    CHARACTER(*) :: kernel

    kernelname = kernel
    SELECT CASE(TRIM(kernelname))
    CASE('gaussian')
       alpha = a0
       support = supp0
    CASE('m4')
       alpha = a4
       support = supp4
    CASE('m5')
       alpha = a5
       support = supp5
    CASE('m6')
       alpha = a6
       support = supp6
    CASE DEFAULT
       alpha = a4
       support = supp4
    END SELECT

  END SUBROUTINE kernel_select
  
  FUNCTION kernel_support(h)
  
    REAL(DP), INTENT(IN) :: h
    REAL(DP) :: kernel_support
  
    kernel_support = support*h
    
  END FUNCTION kernel_support

  FUNCTION f(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: f

    SELECT CASE(TRIM(kernelname))
    CASE('gaussian')
       f = m0(x,h)
    CASE('m4')
       f = m4(x,h)
    CASE('m5')
       f = m5(x,h)
    CASE('m6')
       f = m6(x,h)
    CASE DEFAULT
       f = m4(x,h)
    END SELECT  

  END FUNCTION f

 FUNCTION df(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: df

    SELECT CASE(TRIM(kernelname))
    CASE('gaussian')
       df = dm0(x,h)
    CASE('m4')
       df = dm4(x,h)
    CASE('m5')
       df = dm5(x,h)
    CASE('m6')
       df = dm6(x,h)
    CASE DEFAULT
       df = dm4(x,h)
    END SELECT  

  END FUNCTION df

 FUNCTION d2f(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: d2f

    SELECT CASE(TRIM(kernelname))
    CASE('gaussian')
       d2f = d2m0(x,h)
    CASE('m4')
       d2f = d2m4(x,h)
    CASE('m5')
       d2f = d2m5(x,h)
    CASE('m6')
       d2f = d2m6(x,h)
    CASE DEFAULT
       d2f = d2m4(x,h)
    END SELECT  

  END FUNCTION d2f

  FUNCTION w(x,h)

    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: h
    REAL(DP) :: w

    INTEGER(I4B) :: ndim
    REAL(DP) :: r

    ndim = SIZE(x,1)
    r = SQRT(SUM(x**2))
    w = f(r,h)*alpha(ndim)/h**ndim

  END FUNCTION w

  FUNCTION dw(i,x,h)

    INTEGER(I4B), INTENT(IN) :: i
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: h
    REAL(DP) :: dw

    INTEGER(I4B) :: ndim
    REAL(DP) :: r

    ndim = SIZE(x,1)
    r = SQRT(SUM(x**2))
    IF (r /= 0.0_dp) THEN
       dw = x(i)/r*df(r,h)*alpha(ndim)/h**ndim
    ELSE
       dw = 0.0_dp
    END IF

  END FUNCTION dw
  
  FUNCTION kernel_grad(x,h) RESULT(grad)
  
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: h
    REAL(DP), DIMENSION(SIZE(x,1)) :: grad

    INTEGER(I4B) :: ndim
    REAL(DP) :: r

    ndim = SIZE(x,1)
    r = SQRT(SUM(x**2))
    IF (r /= 0.0_dp) THEN
       grad = 1.0_dp/r*df(r,h)*alpha(ndim)/h**ndim*x
    ELSE
       grad = 0.0_dp
    END IF

  END FUNCTION kernel_grad

  FUNCTION d2w(i,j,x,h)

    INTEGER(I4B), INTENT(IN) :: i,j
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: h
    REAL(DP) :: d2w

    REAL(DP), PARAMETER :: tol = 1.D-14
    INTEGER(I4B) :: ndim
    REAL(DP) :: r,tmp,a

    ndim = SIZE(x,1)
    r = SQRT(SUM(x**2))
    a = alpha(ndim)/h**ndim
    IF (r >= tol) THEN
       tmp = a*df(r,h)/r
       IF (i == j) THEN
          d2w = tmp
       ELSE
          d2w = 0.0_dp
       END IF
       d2w = d2w + x(i)*x(j)/r**2*(a*d2f(r,h) - tmp)
    ELSE
       IF (i == j) THEN
          d2w = a*d2f(r,h)
       ELSE
          d2w = 0.0_dp
       END IF
    END IF

  END FUNCTION d2w

  FUNCTION m4(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: m4

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       m4 = ((2.0_dp - q)**3 - 4.0_dp*(1.0_dp - q)**3)/6.0_dp
    CASE(1)
       m4 = ((2.0_dp - q)**3)/6.0_dp
    CASE DEFAULT
       m4 = 0.0_dp
    END SELECT

  END FUNCTION m4

  FUNCTION dm4(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: dm4

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       dm4 = (-3.0_dp*(2.0_dp - q)**2 + 12.0_dp*(1.0_dp - q)**2)/6.0_dp
    CASE(1)
       dm4 = (-3.0_dp*(2.0_dp - q)**2)/6.0_dp
    CASE DEFAULT
       dm4 = 0.0_dp
    END SELECT
    dm4 = dm4*SIGN(1.0_dp/h,x)

  END FUNCTION dm4

  FUNCTION d2m4(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: d2m4

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       d2m4 = (6.0_dp*(2.0_dp - q) - 24.0_dp*(1.0_dp - q))/6.0_dp
    CASE(1)
       d2m4 = (6.0_dp*(2.0_dp - q))/6.0_dp
    CASE DEFAULT
       d2m4 = 0.0_dp
    END SELECT
    d2m4 = d2m4/h**2

  END FUNCTION d2m4

  FUNCTION m5(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: m5

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(2.0_dp*q))
    CASE(0)
       m5 = ((2.5_dp - q)**4 - 5.0_dp*(1.5_dp - q)**4 &
            + 10.0_dp*(0.5_dp - q)**4)
    CASE(1:2)
       m5 = ((2.5_dp - q)**4 - 5.0_dp*(1.5_dp - q)**4)
    CASE(3:4)
       m5 = ((2.5_dp - q)**4)
    CASE DEFAULT
       m5 = 0.0_dp
    END SELECT

  END FUNCTION m5

  FUNCTION dm5(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: dm5

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(2.0_dp*q))
    CASE(0)
       dm5 = (-4.0_dp*(2.5_dp - q)**3 + 20.0_dp*(1.5_dp - q)**3 &
            - 40.0_dp*(0.5_dp - q)**3)
    CASE(1:2)
       dm5 = (-4.0_dp*(2.5_dp - q)**3 + 20.0_dp*(1.5_dp - q)**3)
    CASE(3:4)
       dm5 = (-4.0_dp*(2.5_dp - q)**3)
    CASE DEFAULT
       dm5 = 0.0_dp
    END SELECT
    dm5 = dm5*SIGN(1.0_dp/h,x)

  END FUNCTION dm5

  FUNCTION d2m5(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: d2m5

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(2.0_dp*q))
    CASE(0)
       d2m5 = (12.0_dp*(2.5_dp - q)**2 - 60.0_dp*(1.5_dp - q)**2 &
            + 120.0_dp*(0.5_dp - q)**2)
    CASE(1:2)
       d2m5 = (12.0_dp*(2.5_dp - q)**2 - 60.0_dp*(1.5_dp - q)**2)
    CASE(3:4)
       d2m5 = (12.0_dp*(2.5_dp - q)**2)
    CASE DEFAULT
       d2m5 = 0.0_dp
    END SELECT
    d2m5 = d2m5/h**2

  END FUNCTION d2m5

  FUNCTION m6(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: m6

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       m6 = (3.0_dp - q)**5 - 6.0_dp*(2.0_dp - q)**5 + 15.0_dp*(1.0_dp - q)**5
    CASE(1)
       m6 = (3.0_dp - q)**5 - 6.0_dp*(2.0_dp - q)**5
    CASE(2)
       m6 = (3.0_dp - q)**5
    CASE DEFAULT
       m6 = 0.0_dp
    END SELECT

  END FUNCTION m6

  FUNCTION dm6(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: dm6

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       dm6 = -5.0_dp*(3.0_dp - q)**4 + 30.0_dp*(2.0_dp - q)**4 - 75.0_dp*(1.0_dp - q)**4
    CASE(1)
       dm6 = -5.0_dp*(3.0_dp - q)**4 + 30.0_dp*(2.0_dp - q)**4
    CASE(2)
       dm6 = -5.0_dp*(3.0_dp - q)**4
    CASE DEFAULT
       dm6 = 0.0_dp
    END SELECT
    dm6 = dm6*SIGN(1.0_dp/h,x)

  END FUNCTION dm6

  FUNCTION d2m6(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: d2m6

    REAL(DP) :: q

    q = ABS(x)/h
    SELECT CASE(INT(q))
    CASE(0)
       d2m6 = 20.0_dp*(3.0_dp - q)**3 - 120.0_dp*(2.0_dp - q)**3 + 300.0_dp*(1.0_dp - q)**3
    CASE(1)
       d2m6 = 20.0_dp*(3.0_dp - q)**3 - 120.0_dp*(2.0_dp - q)**3
    CASE(2)
       d2m6 = 20.0_dp*(3.0_dp - q)**3
    CASE DEFAULT
       d2m6 = 0.0_dp
    END SELECT
    d2m6 = d2m6/h**2

  END FUNCTION d2m6

  FUNCTION m0(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: m0

    REAL(DP) :: q

    q = ABS(x)/h
    m0 = EXP(-q**2)

  END FUNCTION m0

  FUNCTION dm0(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: dm0

    REAL(DP) :: q

    q = ABS(x)/h
    dm0 = -2.0_dp*q*EXP(-q**2)
    dm0 = dm0*SIGN(1.0_dp/h,x)

  END FUNCTION dm0

  FUNCTION d2m0(x,h)

    REAL(DP), INTENT(IN) :: x,h
    REAL(DP) :: d2m0

    REAL(DP) :: q

    q = ABS(x)/h
    d2m0 = (4.0_dp*q**2 - 2.0_dp)*EXP(-q**2)/h**2

  END FUNCTION d2m0
    
END MODULE kernels
