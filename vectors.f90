MODULE vectors

  USE types
  IMPLICIT NONE

  INTERFACE vector_assign
    MODULE PROCEDURE vectori4b_assign,vectori4b_assign_vector
  END INTERFACE
  
  INTERFACE vector_cat
    MODULE PROCEDURE vectori4b_cat_vector,vectori4b_cat_scalar
  END INTERFACE 
  
  TYPE vectori4b
     INTEGER(I4B) :: ndat,nelem
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: dat
  END TYPE vectori4b
  
CONTAINS

  SUBROUTINE vectori4b_create(n,x)

    INTEGER(I4B), INTENT(IN) :: n
    TYPE(vectori4b), INTENT(INOUT) :: x
    
    x%ndat = 0
    x%nelem = n
    ALLOCATE(x%dat(n))
    
  END SUBROUTINE vectori4b_create
  
  SUBROUTINE vectori4b_destroy(x)

    TYPE (vectori4b), INTENT(INOUT) :: x
    
    DEALLOCATE(x%dat)
    
  END SUBROUTINE vectori4b_destroy
  
  SUBROUTINE vectori4b_double(x)
  
    TYPE(vectori4b), INTENT(INOUT) :: x
    
    TYPE(vectori4b) :: tmp
    INTEGER(I4B) :: nelem,ndat
    
    nelem = x%nelem
    ndat = x%ndat
    CALL vectori4b_create(nelem,tmp)
    tmp%ndat = ndat
    tmp%dat(1:ndat) = x%dat(1:ndat)
    CALL vectori4b_destroy(x)
    CALL vectori4b_create(2*nelem,x)
    x%ndat = ndat
    x%dat(1:ndat) = tmp%dat(1:ndat)
    CALL vectori4b_destroy(tmp)
  
  END SUBROUTINE vectori4b_double
  
  SUBROUTINE vectori4b_cat_scalar(x,a)
  
    TYPE(vectori4b), INTENT(INOUT) :: x 
    INTEGER(I4B), INTENT(IN) :: a
    
    INTEGER(I4B) :: ndat
    
    ndat = x%ndat + 1
    IF (ndat > x%nelem) THEN
       CALL vectori4b_double(x)
    END IF
    x%ndat = ndat
    x%dat(ndat) = a
    
  END SUBROUTINE vectori4b_cat_scalar
  
  SUBROUTINE vectori4b_cat_vector(x,a)
  
    TYPE(vectori4b), INTENT(INOUT) :: x 
    TYPE(vectori4b), INTENT(IN) :: a
    
    INTEGER(I4B) :: ndat
    
    ndat = x%ndat + a%ndat
    DO WHILE (ndat > x%nelem)
       CALL vectori4b_double(x)
    END DO
    x%dat(x%ndat+1:ndat) = a%dat(1:a%ndat)
    x%ndat = ndat
    
  END SUBROUTINE vectori4b_cat_vector

  SUBROUTINE vectori4b_assign(x,a)
  
    TYPE(vectori4b), INTENT(INOUT) :: x 
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a
    
    INTEGER(I4B) :: ndat
    
    ndat = SIZE(a,1)
    DO WHILE (ndat > x%nelem)
       CALL vectori4b_double(x)
    END DO
    x%dat(1:ndat) = a
    x%ndat = ndat
    
  END SUBROUTINE vectori4b_assign

  SUBROUTINE vectori4b_assign_vector(x,a)
  
    TYPE(vectori4b), INTENT(INOUT) :: x 
    TYPE(vectori4b), INTENT(IN) :: a
    
    INTEGER(I4B) :: ndat
    
    IF (a%nelem /= x%nelem) THEN
       DEALLOCATE(x%dat)
       ALLOCATE(x%dat(a%nelem))
       x%nelem = a%nelem
    END IF
    x%ndat = a%ndat
    x%dat = a%dat
    
  END SUBROUTINE vectori4b_assign_vector
  
  SUBROUTINE vectori4b_clear(x)
  
    TYPE(vectori4b), INTENT(INOUT) :: x
    
    x%ndat = 0
    
  END SUBROUTINE vectori4b_clear
  
  FUNCTION vectori4b_size(x) RESULT(ndat)
  
    TYPE(vectori4b), INTENT(INOUT) :: x 
    INTEGER(I4B) :: ndat
    
    ndat = x%ndat
    
  END FUNCTION vectori4b_size
    

END MODULE vectors
