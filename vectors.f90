module vectors

  use types
  implicit none

  interface vector_assign
    module procedure vectori4b_assign, vectori4b_assign_vector
  end interface
  
  interface vector_cat
    module procedure vectori4b_cat_vector, vectori4b_cat_scalar
  end interface 
  
  type vectori4b
     integer(i4b) :: ndat, nelem
  !   integer(i4b), dimension(:), allocatable :: dat
     integer(i4b), dimension(:), pointer :: dat
  end type vectori4b
  
contains

  subroutine vectori4b_create(n, x)

    integer(i4b), intent(in)       :: n
    type(vectori4b), intent(inout) :: x
    
    x%ndat = 0
    x%nelem = n
    allocate(x%dat(n))
    
  end subroutine vectori4b_create
  
  subroutine vectori4b_destroy(x)

    type (vectori4b), intent(inout) :: x
    
    deallocate(x%dat)
    
  end subroutine vectori4b_destroy
  
  subroutine vectori4b_double(x)
  
    type(vectori4b), intent(inout) :: x
    
    type(vectori4b) :: tmp
    integer(i4b)    :: nelem, ndat
    
    nelem = x%nelem
    ndat = x%ndat
    call vectori4b_create(nelem,tmp)
    tmp%ndat = ndat
    tmp%dat(1:ndat) = x%dat(1:ndat)
    call vectori4b_destroy(x)
    call vectori4b_create(2*nelem,x)
    x%ndat = ndat
    x%dat(1:ndat) = tmp%dat(1:ndat)
    call vectori4b_destroy(tmp)
  
  end subroutine vectori4b_double
  
  subroutine vectori4b_cat_scalar(x, a)
  
    type(vectori4b), intent(inout) :: x 
    integer(i4b), intent(in)       :: a
    
    integer(i4b) :: ndat
    
    ndat = x%ndat + 1
    if (ndat > x%nelem) then
       call vectori4b_double(x)
    end if
    x%ndat = ndat
    x%dat(ndat) = a
    
  end subroutine vectori4b_cat_scalar
  
  subroutine vectori4b_cat_vector(x, a)
  
    type(vectori4b), intent(inout) :: x 
    type(vectori4b), intent(in)    :: a
    
    integer(i4b) :: ndat
    
    ndat = x%ndat + a%ndat
    do while (ndat > x%nelem)
       call vectori4b_double(x)
    end do
    x%dat(x%ndat+1:ndat) = a%dat(1:a%ndat)
    x%ndat = ndat
    
  end subroutine vectori4b_cat_vector

  subroutine vectori4b_assign(x, a)
  
    type(vectori4b), intent(inout)         :: x 
    integer(i4b), dimension(:), intent(in) :: a
    
    integer(i4b) :: ndat
    
    ndat = size(a,1)
    do while (ndat > x%nelem)
       call vectori4b_double(x)
    end do
    x%dat(1:ndat) = a
    x%ndat = ndat
    
  end subroutine vectori4b_assign

  subroutine vectori4b_assign_vector(x, a)
  
    type(vectori4b), intent(inout) :: x 
    type(vectori4b), intent(in)    :: a
    
    integer(i4b) :: ndat
    
    if (a%nelem /= x%nelem) then
       deallocate(x%dat)
       allocate(x%dat(a%nelem))
       x%nelem = a%nelem
    end if
    x%ndat = a%ndat
    x%dat = a%dat
    
  end subroutine vectori4b_assign_vector
  
  subroutine vectori4b_clear(x)
  
    type(vectori4b), intent(inout) :: x
    
    x%ndat = 0
    
  end subroutine vectori4b_clear
  
  function vectori4b_size(x) result(ndat)
  
    type(vectori4b), intent(inout) :: x 
    
    integer(i4b) :: ndat
    
    ndat = x%ndat
    
  end function vectori4b_size
    

end module vectors
