module grid

  use types
  use vectors
  implicit none
  
  interface grid_neighbours
    module procedure grid_neighbours_cell, grid_neighbours_point
  end interface

! if you create a type with all this data in, then you can have multiple data sets
  
  integer(i4b)                               :: ndim,n
  integer(i4b), dimension(:), allocatable    :: nx, periodic
  real(dp), dimension(:), allocatable        :: xmin, xmax, xdiff
  real(dp)                                   :: dx
  type(vectori4b), dimension(:), allocatable :: map
  
  ! map contains list of particles in each cell

contains

  subroutine grid_create(x1, x2, h, p)
  
    real(dp), dimension(:), intent(in)     :: x1, x2
    real(dp), intent(in)                   :: h
    integer(i4b), dimension(:), intent(in) :: p
    
    integer(i4b) :: i
    
    ndim = size(x1,1)
    allocate(nx(ndim), xmin(ndim), xmax(ndim), xdiff(ndim), periodic(ndim))
    xmin = x1
    xmax = x2
    xdiff = xmax - xmin
    dx = h
    nx = ceiling(xdiff/dx)
    periodic = p
    
    n = nx(1)
    do i=2,ndim
       n = n*nx(i)
    end do
    allocate(map(n))
    do i=1,n
       call vectori4b_create(1, map(i))
    end do
    
  end subroutine grid_create
  
  subroutine grid_destroy
  
    integer(i4b) :: i
  
    do i=1,n
       call vectori4b_destroy(map(i))
    end do
    deallocate(nx, xmin, xmax, xdiff, periodic, map)
    
  end subroutine grid_destroy
  
  function grid_index(x) result(i)
  
    real(dp), dimension(:), intent(in) :: x
    
    integer(i4b), dimension(size(x,1)) :: i
    
    i = ceiling((x - xmin)/dx)
    
  end function grid_index
  
  function gridtomap_index(i) result(j)
  
    integer(i4b), dimension(ndim), intent(in) :: i
    
    integer(i4b) :: j, k
    
    j = i(ndim)
    do k=ndim-1,1,-1
       j = (j - 1)*nx(k) + i(k)
    end do 

  end function gridtomap_index
  
  subroutine maptogrid_index(j, i)
  
    integer(i4b), intent(in)                   :: j
    integer(i4b), dimension(ndim), intent(out) :: i
    
    integer(i4b), dimension(ndim) :: p
    integer(i4b)                  :: k, q

    q = j - 1
    do k=1,ndim
       i(k) = mod(q,nx(k)) + 1
       q = (q - i(k) + 1)/nx(k)
    end do
    
    
!    p(1) = j - 1
!    do k=2,ndim
!       p(k) = p(k-1)/nx(k)
!    end do
!    i(ndim) = p(ndim) + 1
!    do k=ndim-1,1,-1
!       i(k) = p(k) - p(k+1)*nx(k) + 1
!    end do

  end subroutine maptogrid_index
  
  subroutine grid_insert(x)
  
    real(dp), dimension(:,:), intent(in) :: x
    
    integer(i4b) :: ndim,np,i,j
    
    ndim = size(x,1)
    np = size(x,2)

    do i=1,np
       j = gridtomap_index(grid_index(x(:,i)))
       call vector_cat(map(j),i)
    end do
    
  end subroutine grid_insert
  
  subroutine grid_clear()
  
    integer(i4b) :: i
  
    do i=1,n
       call vectori4b_clear(map(i))
    end do
    
  end subroutine grid_clear
  
  subroutine grid_neighbours_point(x, list)
  
    real(dp), dimension(ndim), intent(in) :: x ! particle position to find in grid
    type(vectori4b), intent(inout)        :: list
    
    call vectori4b_clear(list)
    call grid_neighbours_cell(1, grid_index(x), list)
    
  end subroutine grid_neighbours_point
    
  recursive subroutine grid_neighbours_cell(i, igrid, list)
  
! find all neighbouring cells and add particles in those cells to list
    
    integer(i4b), intent(in)                  :: i
    integer(i4b), dimension(ndim), intent(in) :: igrid
    type(vectori4b), intent(inout)            :: list
      
    integer(i4b)                  :: j,jmin,jmax
    integer(i4b), dimension(ndim) :: jgrid
   
    jgrid = igrid
    jmin = igrid(i) - 1
    jmax = igrid(i) + 1
    if (periodic(i) == 1) then
       if (jmin == 0) jmin = -1
       if (jmax == nx(i)) jmax = nx(i) + 1
    end if
    do j=jmin,jmax
       if (j < 1) then
          if (periodic(i) == 0) then
             cycle
          else
             jgrid(i) = j + nx(i)
          end if
       else if (j > nx(i)) then
          if (periodic(i) == 0) then
             cycle
          else
             jgrid(i) = j - nx(i)
          end if
       else
          jgrid(i) = j   
       end if
       if (i == ndim) then
          call vector_cat(list, map(gridtomap_index(jgrid)))
       else        
          call grid_neighbours_cell(i+1, jgrid, list)
       end if    
    end do
  
  end subroutine grid_neighbours_cell
  
  subroutine grid_loop(interact)
  
! here we loop over all cells. send to func lists of particles in a cell and neighbour list

    interface
      subroutine interact(a,b)
        use types
        integer(i4b), dimension(:), intent(in) :: a, b
      end subroutine interact
    end interface

    integer(i4b)                  :: i
    type(vectori4b)               :: list
    integer(i4b), dimension(ndim) :: igrid
    
    !$omp parallel private(list, igrid)
    call vectori4b_create(1,list)
    !$omp do
    do i=1,n
      call maptogrid_index(i,igrid)
  !    print *,i,nx
  !    print *,i,igrid,gridtomap_index(igrid)
      call vectori4b_clear(list)
      call grid_neighbours_cell(1, igrid,list)
      call interact(map(i)%dat(1:map(i)%ndat), list%dat(1:list%ndat))
    end do
    !$omp end do
    call vectori4b_destroy(list)
    !$omp end parallel

  end subroutine grid_loop
    
  subroutine grid_periodic_map(x)
  
    real(dp), dimension(:,:), intent(inout) :: x
    
    integer(i4b) :: i,np
    
    np = size(x,2)
    
    do i=1,np
!       x(:,i) = modulo(x(:,i)-xmin,xdiff) + xmin
       where (x(:,i) > xmax)
          x(:,i) = x(:,i) - xdiff
       else where (x(:,i) < xmin)
          x(:,i) = x(:,i) + xdiff
       end where
    end do
  
  end subroutine grid_periodic_map

end module grid
