module mesh

  use types
  implicit none

  integer(i4b)                            :: ndim, n
  real(dp)                                :: h
  integer(i4b), dimension(:), allocatable :: nx
  real(dp), dimension(:), allocatable     :: dx, u
  real(dp), dimension(:,:), allocatable   :: x

contains

  subroutine mesh_create(h_, xmin, xmax)

    real(dp), intent(in)               :: h_
    real(dp), dimension(:), intent(in) :: xmin, xmax

    real(dp), dimension(:), allocatable     :: xdiff
    integer(i4b), dimension(:), allocatable :: i
    integer(i4b)                            :: j,k,q

    h = h_
    ndim = size(xmin,1)
    allocate(nx(ndim), dx(ndim), i(ndim), xdiff(ndim))
    xdiff = xmax - xmin - 2.d-8
    nx = ceiling(xdiff/h)
    dx = xdiff/(nx-1)
    n = nx(1)
    do j=2,ndim
       n = n*nx(j)
    end do
    allocate(x(ndim,n), u(n))
    
    do j=1,n                          ! loop through all points
       call mesh_index(j, i)
       x(:,j) = dx*(i - 1) + xmin + 1.d-8    ! calculate point in grid
    end do

    deallocate(i, xdiff)

  end subroutine mesh_create
  
  subroutine mesh_destroy()
  
    deallocate(x, u, nx, dx)
    
  end subroutine mesh_destroy
  
  subroutine mesh_index(j, i)
  
    integer(i4b), intent(in)                   :: j
    integer(i4b), dimension(ndim), intent(out) :: i
    
    integer(i4b) :: k, q

    q = j - 1
    do k=1,ndim
       i(k) = mod(q,nx(k)) + 1
       q = (q - i(k) + 1)/nx(k)
    end do

  end subroutine mesh_index
  
  subroutine mesh_interp(xp, mp, rhop, up)
  
    use vectors
    use grid, only : grid_clear, grid_insert, grid_neighbours_cell, maptogrid_index,map
  
    real(dp), dimension(:,:), intent(in) :: xp
    real(dp), dimension(:), intent(in)   :: mp, rhop, up
    
    integer(i4b)                               :: i
    type(vectori4b)                            :: list
    integer(i4b), dimension(ndim)              :: igrid
    type(vectori4b), dimension(:), allocatable :: meshmap

! map mesh points into grid

    call grid_clear
    call grid_insert(x)

    allocate(meshmap(size(map,1))) ! copy map to meshmap
    do i=1,size(map,1)
       call vectori4b_create(1, meshmap(i))
       call vector_assign(meshmap(i), map(i))
    end do

! map particles into grid

    call grid_clear
    call grid_insert(xp)   

! loop through grid cells
 
    call vectori4b_create(1, list)
    do i=1,size(map,1)
      call maptogrid_index(i, igrid)
      call vectori4b_clear(list)
      call grid_neighbours_cell(1, igrid, list)
      call mesh_interact(meshmap(i)%dat(1:meshmap(i)%ndat), list%dat(1:list%ndat), xp, mp, rhop, up)
    end do
    call vectori4b_destroy(list)

! tidy up

    do i=1,size(meshmap, 1)
       call vectori4b_destroy(meshmap(i))
    end do
    deallocate(meshmap)
    
  end subroutine mesh_interp
  
  subroutine mesh_interact(a, b, xb, mb, rhob, ub)
  
    use kernels, only : kernel_support, kernel_eval
    use grid, only : xdiff
  
    integer(i4b), dimension(:), intent(in) :: a, b
    real(dp), dimension(:,:), intent(in)   :: xb
    real(dp), dimension(:), intent(in)     :: mb, rhob, ub
    
    integer(i4b)              :: na, nb, i, j
    real(dp), dimension(ndim) :: rab
    real(dp)                  :: supp, twosupp, supp2, rab2
    
    na = size(a, 1)
    nb = size(b, 1)
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
    do i=1,na
       u(a(i)) = 0.0_dp
       do j=1,nb
          rab = x(:,a(i)) - xb(:,b(j))
          where (rab > twosupp)
             rab = rab - xdiff
          elsewhere (rab < -twosupp)
             rab = rab + xdiff
          end where   
          rab2 = sum(rab**2)
          if (rab2 < supp2) then
             u(a(i)) = u(a(i)) + mb(b(j))/rhob(b(j))*ub(b(j))*kernel_eval(rab,h)
          end if
       end do
    end do
  
  end subroutine mesh_interact

  subroutine mesh_write(filename)

    character(*), intent(in) :: filename

    integer(i4b), parameter     :: dbg=0, dbl=1
    integer(i4b)                :: i, j, m(3)
    character(len=1), parameter :: nullc = char(0)
    character(len=80)           :: tmpc
    real(dp), dimension(ndim+1) :: tmp
    integer(i4b)                :: tecini, teczne, tecdat, tecend, ierr

    write(tmpc, '(a,i1)') 'x', 1
    do i=2,ndim
       write(tmpc, '(a,a,i1)') trim(tmpc), ' x', i
    end do
    write(tmpc, '(a,a)') trim(tmpc), ' u'
    ierr = tecini('sph'//nullc, trim(tmpc)//nullc, trim(filename)//'.plt'//nullc, '.'//nullc, dbg, dbl)
    m = 1
    m(1:ndim) = nx
    ierr = teczne('sph'//nullc, m(1), m(2), m(3), 'point'//nullc, nullc)
    do j=1,n
       do i=1,ndim
          tmp(i) = x(i,j)
       end do
       tmp(ndim+1) = u(j)
       ierr = tecdat(ndim+1,tmp,dbl)
    end do
    ierr = tecend()

  end subroutine mesh_write

end module mesh
