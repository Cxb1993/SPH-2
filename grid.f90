MODULE grid

  USE types
  USE vectors
  IMPLICIT NONE
  
  INTERFACE grid_neighbours
    MODULE PROCEDURE grid_neighbours_cell, grid_neighbours_point
  END INTERFACE

! if you create a type with all this data in, then you can have multiple data sets
  
  INTEGER(I4B) :: ndim,n
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nx,periodic
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xmin,xmax,xdiff
  REAL(DP) :: dx
  TYPE(vectori4b), DIMENSION(:), ALLOCATABLE :: map
  
  ! map contains list of particles in each cell

CONTAINS

  SUBROUTINE grid_create(x1,x2,h,p)
  
    REAL(DP), DIMENSION(:), INTENT(IN) :: x1,x2
    REAL(DP), INTENT(IN) :: h
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: p
    
    INTEGER(I4B) :: i
    
    ndim = SIZE(x1,1)
    ALLOCATE(nx(ndim),xmin(ndim),xmax(ndim),xdiff(ndim),periodic(ndim))
    xmin = x1
    xmax = x2
    xdiff = xmax - xmin
    dx = h
    nx = CEILING(xdiff/dx)
    periodic = p
    
    n = nx(1)
    DO i=2,ndim
       n = n*nx(i)
    END DO
    ALLOCATE(map(n))
    DO i=1,n
       CALL vectori4b_create(1,map(i))
    END DO
    
  END SUBROUTINE grid_create
  
  SUBROUTINE grid_destroy
  
    INTEGER(I4B) :: i
  
    DO i=1,n
       CALL vectori4b_destroy(map(i))
    END DO
    DEALLOCATE(nx,xmin,xmax,xdiff,periodic,map)
    
  END SUBROUTINE grid_destroy
  
  FUNCTION grid_index(x) RESULT(i)
  
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    
    INTEGER(I4B), DIMENSION(SIZE(x,1)) :: i
    
    i = CEILING((x - xmin)/dx)
    
  END FUNCTION grid_index
  
  FUNCTION gridtomap_index(i) RESULT(j)
  
    INTEGER(I4B), DIMENSION(ndim), INTENT(IN) :: i
    INTEGER(I4B) :: j
    
    INTEGER(I4B) :: k
    
    j = i(ndim)
    DO k=ndim-1,1,-1
       j = (j - 1)*nx(k) + i(k)
    END DO 

  END FUNCTION gridtomap_index
  
  SUBROUTINE maptogrid_index(j,i)
  
    INTEGER(I4B), INTENT(IN) :: j
    INTEGER(I4B), DIMENSION(ndim), INTENT(OUT) :: i
    
    INTEGER(I4B), DIMENSION(ndim) :: p
    INTEGER(I4B) :: k,q

    q = j - 1
    DO k=1,ndim
       i(k) = MOD(q,nx(k)) + 1
       q = (q - i(k) + 1)/nx(k)
    END DO
    
    
!    p(1) = j - 1
!    DO k=2,ndim
!       p(k) = p(k-1)/nx(k)
!    END DO
!    i(ndim) = p(ndim) + 1
!    DO k=ndim-1,1,-1
!       i(k) = p(k) - p(k+1)*nx(k) + 1
!    END DO

  END SUBROUTINE maptogrid_index
  
  SUBROUTINE grid_insert(x)
  
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
    
    INTEGER(I4B) :: ndim,np,i,j
    
    ndim = SIZE(x,1)
    np = SIZE(x,2)

    DO i=1,np
       j = gridtomap_index(grid_index(x(:,i)))
       CALL vector_cat(map(j),i)
    END DO
    
  END SUBROUTINE grid_insert
  
  SUBROUTINE grid_clear()
  
    INTEGER(I4B) :: i
  
    DO i=1,n
       CALL vectori4b_clear(map(i))
    END DO
    
  END SUBROUTINE grid_clear
  
  SUBROUTINE grid_neighbours_point(x,list)
  
    REAL(DP), DIMENSION(ndim), INTENT(IN) :: x ! particle position to find in grid
    TYPE(vectori4b), INTENT(INOUT) :: list
    
    CALL vectori4b_clear(list)
    CALL grid_neighbours_cell(1,grid_index(x),list)
    
  END SUBROUTINE grid_neighbours_point
    
  RECURSIVE SUBROUTINE grid_neighbours_cell(i,igrid,list)
  
! find all neighbouring cells and add particles in those cells to list
    
    INTEGER(I4B), INTENT(IN) :: i
    INTEGER(I4B), DIMENSION(ndim), INTENT(IN) :: igrid
    TYPE(vectori4b), INTENT(INOUT) :: list
      
    INTEGER(I4B) :: j,jmin,jmax
    INTEGER(I4B), DIMENSION(ndim) :: jgrid
   
    jgrid = igrid
    jmin = igrid(i) - 1
    jmax = igrid(i) + 1
    IF (periodic(i) == 1) THEN
       IF (jmin == 0) jmin = -1
       IF (jmax == nx(i)) jmax = nx(i) + 1
    END IF
    DO j=jmin,jmax
       IF (j < 1) THEN
          IF (periodic(i) == 0) THEN
             CYCLE
          ELSE
             jgrid(i) = j + nx(i)
          END IF
       ELSE IF (j > nx(i)) THEN
          IF (periodic(i) == 0) THEN
             CYCLE
          ELSE
             jgrid(i) = j - nx(i)
          END IF
       ELSE
          jgrid(i) = j   
       END IF
       IF (i == ndim) THEN
          CALL vector_cat(list,map(gridtomap_index(jgrid)))
       ELSE        
          CALL grid_neighbours_cell(i+1,jgrid,list)
       END IF    
    END DO
  
  END SUBROUTINE grid_neighbours_cell
  
  SUBROUTINE grid_loop(interact)
  
! Here we loop over all cells. Send to func lists of particles in a cell and neighbour list

    INTERFACE
      SUBROUTINE interact(a,b)
        USE types
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
      END SUBROUTINE interact
    END INTERFACE

    INTEGER(I4B) :: i
    TYPE(vectori4b) :: list
    INTEGER(I4B), DIMENSION(ndim) :: igrid
    
    CALL vectori4b_create(1,list)
    DO i=1,n
      CALL maptogrid_index(i,igrid)
  !    print *,i,nx
  !    print *,i,igrid,gridtomap_index(igrid)
      CALL vectori4b_clear(list)
      CALL grid_neighbours_cell(1,igrid,list)
      CALL interact(map(i)%dat(1:map(i)%ndat),list%dat(1:list%ndat))
    END DO
    CALL vectori4b_destroy(list)
  
  END SUBROUTINE grid_loop
    
  SUBROUTINE grid_periodic_map(x)
  
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: x
    
    INTEGER(I4B) :: i,np
    
    np = SIZE(x,2)
    
    DO i=1,np
!       x(:,i) = MODULO(x(:,i)-xmin,xdiff) + xmin
       WHERE (x(:,i) > xmax)
          x(:,i) = x(:,i) - xdiff
       ELSE WHERE (x(:,i) < xmin)
          x(:,i) = x(:,i) + xdiff
       END WHERE
    END DO
  
  END SUBROUTINE grid_periodic_map

END MODULE grid
