MODULE mesh

  USE types
  IMPLICIT NONE

  INTEGER(I4B) :: ndim,n
  REAL(DP) :: h
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: nx
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx,u
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: x

CONTAINS

  SUBROUTINE mesh_create(h_,xmin,xmax)

    REAL(DP), INTENT(IN) :: h_
    REAL(DP), DIMENSION(:), INTENT(IN) :: xmin,xmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: xdiff
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: i
    INTEGER(I4B) :: j,k,q

    h = h_
    ndim = SIZE(xmin,1)
    ALLOCATE(nx(ndim),dx(ndim),i(ndim),xdiff(ndim))
    xdiff = xmax - xmin - 2.D-8
    nx = CEILING(xdiff/h)
    dx = xdiff/(nx-1)
    n = nx(1)
    DO j=2,ndim
       n = n*nx(j)
    END DO
    ALLOCATE(x(ndim,n),u(n))
    
    DO j=1,n                          ! loop through all points
       CALL mesh_index(j,i)
       x(:,j) = dx*(i - 1) + xmin + 1.D-8    ! calculate point in grid
    END DO

    DEALLOCATE(i,xdiff)

  END SUBROUTINE mesh_create
  
  SUBROUTINE mesh_destroy()
  
    DEALLOCATE(x,u,nx,dx)
    
  END SUBROUTINE mesh_destroy
  
  SUBROUTINE mesh_index(j,i)
  
    INTEGER(I4B), INTENT(IN) :: j
    INTEGER(I4B), DIMENSION(ndim), INTENT(OUT) :: i
    
    INTEGER(I4B) :: k,q

    q = j - 1
    DO k=1,ndim
       i(k) = MOD(q,nx(k)) + 1
       q = (q - i(k) + 1)/nx(k)
    END DO

  END SUBROUTINE mesh_index
  
  SUBROUTINE mesh_interp(xp,mp,rhop,up)
  
    USE vectors
    USE grid, ONLY : grid_clear,grid_insert,grid_neighbours_cell,maptogrid_index,map
  
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: xp
    REAL(DP), DIMENSION(:), INTENT(IN) :: mp,rhop,up
    
    INTEGER(I4B) :: i
    TYPE(vectori4b) :: list
    INTEGER(I4B), DIMENSION(ndim) :: igrid
    TYPE(vectori4b), DIMENSION(:), ALLOCATABLE :: meshmap

! map mesh points into grid

    CALL grid_clear
    CALL grid_insert(x)

    ALLOCATE(meshmap(SIZE(map,1))) ! copy map to meshmap
    DO i=1,SIZE(map,1)
       CALL vectori4b_create(1,meshmap(i))
       CALL vector_assign(meshmap(i),map(i))
    END DO

! map particles into grid

    CALL grid_clear
    CALL grid_insert(xp)   

! loop through grid cells
 
    CALL vectori4b_create(1,list)
    DO i=1,SIZE(map,1)
      CALL maptogrid_index(i,igrid)
      CALL vectori4b_clear(list)
      CALL grid_neighbours_cell(1,igrid,list)
      CALL mesh_interact(meshmap(i)%dat(1:meshmap(i)%ndat),list%dat(1:list%ndat),xp,mp,rhop,up)
    END DO
    CALL vectori4b_destroy(list)

! tidy up

    DO i=1,SIZE(meshmap,1)
       CALL vectori4b_destroy(meshmap(i))
    END DO
    DEALLOCATE(meshmap)
    
  END SUBROUTINE mesh_interp
  
  SUBROUTINE mesh_interact(a,b,xb,mb,rhob,ub)
  
    USE kernels, ONLY : kernel_support,kernel_eval
    USE grid, ONLY : xdiff
  
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: xb
    REAL(DP), DIMENSION(:), INTENT(IN) :: mb,rhob,ub
    
    INTEGER(I4B) :: na,nb,i,j
    REAL(DP), DIMENSION(ndim) :: rab
    REAL(DP) :: supp,twosupp,supp2,rab2
    
    na = SIZE(a,1)
    nb = SIZE(b,1)
    supp = kernel_support(h)
    twosupp = 2.0_dp*supp
    supp2 = supp**2
    DO i=1,na
       u(a(i)) = 0.0_dp
       DO j=1,nb
          rab = x(:,a(i)) - xb(:,b(j))
          WHERE (rab > twosupp)
             rab = rab - xdiff
          ELSEWHERE (rab < -twosupp)
             rab = rab + xdiff
          END WHERE   
          rab2 = SUM(rab**2)
          IF (rab2 < supp2) THEN
             u(a(i)) = u(a(i)) + mb(b(j))/rhob(b(j))*ub(b(j))*kernel_eval(rab,h)
          END IF
       END DO
    END DO
  
  END SUBROUTINE mesh_interact

  SUBROUTINE mesh_write(filename)

    CHARACTER(*), INTENT(IN) :: filename

    INTEGER(I4B), PARAMETER :: dbg=0, dbl=1
    INTEGER(I4B) :: i,j,m(3)
    CHARACTER(LEN=1), PARAMETER :: nullc = CHAR(0)
    CHARACTER(LEN=80) :: tmpc
    REAL(DP), DIMENSION(ndim+1) :: tmp
    INTEGER(I4B) :: tecini,teczne,tecdat,tecend,ierr

    WRITE(tmpc,'(A,I1)')'X',1
    DO i=2,ndim
       WRITE(tmpc,'(A,A,I1)')TRIM(tmpc),' X',i
    END DO
    WRITE(tmpc,'(A,A)')TRIM(tmpc),' U'
    ierr = tecini('SPH'//nullc,TRIM(tmpc)//nullc,TRIM(filename)//'.plt'//nullc,'.'//nullc,dbg,dbl)
    m = 1
    m(1:ndim) = nx
    ierr = teczne('SPH'//nullc,m(1),m(2),m(3),'POINT'//nullc,nullc)
    DO j=1,n
       DO i=1,ndim
          tmp(i) = x(i,j)
       END DO
       tmp(ndim+1) = u(j)
       ierr = tecdat(ndim+1,tmp,dbl)
    END DO
    ierr = tecend()

  END SUBROUTINE mesh_write

END MODULE mesh
