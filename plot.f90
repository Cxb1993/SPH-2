PROGRAM plot

  USE types
  USE vars, ONLY : vars_create,vars_read,vars_write,vars_destroy,vars_params,vars_step,vars_density
  USE vars, ONLY : x,m,rho
  USE kernels, ONLY : kernel_select,kernel_support
  USE grid, ONLY : grid_create,grid_destroy
  USE mesh, ONLY : mesh_create,mesh_destroy,mesh_interp,mesh_write
  IMPLICIT NONE
  
  INTEGER(I4B) :: ndim,np,i,dump,nsteps
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xmin,xmax,grav
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: periodic
  REAL(DP) :: h,p0,rho0,gam,dt,simtime,cputime1,cputime2
  CHARACTER(LEN=80) :: kernelname,dir,job
  CHARACTER(LEN=120) :: filnam
  
  OPEN(10,FILE='input')
  READ(10,'(A)')dir
  READ(10,'(A)')job
  READ(10,*)ndim
  READ(10,*)np
  ALLOCATE(xmin(ndim),xmax(ndim),periodic(ndim),grav(ndim))
  READ(10,*)xmin
  READ(10,*)xmax
  READ(10,*)periodic
  READ(10,*)h
  READ(10,'(A)')kernelname
  READ(10,*)rho0
  READ(10,*)gam
  READ(10,*)p0
  READ(10,*)grav
  READ(10,*)dt
  READ(10,*)nsteps
  READ(10,*)dump
  READ(10,*)simtime
  CLOSE(10)
  
  DO i=0,dump
     WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),i
     PRINT '(A)',TRIM(filnam)
     CALL vars_create(ndim,np)
     CALL vars_read(TRIM(filnam))
     CALL vars_params(h,gam,rho0,p0,grav)
     CALL kernel_select(kernelname)
     CALL grid_create(xmin,xmax,kernel_support(h),periodic)
     CALL mesh_create(h,xmin,xmax)
     CALL mesh_interp(x,m,rho,rho)
     CALL mesh_write(filnam)
     CALL mesh_destroy
     CALL grid_destroy
     CALL vars_destroy
  END DO
  DEALLOCATE(xmin,xmax,periodic,grav)
  
END PROGRAM plot
