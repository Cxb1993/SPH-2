PROGRAM main

  USE types
  USE vars, ONLY : vars_create,vars_read,vars_write,vars_destroy,vars_params,vars_step,vars_density
  USE kernels, ONLY : kernel_select,kernel_support
  USE grid, ONLY : grid_create,grid_destroy
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
  
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),dump
  CALL vars_create(ndim,np)
  CALL vars_read(TRIM(filnam))
  CALL vars_params(h,gam,rho0,p0,grav)
  CALL kernel_select(kernelname)
  CALL grid_create(xmin,xmax,kernel_support(h),periodic)

 ! IF (dump == 0) CALL vars_density
  DO i=1,nsteps
     CALL CPU_TIME(cputime1)
     CALL vars_step(dt)
     CALL CPU_TIME(cputime2)
     PRINT '(A,E16.10)','CPU time: ',cputime2-cputime1
 !    CALL vars_write('temp.dat')
  END DO
  dump = dump + 1
  simtime = simtime + dt*nsteps
  WRITE(filnam,'(4A,I3.3)')TRIM(dir),TRIM(job),'/',TRIM(job),dump
  CALL vars_write(TRIM(filnam))

  OPEN(10,FILE='input')
  READ(10,*)dir
  READ(10,*)job
  READ(10,*)ndim
  READ(10,*)np
  READ(10,*)xmin
  READ(10,*)xmax
  READ(10,*)periodic
  READ(10,*)h
  READ(10,*)kernelname
  READ(10,*)rho0
  READ(10,*)gam
  READ(10,*)p0
  READ(10,*)grav
  READ(10,*)dt
  READ(10,*)nsteps
  WRITE(10,*)dump
  WRITE(10,*)simtime
  CLOSE(10)

  CALL grid_destroy
  CALL vars_destroy
  DEALLOCATE(xmin,xmax,periodic,grav)
  
END PROGRAM main  
