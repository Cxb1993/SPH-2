program plot

  use types
  use vars, only : vars_create,vars_read,vars_write,vars_destroy,vars_params,vars_step,vars_density
  use vars, only : x,m,rho
  use kernels, only : kernel_select,kernel_support
  use grid, only : grid_create,grid_destroy
  use mesh, only : mesh_create,mesh_destroy,mesh_interp,mesh_write
  implicit none
  
  integer(i4b) :: ndim,np,i,dump,nsteps
  real(dp), dimension(:), allocatable :: xmin,xmax,grav
  integer(i4b), dimension(:), allocatable :: periodic
  real(dp) :: h,p0,rho0,gam,dt,simtime,cputime1,cputime2
  character(len=80) :: kernelname,dir,job
  character(len=120) :: filnam
  
  open(10,file='input')
  read(10,'(a)')dir
  read(10,'(a)')job
  read(10,*)ndim
  read(10,*)np
  allocate(xmin(ndim),xmax(ndim),periodic(ndim),grav(ndim))
  read(10,*)xmin
  read(10,*)xmax
  read(10,*)periodic
  read(10,*)h
  read(10,'(a)')kernelname
  read(10,*)rho0
  read(10,*)gam
  read(10,*)p0
  read(10,*)grav
  read(10,*)dt
  read(10,*)nsteps
  read(10,*)dump
  read(10,*)simtime
  close(10)
  
  do i=0,dump
     write(filnam,'(4a,i3.3)')trim(dir),trim(job),'/',trim(job),i
     print '(a)',trim(filnam)
     call vars_create(ndim,np)
     call vars_read(trim(filnam))
     call vars_params(h,gam,rho0,p0,grav)
     call kernel_select(kernelname)
     call grid_create(xmin,xmax,kernel_support(h),periodic)
     call mesh_create(h,xmin,xmax)
     call mesh_interp(x,m,rho,rho)
     call mesh_write(filnam)
     call mesh_destroy
     call grid_destroy
     call vars_destroy
  end do
  deallocate(xmin,xmax,periodic,grav)
  
end program plot
