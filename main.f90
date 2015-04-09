program main

  use types
  use vars, only : vars_create, vars_read, vars_write, vars_destroy
  use vars, only : vars_params, vars_step, vars_density
  use kernels, only : kernel_select, kernel_support
  use grid, only : grid_create, grid_destroy
  implicit none
  
  integer(i4b)                            :: ndim, np, i, dump, nsteps
  real(dp), dimension(:), allocatable     :: xmin, xmax, grav
  integer(i4b), dimension(:), allocatable :: periodic
  real(dp)                                :: h, p0, rho0, gam, dt, simtime
  integer                                 :: tic, toc, ticrate
  character(len=80)                       :: kernelname, dir, job
  character(len=120)                      :: filnam
  
  open(10,file='input')
  read(10,'(a)')dir
  read(10,'(a)')job
  read(10,*)ndim
  read(10,*)np
  allocate(xmin(ndim), xmax(ndim), periodic(ndim), grav(ndim))
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
  
  write(filnam, '(4a,i3.3)') trim(dir), trim(job), '/', trim(job),dump
  call vars_create(ndim, np)
  call vars_read(trim(filnam))
  call vars_params(h, gam, rho0, p0, grav)
  call kernel_select(kernelname)
  call grid_create(xmin, xmax, kernel_support(h), periodic)

 ! if (dump == 0) call vars_density
  do i=1,nsteps
     call system_clock(tic, ticrate)
     call vars_step(dt)
     call system_clock(toc, ticrate)
     print '(a,e16.10)','wall clock time: ',real(toc - tic)/ticrate
 !    call vars_write('temp.dat')
  end do
  dump = dump + 1
  simtime = simtime + dt*nsteps
  write(filnam, '(4a,i3.3)') trim(dir), trim(job), '/', trim(job), dump
  call vars_write(trim(filnam))

  open(10,file='input')
  read(10,*)dir
  read(10,*)job
  read(10,*)ndim
  read(10,*)np
  read(10,*)xmin
  read(10,*)xmax
  read(10,*)periodic
  read(10,*)h
  read(10,*)kernelname
  read(10,*)rho0
  read(10,*)gam
  read(10,*)p0
  read(10,*)grav
  read(10,*)dt
  read(10,*)nsteps
  write(10,*)dump
  write(10,*)simtime
  close(10)

  call grid_destroy
  call vars_destroy
  deallocate(xmin, xmax, periodic, grav)
  
end program main  
