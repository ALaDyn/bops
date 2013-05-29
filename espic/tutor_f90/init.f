!
!   initialise main particle and grid variables
!
subroutine init

include 'es.h'


 ntrun = 500     ! # timesteps
 nx = 10        ! # grid points
 ne = 5000       ! # electrons
 ni = 0          ! # ions (fixed)


 grid_length = 10.   ! size of spatial grid

 dt = 0.2        ! normalised timestep
     
 q_over_me=-1.   ! electron charge:mass ratio

 rho0 = 1.0      ! background ion density
 vte = 0.1       ! thermal velocity

 ilas = 0        ! laser switch - 0 = off, 1 = uniform sine wave
 w0=1.0          ! laser frequency
 a0 = 0.1        ! laser amplitude

 bc_field = 2     ! field boundary conditions:  1 = periodic, 2 = reflective
 bc_particle = 2  ! particle BCs: 1 = periodic, 2 = reflective, 3 = thermal

 ihist = 10       ! frequency of time-histories output
 igraph = 500      ! freq. of graphical snapshots
 iout = 10       ! freq. of printed diags.

 itime = 0       ! initialise time counter

!  Open file for history plots

 open(60,file='hist.data')

!  ------------------
!  derived parameters
!  ------------------

 dx = grid_length/nx       ! mesh size

 omega_p = sqrt(rho0)      ! plasma frequency
 x_debye = vte/omega_p     ! Debye length
 xdodx = xdebye/dx         ! ratio
 ncell = ne/nx             ! # particles per cell
  
 write(6,*) '# particles = ',ne
 write(6,*) '# mesh points = ',nx
 write(6,*) '# particles/cell = ',ncell

 write(6,*) 'grid length = ',grid_length
 write(6,*) 'thermal velocity: ',vte
 write(6,*) 'mesh size = ',dx
 write(6,*) 'Debye length = ',x_debye

 write(6,*) 'timestep = ',dt
 write(6,*) '# timesteps = ',ntrun
 write(6,*) 'run time = ',dt*ntrun


end
