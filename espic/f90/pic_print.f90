!  *************************************************
!
!    Electrostatic PIC code for Heraeus School 
!    laser-plasma simulation practical
!
!    1d1v,  slab geometry, fortran-90 version
!  
!    Paul Gibbon, August 2000
!
!
!  *************************************************

program pices

include 'es.h'			! common variables

 call init			! initialisation

 call loadx			! load particles onto grid
 call loadv			! define velocity distribution

 call density			! compute initial density from particles

 call field			! compute initial electric field

 call diagnostics		! output initial conditions

 
 do itime=1,ntrun

    call push			! Push particles
    call boundaries		! particle boundary conditions

    call density		! compute density

    call field			! compute electric field (Poisson)

    call diagnostics		! output snapshots and time-histories

 end do

!  close time-hist files
 close(60)
end
!
!   initialise main particle and grid variables
!
subroutine init

include 'es.h'


 ntrun = 200     ! # timesteps
 nx = 100        ! # grid points
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

 bc_field = 1     ! field boundary conditions:  1 = periodic, 2 = reflective
 bc_particle = 1  ! particle BCs: 1 = periodic, 2 = reflective, 3 = thermal

 ihist = 1       ! frequency of time-history output
 igraph = 100      ! freq. of graphical snapshots
 iout = 10        ! freq. of printed diags.

 itime = 0       ! initialise time counter

!  Open file for history plots

 open(60,file='hist.dat')

!  ------------------
!  derived parameters
!  ------------------

 dx = grid_length/nx     ! mesh size

 omega_p = sqrt(rho0)    ! plasma frequency
 x_debye = vte/omega_p        ! Debye length
 xdodx = xdebye/dx       ! ratio
 ncell = ne/nx  ! # particles per cell
  
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

!  =======================================
!
!    Particle loading for uniform density
!
!  =======================================

subroutine loadx

include 'es.h'


 if (bc_particle.eq.2) then

!  for reflective bcs to work properly, must place particle
!  boundaries half a mesh spacing inside field boundaries

      wall_left = dx/2.
      wall_right = grid_length-dx/2.
 else 

!  periodic boundaries

      wall_left = 0.
      wall_right = grid_length
 endif

 xload = wall_right - wall_left     !  length for particle loading
 dpx = xload/ne                     !  particle spacing

!  pseudo-particle charge normalised to give ncrit=1 (rhoc=-1)

 qe = -rho0*dpx

!  pseudo-particle mass (need for kinetic energy diagnostic)

 e_mass = qe / q_over_me

!  set up initial positions

 do i = 1,ne
    x(i) = wall_left + dpx*(i-0.5)
 end do
      
end
!  ======================================
!
!   Set up particle velocity distribution
!
!  ======================================

subroutine loadv

include 'es.h'

 iseed = 76523; idum1 = 137; idum2 = 45126

 do l=1,ne

!  inverted 2v-distribution - amplitude

    vm=vte*(-2.*alog((l-0.5)/ne))**0.5

!  random angle
    rs = ran(iseed)
    theta=2*pi*rs

!  x-component of v
    vx(l)=vm*sin(theta)

 end do


!  scramble particle indicies to remove correlations
!  between x and vx


!        call warn('Rand. particles with clock time')
!	write(15,*) time()
!        idum1 = -time()
!	idum2 = -time()-23
!	idum3 = -time()-47

!  exclude odd one out

      n1=ne
      if (mod(ne,2).ne.0) then
	n1=n1-1
      endif

      do i=1,n1
	j=n1*ran(idum1)+1

	temp1 = vx(i)        ! switch i,j
	vx(i) = vx(j)
	vx(j) = temp1

      end do

end

!  ======================================
!
!    Electron density
!
!  ======================================


subroutine density

include 'es.h'

 re=qe/dx    !  charge weighting factor
     
 rhoe = 0.   ! zero density array

!  gather density onto grid from particles

 do i=1,ne
      xa = x(i)/dx
      j1 = xa
      j2 = j1+1
      f2 = xa-j1
      f1 = 1.-f2
      rhoe(j1) = rhoe(j1) + re*f1
      rhoe(j2) = rhoe(j2) + re*f2
 end do

 if (bc_field.eq.1) then
!   periodic boundaries 
      rhoe(0) = rhoe(0) + rhoe(nx)
      rhoe(nx) = rhoe(0)

 else if (bc_field.eq.2) then
!   reflective - 1st and last (ghost) cells folded back onto physical grid
	iwl=wall_left/dx
	rhoe(iwl+1)=rhoe(iwl+1)+rhoe(iwl)
        rhoe(iwl)=0.	

	iwr=wall_right/dx
	rhoe(iwr)=rhoe(iwr)+rhoe(iwr+1)
        rhoe(iwr+1)=rhoe(iwr)     

 endif

end

!  ============================================
!
!    Direct field integration - slab geometry
!
!  ============================================

subroutine field

include 'es.h'

real :: rhot(0:nxm)     !  net density

!  Add neutral background to get net density
 
 rhot(1:nx) = rhoe(1:nx) + rho0
      

!   integrate div.E=rho directly (trapezium approx)


!   end point - ex=0 mirror at right wall
      Ex(nx+1)=0.
      edc = 0.

!  integrate from right to left

      do j=nx,1,-1
         Ex(j) = Ex(j+1) - 0.5*(rhot(j) + rhot(j+1))*dx
         edc = edc + Ex(j)
      end do


 if (bc_field.eq.1) then

!  periodic fields - remove DC component 
!  -- need this for consistency with charge conservation
      
      do i=1,nx
         ex(i)=ex(i) - edc/nx
      end do
      
!  end points periodic
      ex(0) = ex(nx)

 endif

!   potential - same again
!   - integration const. is arbitrary - phi not used for anything at present

      phi(nx+1)=0.

      do j=nx,1,-1
         phi(j) = phi(j+1) + 0.5*(Ex(j) + Ex(j+1))*dx
      end do 

 
end

!  ============================================
!
!    Diagnostics - graphics and printed output
!
!  ============================================

subroutine diagnostics

include 'es.h'


!  write run information to terminal

 if (mod(itime,iout).eq.0) then
   write (6,*) 'timestep:', itime
 endif

!  do graphics snapshots

 if (mod(itime,igraph).eq.0) then
    call plots
 endif

!  write out time-dep. quantities to file
 if (mod(itime,ihist).eq.0) then
    call histories
 endif

end

!  ===================================
!
!    Electrostatic particle pusher
!
!  ===================================  

subroutine push

include 'es.h'

  do i=1,ne
 
!   interpolate field Ex from grid to particle

      xa = x(i)/dx
      j1 = xa
      j2 = j1+1
      b2 = xa-j1
      b1 = 1.-b2
      exi = b1*Ex(j1) + b2*Ex(j2)
      
!  update velocities

      vx(i) = vx(i) + q_over_me*dt*exi 

 end do



!  update positions (2nd half of leap-frog)

 do i=1,ne
      x(i) = x(i) + dt*vx(i)
 end do

end

!  ===================================
!
!    Particle boundary conditions
!
!  ===================================  

subroutine boundaries

include 'es.h'
     
 iseed1 = 28631        ! random number seed
 iseed2 = 1631         ! random number seed

 xl = grid_length     !  shorthand  
 wl = wall_left
 wr = wall_right


!  loop over all particles to see if any have
!  left simulation region: if so, we put them back again
!  according to the switch 'bc_particle'

 do i=1,ne

!  periodic

      if (bc_particle.eq.1) then

         if (x(i).lt.0.0) x(i) = x(i) + xl
         if (x(i).ge.xl) x(i) = x(i) - xl



! reflective at x = wall_left and x = wall_right

      else if (bc_particle.eq.2) then

         if (x(i).le.wl) then
	    x(i) = 2*wl - x(i)
	    vx(i) = -vx(i)
         endif

	 if (x(i).ge.wr) then
	    x(i) = 2*wr - x(i)
	    vx(i) = -vx(i)
	 endif 



!   reflect at LH boundary
!   re-inject with thermal velocity at RH (high density) boundary

      else if (bc_particle.eq.3) then

         if (x(i).le.wl) then
	    x(i) = 2*wl - x(i) 
	    vx(i) = -vx(i)
	 endif

	 if (x(i).ge.wr) then
	    x(i) = 2*wr - x(i)

!  find new random velocity in cold Maxwellian
	    i_inject = ne*ran(iseed1) + 1
	    vm = vte*(-2.*alog((i_inject-0.5)/ne))**0.5
	    theta = 2*pi*ran(iseed2)
	    vx(i) = -vm*abs(cos(theta))

	  endif
	endif

 end do

end

!  ============================================
!
!    Time-histories
!
!  ============================================

subroutine histories
include 'es.h'



!   kinetic energy
 
 uth = 0.
 do i=1,ne
      uth = uth + 0.5*e_mass*vx(i)**2
 end do


!  write energies out to file

 write(60,'(2f12.6)') itime*dt,uth

 end

!  =========================================
!
!   Graphical snapshots
!
!  =========================================

subroutine plots
include 'es.h'

 real :: xgrid(0:nxm)             ! grid work array

 character*40 cfile
 character*1 cid

 integer :: isnap                ! counts number of calls to routine 'plots'

 data isnap/0/ 
 save isnap                      ! retains value of isnap for next call

 xgrid(0:nx) = (/(i*dx,i=0,nx)/)   ! set up x-values for grid plots


 cid = achar(mod(isnap,10)+48)   !  convert counter to ASCII character '0-9'


!  electron density

 cfile = 'edens'//cid//'.dat'    ! build filename from components
 open (50,file=cfile(1:10))      !  open data file

 write (50,'(2f10.5)') (xgrid(j),-rhoe(j),j=0,nx)
 close(50)


!  electrostatic field

 cfile = 'field'//cid//'.dat'    ! build filename from components
 open (50,file=cfile(1:10))      !  open data file

 write (50,'(2f10.5)') (xgrid(j),Ex(j),j=0,nx)
 close(50)


 isnap = isnap + 1
     
end
!  ===================================
!
!    es.h  Header file for ES PIC code
!
!  ===================================  

parameter(pi=3.141592654,c=3.e8, &
	  nxm=500,               &
	  npm=5000,              &           !  constants
	  nvxm=200,              &
	  ntm=300                & 
	  )

!  Variable declarations which deviate from defaults (integer i-n, otherwise real)

integer :: bc_particle, bc_field


!  ** particle arrays **

common/part/               &         
          x(0:npm), vx(0:npm)

!  ** grid arrays **

common/grid/  & 
          rhoe(0:nxm), Ex(0:nxm), phi(0:nxm), a(0:nxm)


!  ** particle constants/ parameters **

common/phys/           &  
         qe,   vte,   q_over_me,  e_mass, ne,   &
         xload, wall_left, wall_right, rho0, grid_length


!  ** grid variables **

common/grid/            &  
         nt,    nx,  nxo2,   dx,   xn0,   rdx,   dkx,   &
         dt, itime,    yi,   nvx,   nvy,   dvx,   dvy,   vxm,  &
         bc_particle,  bc_field, ntrun


!  ** laser parameters **

common/las/             &  
         a0,    w0, trise,  ilas,  tdel,    tp, xlam

!  ** diagnostic parameters **

common/diag/             &  
       iout,  ihist,   igraph
