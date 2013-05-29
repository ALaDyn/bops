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
    call particle_bc		! particle boundary conditions

    call density		! compute density

    call field			! compute electric field (Poisson)

    call diagnostics		! output snapshots and time-histories

 end do

!  close time-hist files
 close(60)
end
