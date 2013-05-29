
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
