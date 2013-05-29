
!  ===================================
!
!    Particle boundary conditions
!
!  ===================================  

subroutine particle_bc

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
