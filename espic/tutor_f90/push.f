
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
