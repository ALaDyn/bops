
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

!   potential energy
 
 upot = 0.
 do i=1,nx
      upot = upot + 0.5*Ex(i)**2*dx
 end do
 
! Total
 utot = upot + uth

!  write energies out to file

 write(60,'(4f12.6)') itime*dt,uth,upot,utot

 end
