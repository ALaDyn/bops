
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
