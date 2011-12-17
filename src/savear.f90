

!  ==========================================
!
!   Dump particle and field data for restart
!
!  ==========================================

subroutine savear
  use bopsvars
  implicit none

  integer i
! TODO: add new particle, field variables
  open(8,file='particle_dump')
  write(*,*) 'Doing particle dump at timestep',itime
  write (8,101) npart,ne,ni,np,itime
101 format(4i8)
  write (8,102) (xn(i),ux(i),uy(i),q(i),m(i),i=1,npart)
102 format(5(1pe16.8))
  write (8,103) nx,(rhoi(i),ex(i),ey(i),bz(i),i=1,nx)
  close(8)

103 format(i8/4(1pe16.8))
end subroutine savear
