
!  ===================================================
!
!   Read in particle and field data from restart file
!
!  ===================================================

subroutine restrt
  use bopsvars
  implicit none

  integer i

! TODO: needs additional particle, field variables 

  read (9,*,end=99) npart,ne,ni,ion1,itim0
  read (9,*,end=99) (xo(i),xn(i),ux(i),uy(i),gamma(i),i=1,npart)
  read (9,*,end=99) nx,(rhoi(i),ex(i),ey(i),bz(i),i=1,nx)
  !  recover forward and backward waves
  do i=1,nx+1
     ff(i)=0.5*(ey(i)+bz(i))
     fb(i)=0.5*(ey(i)-bz(i))
  end do
  return
99 call warn('error in restart data         ')
  stop
end subroutine restrt
