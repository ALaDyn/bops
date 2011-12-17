
!     ==========

subroutine empert
  use bopsvars
  implicit none
  integer i
  w0=sqrt(1.0+wp*wp)
  do i=1,nx+1
     ey(i)=-a0*sin(i*dx-dx)
     bz(i)=ey(i)/w0
     ff(i)=0.5*(ey(i)+bz(i))
     fb(i)=0.5*(ey(i)-bz(i))
  end do
end subroutine empert
