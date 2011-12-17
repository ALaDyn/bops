
!     ==========

subroutine bunif
  use bopsvars
  implicit none
  integer i
  do  i=1,nx+1
     bz(i)=a0
  end do
end subroutine bunif
