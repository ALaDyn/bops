
!  =================================
!
!  Copy new positions into old array
!
!  =================================

subroutine pshif
  use bopsvars
  implicit none

  integer l

  do l=1,npart
     xo(l)=xn(l)
  end do
end subroutine pshif
