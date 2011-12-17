!  =======================
!  SLAB - set up flat slab
!  ======================

subroutine slab(n,x,x_start,x_finish)
  implicit none

  real*8 :: x(n)  ! Array returning particle positions
  real*8, intent(in) :: x_start, x_finish
  integer, intent(in) :: n  ! # particles
  real*8 :: dpx
  integer :: i

  if (n.eq.0) return

! set up flat slab
   dpx = (x_finish - x_start)/n
   do i=1,n
      x(i) = x_start + i*dpx - dpx/2.
   end do

end subroutine slab
