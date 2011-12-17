
! =============================================================
!
!            FILTER1
!
!    3-point filter - returns smoothed function in same array
!
! =============================================================


      subroutine filter1(w1,n)
      implicit none
      integer, intent(in) :: n
      real*8, intent(inout) :: w1(n)
      real*8 :: w2(n) ! work array
      integer i

      do i=2,n-1
        w2(i)=0.25*(w1(i-1)+w1(i+1))+0.5*w1(i)
      end do

!  end points
      w2(1)=0.75*w1(1)+0.25*w1(2)
      w2(n)=0.75*w1(n)+0.25*w1(n-1)

      do i=1,n
        w1(i)=w2(i)
      end do
      end
