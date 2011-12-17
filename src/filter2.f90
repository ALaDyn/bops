
! =============================================================
!
!            FILTER2
!
!    3-point filter - returns smoothed function in 2nd array
!
! =============================================================

      subroutine filter2(w1,w2,nmax,n)
      implicit none
      integer nmax,n
      real*8 w1(0:nmax),w2(0:nmax)
      integer i

      do i=2,n-1
        w2(i)=0.25*(w1(i-1)+w1(i+1))+0.5*w1(i)
      end do

      w2(1)=0.75*w1(1)+0.25*w1(2)
      w2(n)=0.75*w1(n)+0.25*w1(n-1)
      end
