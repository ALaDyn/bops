
!  ========================================================
!
!       Binary quiet start
!
!  ========================================================

subroutine ixq(i1,n)
  use bopsvars
  implicit none

  integer i1, n, l, ip, ipn
  real*8 rs, rsi

  if (n.eq.0) return

  !   binary reversal for indices
  rs=0.
  do l=1,n
     ip=l+i1-1
     ipn=rs*n+i1
     xn(ipn)=xo(ip)
     rsi=1.0
10   rsi=rsi/2.
     rs=rs-rsi
     if (rs.ge.0.) goto 10
     rs=rs+2*rsi
  end do

  do l=1,n
     ip=l+i1-1
     xo(ip)=xn(ip)
  end do
end subroutine ixq
