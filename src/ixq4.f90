
!  ========================================================
!
!       Quiet start for rel. velocity load
!
!  ========================================================

subroutine ixq4(i1,n)
  use bopsvars
  implicit none

  integer n, no4, l, ip, i1, ipn
  real*8 rs, rsi

  if (n.eq.0) return

  !   binary reversal for indices
  !   swap 4 particles at a time (in indep. vel. quadrants)
  rs=0.
  no4=n/4
  do l=0,no4-1
     ip=4*l+i1
     ipn=4*rs*no4+i1
     xn(ipn)=xo(ip)
     xn(ipn+1)=xo(ip+1)
     xn(ipn+2)=xo(ip+2)
     xn(ipn+3)=xo(ip+3)
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
end subroutine ixq4
