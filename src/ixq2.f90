
!  ========================================================
!
!       Quiet start for rel. velocity load
!
!  ========================================================

subroutine ixq2(i1,n)
  use bopsvars
  implicit none

  integer l, ip, i1, ipn, n, no2
  real*8 rsi, rs

  if (n.eq.0) return

  !   binary reversal for indices
  !   swap 4 particles at a time (in indep. vel. quadrants)
  rs=0.
  no2=n/2
  do l=0,no2-1
     ip=2*l+i1
     ipn=2*rs*no2+i1
     xn(ipn)=xo(ip)
     xn(ipn+1)=xo(ip+1)
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
end subroutine ixq2
