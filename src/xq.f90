!
!     ==========
!
subroutine xq(i1,n,xstart,xlo)
  use bopsvars
  implicit none

  real*8 :: xstart, xlo
  integer :: n, l, ip, i1
  real*8 :: dpx, rs, rsi

  if (n.eq.0) return
  !
  !   binary reversal for positions
  dpx=xlo/n
  rs=0.
  do l=1,n
     ip=l+i1-1
     xo(ip)=dpx/2.+rs*xlo+xstart
     rsi=1.0
10   rsi=rsi/2.
     rs=rs-rsi
     if (rs.ge.0.) goto 10
     rs=rs+2*rsi
  end do

end subroutine xq
