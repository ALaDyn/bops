!
!     ==========
!
subroutine rand(i1,n,vt)
  use bopsvars
  implicit none

  integer idum, l, n, ip, i1
  real*8 rano, vm, vt, theta
  data idum/-7/
  save idum
  do l=1,n
     ip=l+i1-1
     xo(ip)=rano(idum)*xl
     vm=vt*(-2.*alog((l-0.5)/n))**0.5
     theta=2*pi*rano(idum)
     ux(ip)=vm*cos(theta)
     uy(ip)=vm*sin(theta)
  end do
end subroutine rand
