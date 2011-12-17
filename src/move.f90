
!  ================
!
!  Update positions
!
!  ================

subroutine move(ip1,n,dts)
  use bopsvars
  implicit none

  integer l, n, ip, ip1
  real*8 dts

  ! OpenMP parallel loop

  !$omp parallel do
  !$omp& default(shared)
  !$omp& private(l,ip)

  do l=1,n
     ip=ip1+l-1
     xn(ip)=xo(ip)+dts*ux(ip)/gamma(ip)
  end do

  !$omp end parallel do

end subroutine move
