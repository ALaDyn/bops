
!  ==========================
!
!  Half-shift for diagnostics
!
!  ==========================

subroutine hmove(dtm)
  use bopsvars
  implicit none

  integer ip
  real*8 dtm
  do ip=1,npart
     xn(ip)=xn(ip)+dt*dtm*ux(ip)/gamma(ip)
  end do
  !      call pbcs
end subroutine hmove
