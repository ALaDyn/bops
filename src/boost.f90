
!  =============================
!
!    Relativistic boost parallel
!     to target
!
!  =============================

subroutine boost
  use bopsvars
  implicit none

  integer ip
  real*8 uylab,glab

  !   add drift vel relativistically (from theta0)
  !   ux invariant, uy, gamma from lorentz transfs

  vyav=0.
  do ip=1,npart
     uylab = uy(ip)
     glab = gamma(ip)
     uy(ip) = gam0*( uylab + vy0*glab )
     gamma(ip) = gam0*( glab + vy0*uylab )
     vyav=vyav+uy(ip)/gamma(ip)
  end do

  vyav=vyav/npart
end subroutine boost
