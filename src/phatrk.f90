

!  ==========================================================
!
!        Catch exit phase of particles
!
!  ==========================================================

subroutine phatrk
  use bopsvars
  implicit none

  integer i, l
  real*8 gam

  do i=1,ntrack
     l=itrack(i)
     !  store exit phase of particle
     if (phat(i).eq.0.and.xn(l).lt.xsol) then
        phat(i)=mod(dt*itime,2*pi)
        !  ... then catch energy before it exits!
     else if (uest(i).eq.0.and.xn(l).ge.xl-5*dx.and.ux(l).gt.0 &
          .and.phat(i).gt.0) then

        gam = gam0*(gamma(l) - vy0*uy(l))
        uest(i)=gam-1.
     endif
  end do

end subroutine phatrk



