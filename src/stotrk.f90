

!  ============================================================
!
!     Store data on tracked particles
!
!  ============================================================

subroutine stotrk
  use bopsvars
  implicit none

  integer i, itr

  if (itime.lt.itstart.or.itime.gt.itend) return

  !  Initialise particle tracking
  if (itime.eq.itstart) call trkini
  if (mod(itime,itsk).ne.0) return

  ittrk=ittrk+1

  do i=1,ntrack
     itr = itrack(i)
     uxtrk(i,ittrk)=ux(itr)
     uytrk(i,ittrk)=uy(itr)
     xtrk(i,ittrk)=xn(itr)
     axtrk(i,ittrk)=ax(itr)
     !  rough integration for y-coord
     xytrk(i,ittrk) = xytrk(i,ittrk-1) + dt*itsk*uy(itr)/gamma(itr)
  end do
end subroutine stotrk
