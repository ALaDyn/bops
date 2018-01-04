

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

!  Select lab or sim frame with ioboost
  do i=1,ntrack
     itr = itrack(i)
     if (ioboost.eq.1) then
       uxtrk(i,ittrk)=ux(itr)
       uytrk(i,ittrk)=uy(itr)
       uztrk(i,ittrk)=uz(itr)
       xtrk(i,ittrk)=xn(itr)
       axtrk(i,ittrk)=ax(itr)
     !  rough integration for y-coord
       ytrk(i,ittrk) = ytrk(i,ittrk-1) + dt*itsk*uy(itr)/gamma(itr)
     else
       uxtrk(i,ittrk)=ux(itr)
       uztrk(i,ittrk)=uz(itr)
       uytrk(i,ittrk) = gam0*( uy(itr) - vy0*gamma(itr) )
       xtrk(i,ittrk) = gam0*xn(itr)/xconv
       axtrk(i,ittrk)=ax(itr) ! TODO: need lab frame force 
     !  rough integration for y-coord
       ytrk(i,ittrk) = ytrk(i,ittrk-1) + dt*itsk*uytrk(i,ittrk)/gamma(itr)
     endif
  end do
end subroutine stotrk
