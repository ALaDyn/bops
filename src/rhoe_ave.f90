
!  ===========================
!
!   Time-averaged electron density
!
!  ===========================

subroutine rhoeav
  use bopsvars
  implicit none

  integer isnap, i, iplas

!  real*8 :: rhoe_last(0:nx+1)
  real*8 :: tcycav, rhel


  if (itime.eq.1) then
    dcrhe(1:nx+1) = (rhoi(1:nx+1) -vy0*jyi(1:nx+1))/gam0  ! Initially set to ion density
    rhoe_last=0.
  endif


  if (mod(itime,itav).eq.0) then
!  Store in permanent array and reset for next cycle
     dcrhe=rhoe_last
     rhoe_last = 0.

  else 
! Add latest density to cumulative DC electron number density (+ve)
    do i=1,nx
      rhel = -(rhoe(i) - vy0*jye(i))/gam0
      rhoe_last(i) = rhoe_last(i) + rhel/itav
    end do

  endif

iplas = 9./dx
!write(*,'(i6,4f12.3)') itime,rhoe_last(iplas),dcrhe(iplas),rhoe(iplas),rhoi(iplas)

end subroutine rhoeav










