
!  ==================================
!
!    Electromagnetic fields
!
!      - Dawson/Langdon algorithm
!
!   23/8/95:  Inclusion of TM fields Ez, By
!
!  ==================================

subroutine emfield
  use bopsvars
  implicit none

  real*8 dto4
  integer i, ib
  real*8 ffo(0:nx+1),gfo(0:nx+1), gam(1:nx+1)

  dto4=0.25*dt

  !  store old forward wave values
  do i=1,nx+1
     ffo(i)=ff(i)
     gfo(i)=gf(i)
  end do

  do i=nx+1,2,-1
     ib=nx-i+2

     !   forward TE wave
     ff(i)=ff(i-1)-dto4*(jm(i-1)+jp(i))
     !   backward TE wave
     fb(ib)=fb(ib+1)-dto4*(jm(ib+1)+jp(ib))

     !   forward TM wave
     gf(i)=gf(i-1)-dto4*(jzm(i-1)+jzp(i))
     !   backward TM wave
     gb(ib)=gb(ib+1)-dto4*(jzm(ib+1)+jzp(ib))
  end do

  !   open bcs
  if (iembc.le.2) then
     fb(nx+1)=0.
     ff(1)=0.
     gb(nx+1)=0.
     gf(1)=0.
     ay(nx+1) = 0.
     az(nx+1) = 0.

  else if (iembc.eq.3) then
     !  reflect wave at RHB
     fb(nx+1)=ff(nx+1)
     gb(nx+1)=gf(nx+1)

  else if (iembc.eq.4) then
     !  reflect wave at LHB
     ff(1)=fb(1)
     gf(1)=gb(1)
  endif

  !   wave launched from x=0
  call laser

  !   fields
  do i=1,nx+1

     !  TE
     ey(i)=ff(i)+fb(i)
     bz(i)=ff(i)-fb(i)

     !  TM - note sign of By!!
     ez(i)=gf(i)+gb(i)
     by(i)=gb(i)-gf(i)
  end do

  !  vector potential: backward sweep
  do i=1,nx
     ay(i) = ay(i+1) - dt*( ffo(i+1) + ff(i) )
     az(i) = az(i+1) - dt*( gfo(i+1) + gf(i) )
  end do

  !  ponderomotive force from vector potential
  ay(0) = ay(1) - .5*dt*( ffo(1) )
  az(0) = az(1) - .5*dt*( gfo(1) )

  do i=1,nx
     gam(i) = sqrt(1.+ay(i)**2 + az(i)**2)
     !  leave off gamma factor - this gets corrected in pusher because
     !  we need to include ux too.  Ep = 1/2 d/dx A^2
     epond(i) = .25/dx*( (ay(i+1)**2-ay(i-1)**2) &
          + (az(i+1)**2-az(i-1)**2) )
  end do

end subroutine emfield
