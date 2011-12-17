
!  =================
!
!   Waveguide mode
!
!  =================

subroutine wguide
  use bopsvars
  implicit none

  integer i, ireset
  logical lreset
  lreset=.false.
  do i=1,10
     ireset=xm1/dt*(2*i+1)
     if (mod(itime,ireset).eq.0) lreset=.true.
  end do

  !  reset plasma to initial condition: fresh wall
  if (lreset) then
     call blank
     write(15,*) 'PLASMA RELOAD at t=',itime*dt,itime,' dt'
     call parload
     call movie
  endif
end subroutine wguide
