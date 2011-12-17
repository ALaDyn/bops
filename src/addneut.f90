
!     ==========
  !   add neutralising background

subroutine addneut
  use bopsvars
  implicit none
  integer :: i,i1,i2,i3
  real*8 :: xplas,xc,xprime,a,xlol,x

if (inprof.eq.6) then
     !  tanh ramp

     xplas=xsol-xm1
     xc = xm1+xplas/2.
     xlol = 2*pi*xlolam/gam0
     a = 2/xlol/rho0*gam0**3
     xprime=x-xc
     i1=xm1/dx+1
     i2=xsol/dx+1
     i3=xm2/dx+1
     rhoi(1:i1-1)=0.
     do i=i1,i2
       xprime=i*dx-xc
       rhoi(i)=nonc/2.*(1.+tanh(a*xprime))
     end do
     rhoi(i2:i3)=nonc
     rhoi(i3:nx+1)=0.
else

!  default - set ion density equal to electron

  do i=1,nx+1
     rhoi(i)=-rhoe(i)
  end do

endif

end subroutine addneut
