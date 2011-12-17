
!     ==================================
!
!     Electrostatic field
!     - FFT Poisson solver
!
!     ==================================

subroutine esf1
  use bopsvars
  implicit none

  integer i, ifail, ipl, imi
  real*8 b,rdx
  real :: rw1(nx), rw2(nx)
  rdx = 1./dx

  do i=1,nx
     rw1(i)=rhoe(i)+rhoi(i)
     rw2(i)=0.0
  end do
  ifail=0
  call c06ecf(rw1,rw2,nx,ifail)
  !     zero dc part
  rw1(1)=0.0
  rw2(1)=0.0
  yrhok(1)=(0.,0.)
  yphik(1)=(0.,0.)

  do ipl=2,nxo2
     imi=nx+2-ipl
     yrhok(ipl)=rw1(ipl)+yi*rw2(ipl)
     yrhok(imi)=rw1(imi)+yi*rw2(imi)
     !     rksq=1.0d0/((ipl-1)*(ipl-1)*dkx*dkx)
     rw1(ipl)=rw1(ipl)*rk2(ipl-1)
     rw2(ipl)=rw2(ipl)*rk2(ipl-1)
     rw1(imi)=rw1(imi)*rk2(ipl-1)
     rw2(imi)=rw2(imi)*rk2(ipl-1)
     yphik(ipl)=rw1(ipl)+yi*rw2(ipl)
     yphik(imi)=rw1(imi)+yi*rw2(imi)
  end do

  yrhok(nxo2+1)=rw1(nxo2+1)+yi*rw2(nxo2+1)
  !     rksq=1.0d0/(nxo2*nxo2*dkx*dkx)
  rw1(nxo2+1)=rw1(nxo2+1)*rk2(nxo2)
  rw2(nxo2+1)=rw2(nxo2+1)*rk2(nxo2)
  yphik(nxo2+1)=rw1(nxo2+1)+yi*rw2(nxo2+1)
  !
  !     ift for potential
  call c06gcf(rw2,nx,ifail)
  call c06ecf(rw1,rw2,nx,ifail)
  call c06gcf(rw2,nx,ifail)
  do i=1,nx
     phi(i)=rw1(i)
  end do
  !
  !     periodic bcs
  if (ifbc.eq.1) then
     phi(nx+1)=phi(1)
     phi(0)=phi(nx)
     !     bounded
  else if (ifbc.eq.2) then
     phi(nx+1)=phi(1)
     phi(0)=phi(nx)
     b=0.5*(phi(nx)-phi(2))*rdx

     do i=0,nx+1
        phi(i)=phi(i)+b*(i-1)*dx
     end do

  endif

  !     es field ex
  do i=1,nx
     ipl=i+1
     imi=i-1
     ex(i)=0.5*(phi(imi)-phi(ipl))*rdx
  end do
  !     bcs
  if (ifbc.eq.1) then
     ex(0)=ex(nx)
     ex(nx+1)=ex(1)
  else if (ifbc.eq.2) then
     ex(nx+1)=0.
  endif
end subroutine esf1
