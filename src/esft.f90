

!     ==========

subroutine esft
  use bopsvars
  implicit none

  integer i, ifail, nk, k, igk, j
  real*8 uekm, xmin, xmax, dax, ukmax, ukmin
  real ze
  real :: rw1(nx), rw2(6*nx+150)  ! FT work arrays
  integer :: iwk(6*nx+150)    !

  real*8 uesk(0:nx+1)
  complex cx(nx)

  !

  do i=1,nx
     rw1(i)=ex(i)
  end do

  ifail=0

  call fftrc(rw1,nx,cx,iwk,rw2)

  !  get remainder of coefficients
  nk=nx/2

  do k=2,nx/2
     cx(nx+2-k) = conjg(cx(k))
  end do

  uekm = 0.
  !  Power density
  do i=2,nk
     uesk(i)=(abs(cx(i))**2 + abs(cx(nx+2-i))**2)/nk**2
     uekm=max(uekm,uesk(i))
  end do


  uekm=0.
  uesk(0)=0.
  xmin = 0.
  xmax = pi/dx
  dax = dkx

  igk=1+nk/2000

  ukmax=ze(uekm)+1
  ukmin=ukmax-6

  do j=1,nk
     uesk(j)=max(uesk(j),ukmin)
  end do

  !  plot ft Ek**2 on (manual) log-lin
!#### Anupam & Bin 2009/2010: conforms to this version change 730 to 7300
  call grxy(xk,uesk,nk,7300+idc,igk,-2 & 
       ,'       k       ','     Ues       ','uesk'//ctime(1:12) &
       ,xmin,xmax,dax,ukmin,ukmax,1.)
end subroutine esft









