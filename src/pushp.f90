
!  ==============================================
!
!     2v particle pusher
!
!  TE fields only (p-pol):  (Ex,Ey,0)  (0,0,Bz)
!
!  ==============================================


subroutine pushp(ip1,n,dts)
  use bopsvars
  implicit none

  integer, intent(in) :: ip1  ! start index
  integer, intent(in) :: n    ! # particles of this type
  real*8, intent(in) :: dts   ! timestep

  real*8 :: beta, xa, b1, b2, exi, eyi, bzi
  integer :: l, ip, i1, i2
  real*8 :: gam2,gam1,tt,ss,uxd,uyp,uxp,uxm,uym,uzm
  real*8 :: rdx

  rdx=1./dx

!  beta=qom*dts*0.5


  do l=1,n
     ip=ip1+l-1
     beta=q(ip)/m(ip)*dts*0.5  ! charge/mass constant

     !   interpolate fields

     xa=xo(ip)*rdx
     i1=int(xa)+1
     i2=i1+1
     b2=i1-xa
     b1=1.-b2
     exi=b1*ex(i1)+b2*ex(i2)
     eyi=b1*ey(i1)+b2*ey(i2)
     bzi=b1*bz(i1)+b2*bz(i2)


     !   half-accn

     uxm = ux(ip) + beta*exi
     uym = uy(ip) + beta*eyi
     uzm = uz(ip)

     !   rotation

     gam1=sqrt(1.0+uxm**2+uym**2+uzm**2)
     tt=beta*bzi/gam1
     ss=2.d0*tt/(1.d0+tt**2)

     uxd = uxm + uym*tt
     uyp = uym - uxd*ss
     uxp = uxd + uyp*tt

     !   half-accn

     ux(ip)=uxp+beta*exi
     uy(ip)=uyp+beta*eyi

  end do


  !   get gamma

  do l=1,n
     ip=ip1+l-1
     gamma(ip)=sqrt(1.0+ux(ip)**2+uy(ip)**2+uz(ip)**2)
  end do

end subroutine pushp









