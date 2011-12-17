!  ==============================================
!
!     2v particle pusher with fluid uz
!
!  TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
!
!  ==============================================


subroutine push_pond(ip1,n,dts)
  use bopsvars
  implicit none
  integer, intent(in) :: ip1  ! start index
  integer, intent(in) :: n    ! # particles of this type
  real*8, intent(in) :: dts   ! timestep

  real*8 :: qom, beta, xa, b1, b2, rdx
  real*8 :: exi, eyi, ezi, byi, bzi, uzi
  integer :: l, ip, i1, i2
  real*8 :: gam2,gam1,tt,ss,uxd,uyp,uzp,uxp,uxm,uym,uzm,uzn1

  rdx = 1./dx

  do l=1,n
     ip=ip1+l-1
     qom = q(ip)/m(ip)
     beta=qom*dts*0.5  ! charge/mass constant

     !   interpolate fields to particle positions

     xa=xo(ip)*rdx
     i1=int(xa)+1
     i2=i1+1
     b2=i1-xa
     b1=1.-b2
     exi = b1*ex(i1) + b2*ex(i2)
   
     ! transverse momentum from pz=az including thermal motion from uy
     uzi = -qom*b1*az(i1) -qom* b2*az(i2) + uy(ip)
     byi = b1*by(i1) + b2*by(i2)
     ezi = b1*ez(i1) + b2*ez(i2)


     !   half-accn

     uxm = ux(ip) + beta*exi
     uzm = uzi + beta*ezi   ! intermediate uz
     uym = uy(ip)

     !   ux,uz rotation

     gam1 = sqrt(1.0+uxm**2+uzm**2+uym**2)
     tt = beta*byi/gam1
     ss = 2.d0*tt/(1.0+tt**2)

     uxd = uxm - uzm*tt
     uzp = uzm + uxd*ss
     uxp = uxd - uzp*tt

     !   half-accn

     ux(ip)=uxp+beta*exi
     uz(ip)=uzp+beta*ezi
     !        uz(ip)=uzi
     gamma(ip)=sqrt(1.0+ux(ip)**2+uz(ip)**2+uy(ip)**2)
  end do
end subroutine push_pond








