!  ==============================================
!
!     3V pusher:  (Ex,Ey,Ez), (0,By,Bz)
!
!  ==============================================

subroutine pushc(ip1,n,dts)
  use bopsvars
  implicit none
  integer, intent(in) :: ip1  ! start index
  integer, intent(in) :: n    ! # particles of this type
  real*8, intent(in) :: dts   ! timestep

  real*8 :: dtso2, beta, xa, b1, b2, rdx
  real*8 :: exi, eyi, ezi, byi, bzi, tt1
  integer :: l, ip, i1, i2
  real*8 :: gam2,gam1,ty,tz,sy,sz,uxd,uyp,uxp,uzp,uxm,uym,uzm

  dtso2=dts/2.
  rdx = 1./dx

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
     ezi=b1*ez(i1)+b2*ez(i2)
     byi=b1*by(i1)+b2*by(i2)
     bzi=b1*bz(i1)+b2*bz(i2)

     !   first half-accn
     uxm = ux(ip) + beta*exi
     uym = uy(ip) + beta*eyi
     uzm = uz(ip) + beta*ezi

     !   rotation
     gam1=sqrt(1.0 + uxm**2 + uym**2 + uzm**2)
     ty = beta*byi/gam1
     tz = beta*bzi/gam1
     tt1 = 1.0 + ty**2+tz**2
     sy = 2.0*ty/tt1
     sz = 2.0*tz/tt1

     uxd = uxm + uym*tz - uzm*ty
     uyp = uym - uxd*sz
     uzp = uzm + uxd*sy
     uxp = uxd + uyp*tz - uzp*ty

     !   second half-accn
     ux(ip) = uxp + beta*exi
     uy(ip) = uyp + beta*eyi
     uz(ip) = uzp + beta*ezi

     gamma(ip) = sqrt(1.0 + ux(ip)**2 + uy(ip)**2 + uz(ip)**2)
  end do
end subroutine pushc








