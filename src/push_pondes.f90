!     ==============================================
!
!     2v particle pusher
!
!     TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
!
!     ==============================================


subroutine push_pondes(ip1,n,dts)
  use bopsvars
  implicit none
  integer, intent(in) :: ip1  ! start index
  integer, intent(in) :: n    ! # particles of this type
  real*8, intent(in) :: dts   ! timestep

  real*8 :: beta, xa, b1, b2, rdx
  real*8 :: exi, eyi, ezi, byi, bzi, uzi, ayi, azi, epi, uyi
  integer :: l, ip, i1, i2
  real*8 :: gam2,gam1,tt,ss,uxd,uyp,uzp,uxp,uxm,uym,uzm,uzn1

  rdx = 1./dx

  do l=1,n
     ip=ip1+l-1
     beta=q(ip)/m(ip)*dts*0.5  ! charge/mass constant

     !     interpolate fields
     ! linear weighting scheme W_j=1-|x_i-x_j| 

     xa=xo(ip)*rdx
     i1=int(xa)+1
     i2=i1+1
     b2=i1-xa
     b1=1.-b2
     if (ip.le.ne) then
        !  include fpond on electrons
        uyi = b1*ay(i1) + b2*ay(i2) + gam0*vy0  ! include boost momentum
        uzi = b1*az(i1) + b2*az(i2)
        gam1 = sqrt(1. + ux(ip)**2 + uyi**2 + uzi**2)
        ! corrected pond field including total particle gamma
        epi = (b1*epond(i1) + b2*epond(i2))/gam1
        exi = b1*ex(i1) + b2*ex(i2) + epi  ! total ex at particle

     else
        exi = b1*ex(i1) + b2*ex(i2)
        eyi = 0.
        ezi = 0.
        uyi = gam0*vy0
        uzi = 0.
     endif

     !     ES accn
     ux(ip) = ux(ip) + beta*exi

     uy(ip) = uyi  ! Transverse momenta completely determined by the fields
     uz(ip) = uzi
     gamma(ip)=sqrt(1.0+ux(ip)**2+uz(ip)**2+uy(ip)**2)
  end do
end subroutine push_pondes





