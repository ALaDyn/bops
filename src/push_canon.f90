!     ==============================================
!
!     1v particle pusher
!	exploiting cons. canonical momentum p_perp = A_perp
!     
!
!     ==============================================


subroutine push_canon(ip1,n,dts)
  use bopsvars
  implicit none
  integer, intent(in) :: ip1  ! start index
  integer, intent(in) :: n    ! # particles of this type
  real*8, intent(in) :: dts   ! timestep

  real*8 beta, xa, b1, b2, rdx
  real*8 exi, eyi, ezi, byi, bzi, uzi, ayi, azi, epi
  integer l, ip, i1, i2
  real*8 gam2,gam1,tt,ss,uxd,uyp,uzp,uxp,uxm,uym,uzm,uzn1

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
        !  interpolate transverse vector potentials to particle
        ayi = b1*ay(i1) + b2*ay(i2)
        azi = b1*az(i1) + b2*az(i2)
! electrostatic field at particle
        exi = b1*ex(i1) + b2*ex(i2) + epi  

     else
        exi = b1*ex(i1) + b2*ex(i2)
! ignore transverse field on ions
        ayi = 0.
        azi = 0.
     endif

     !     ES accn

     ux(ip) = ux(ip) + beta*exi
     uy(ip) = ayi + gam0*vy0   ! include boost momentum on uy
     uz(ip) = azi
     gamma(ip)=sqrt(1.0+ux(ip)**2+uz(ip)**2+uy(ip)**2)
  end do
end subroutine push_canon





