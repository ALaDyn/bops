
!  =============================
!
!     2v relativistic maxwellian
!     velocity loading
!
!  =============================

subroutine rel2u(i1,n,vt)
  use bopsvars
  implicit none

  real*8 ute, vt, pio2, a, theta, theta2, theta3, theta4
  real*8 rano, mass, uth0
  integer idum, n, nsteps, k, ip, i1, i, l
  real*8 f0,df,rs,u,du,umax,finf,g,aux,auy
  data idum/-9/
  save idum

  if (n.eq.0) return
  nsteps=100000
  !  dimensionless electron temp
  ute = vt**2
  pio2 = pi/2
  umax=6.*sqrt(ute + 4*ute**2)
  du=umax/nsteps
  aux=0.d0
  auy=0.d0
  !  distn norm factor
  A = n/2/pi/ute/(1+ute)
  f0=0.
  k=1
  ip=i1

  !  load loop
  do i=1,nsteps
     !  cumulative integral of rel. distn function f(u)
     u = (i-0.5d0)*du
     g = sqrt(u**2 + 1.d0)
     df = pio2*A*du*u*exp(max(-10.0d0,-(g-1.)/ute))
     f0 = f0+df

     if(f0.ge.k) then
        !  assign 4 particles in each quadrant to avoid drifts
        theta = pio2*rano(idum)
        theta2 = theta+pio2
        theta3 = theta+pi
        theta4 = theta+3*pi/2
        !  1st quadrant
        ux(ip) = u*cos(theta)
        uy(ip) = u*sin(theta)
        !  2nd quadrant
        ux(ip+1) = u*cos(theta2)
        uy(ip+1) = u*sin(theta2)
        !  3rd quadrant
        ux(ip+2) = u*cos(theta3)
        uy(ip+2) = u*sin(theta3)
        !  4th quadrant
        ux(ip+3) = u*cos(theta4)
        uy(ip+3) = u*sin(theta4)

        aux = aux+ux(ip)+ux(ip+1)+ux(ip+2)+ux(ip+3)
        auy = auy+uy(ip)+uy(ip+1)+uy(ip+2)+uy(ip+3)

        k=k+1
        ip=ip+4
     endif
  end do

  !  correct drifts to ensure 1st moments vanish
  aux = aux/n
  auy = auy/n

  do l=1,n
     ip=i1+l-1
     ux(ip)=ux(ip)-aux
     uy(ip)=uy(ip)-auy
  end do


  !  calculate relativistic gamma factor

  if (i1.eq.1) then
     mass = me
  else
     mass = mi
  endif

  do l=1,n
     ip=i1+l-1
     gamma(ip) = sqrt(1.+ ux(ip)**2 + uy(ip)**2)
     uth0=uth0+mass*(gamma(ip)-1.0)
  end do

end subroutine rel2u
