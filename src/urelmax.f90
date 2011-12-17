
!     =============================
!
!     2v relativistic maxwellian velocity loading
!     - loads one vel. component at a time
!     =============================

subroutine urelmax(uu,n,vt)
  use bopsvars
  implicit none

  integer, intent(in) :: n
  real*8 :: uu(n)

  integer :: idum, nsteps, ip, k, i
  real*8 vt, ute, aux, auy, pio2, a, theta
  real*8 rano

  real*8 f0,df,rs,u,du,umax,finf,g
  data idum/-9/
  save idum
  integer istep

  if (n.eq.0) return
  ! cold particles
  if (vt.eq.0) then
     do i=1,n
        uu(i)=0.
     end do
    return
  endif
  nsteps=n*20

  !     dimensionless electron temp
  ute = vt**2
  umax=6.*sqrt(ute + 4*ute**2)
  du=umax/nsteps
  aux=0.d0
  auy=0.d0
  pio2=0.5*pi

  !     distn norm factor
  A = n/2/pi/ute/(1+ute)
  f0=0.
  k=1
  ip=1
  istep=1
  !     first step of integral of rel. distn function f(u)
  u = (istep-0.50)*du
  g = sqrt(u**2 + 1.0)
  df = 2*pi*A*du*u*exp(max(-10.d0,-(g-1.)/ute))
  f0 = f0+df

  !     load loop
  do k=1,n
     do while (f0.lt.k)
        istep=istep+1
        !     cumulative integral of rel. distn function f(u)
        u = (istep-0.50)*du
        g = sqrt(u**2 + 1.0)
        df = 2*pi*A*du*u*exp(max(-10.d0,-(g-1.)/ute))
        f0 = f0+df
!	write (*,*) k, istep, f0
     end do
     !     store velocity
     theta = pio2*rano(idum)
     uu(ip) = u*cos(theta)
     ip=ip+1
  end do


end subroutine urelmax



