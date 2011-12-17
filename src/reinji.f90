
!  ========================================
!
!     Relativistic particle reinjection - IONS
!    - conserves  background temp.
!
!  ========================================

      subroutine reinji(vt,n,idir,uxi,uyi,uzi,gami)

      implicit none
      real*8 pi, rs
      integer nsamp, idum, isamp, nsteps, i, k, i2, idir, n
      parameter(pi=3.14159265,nsamp=16384)
      real*8 usamp(nsamp)
      real*8 df,f0,A,u,du,g,uxi,uyi,uzi,gami,g0
      real*8 ute, vt, pio2, umax, theta
      real*8 rano
      data idum,isamp/-11,1/
      save idum,isamp,usamp

!  conserve flux in x dirn: maintain const. Te

!  Tabulate momenta on first call
        ute = vt**2
        pio2 = pi/2
      if (isamp.eq.1) then

        nsteps=100000
!  dimensionless electron/ion temp
        umax=6.*sqrt(ute + 4*ute**2)
        du=umax/nsteps
!  distn norm factor
        A = 1.0001*nsamp/2./pi/ute/(1+ute)
        f0=0.
        k=1

!  load loop
        do i=1,nsteps

!  cumulative integral of rel. distn function f(u)
        u = i*du
          g = sqrt(u**2 + 1.d0)
        df = 2*pi*A*du*u*exp((1.-g)/ute)

        f0 = f0+df

        if(f0.ge.k .and. k.le.nsamp) then
!  store u corresponding to integer values of f0
            usamp(k) = u
          k=k+1
        endif
        end do
      endif

!  Flux in x-dirn: invert Int (vx f(ux) dux) directly

      rs=rano(idum)
      g0 = max(1.d0,1.d0-ute*log(rs))
      uxi = idir*sqrt(g0**2-1.d0)

!  pick random  u from sample:
!  usamp contains integrated momentum flux Int(u f(u) du)

      i2=min(1.d0*nsamp,nsamp*rano(idum)+1)
      theta=2*pi*rano(idum)
!  components uy,uz
      uyi = usamp(i2)*sin(theta)
      uzi = usamp(i2)*sin(theta)
      gami = sqrt(1. + uxi**2 + uyi**2 + uzi**2)

!  increment sample index
      isamp = isamp + 1
!      if (isamp.gt.nsamp) isamp = 1

      end





