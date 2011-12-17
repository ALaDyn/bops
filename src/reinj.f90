

!  ========================================
!
!     Relativistic particle reinjection
!    - conserves  background temp.
!    - pairwise injection for uy,uz to give zero transverse current
!
!  ========================================

      subroutine reinj(vt,n,idir,uxi,uyi,uzi,gami)

      implicit none
      real*8 vt
      integer n, nsamp, idir
      real*8 pi
      parameter(pi=3.14159265,nsamp=16384)
      real*8 usamp(nsamp+1)
      real*8 df,f0,A,u,du,g,uxi,uyi,uzi,gami,g0
      real*8 ute, pio2, umax, rs, theta
      real*8 rano
      integer idum, isamp, i2, nsteps, k, i
      data idum,isamp,i2/-13,1,1/
      save idum,isamp,i2,theta,usamp

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

! Uncomment to debug reinjection sample
!	open(85,file='usample.data')
! 	write (85,'(i6,1pe12.3)') (i,usamp(i),i=1,nsamp)
!	close(85)

      endif

!  Flux in x-dirn: invert Int (vx f(ux) dux) directly

      rs=rano(idum)
      g0 = dmax1(1.d0,1.d0-ute*log(rs))
      uxi = idir*sqrt(g0**2-1.d0)

!  pick random  u from sample:
!  usamp contains integrated momentum flux Int(u f(u) du)
      if (mod(isamp,2).eq.0) then
       i2=min(1.d0*nsamp,nsamp*rano(idum)+1)
!  components for uy,uz
       theta=2*pi*rano(idum)
       uyi = usamp(i2)*cos(theta)
       uzi = usamp(i2)*sin(theta)
      else
!  use -ve of previous sample
       uyi = -usamp(i2)*cos(theta)
       uzi = -usamp(i2)*sin(theta)
      endif

! Uncomment to debug reinjection sampling
! write(*,'(a,i6,a,i6,a,5(1pe12.3))') 'sample',isamp,'index',i2,'u,g:',usamp(i2),uxi,uyi,uzi,g0

      gami = sqrt(1.d0 + uxi**2 + uyi**2 + uzi**2)

!  increment sample index
      isamp = isamp + 1
      if (isamp.gt.nsamp) isamp = 1

      end


