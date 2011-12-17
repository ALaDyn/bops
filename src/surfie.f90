
!     ==========

      subroutine surfie
      use bopsvars
      implicit none

      integer isurf,i
      real*8 :: rhoimax, rhoemax
      isurf=max(xm1,dx)/dx
      exsurf(ntc)=ex(isurf)
      eysurf(ntc)=ey(isurf)
      bzsurf(ntc)=bz(isurf)

      rhoimax=0.
      rhoemax=0.
      do i=1,nx
        !  lab frame electron and ion densities
        rhoemax = max(rhoemax,-(rhoe(i) - vy0*jye(i))/gam0)
        rhoimax = max(rhoimax,(rhoi(i) - vy0*jyi(i))/gam0)
      end do
      idenmax(ntc) = rhoimax
      edenmax(ntc) = rhoemax

      end
