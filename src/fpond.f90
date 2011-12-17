
!     ==================================
!
!     Electromagnetic fields
!
!     - ponderomotive laser model for step-profile
!
!
!     ==================================

!     fpond(t,tpulse,sigma0,vosc,omega,x,y,z,epon_x,epon_y,epon_z,phipon)

      subroutine fpond
      use bopsvars
      implicit none

      integer i
      real*8 pha, gamma_c, wp_r, xc, chi_vac, ezhelm, byhelm
      real*8 t                    ! time
      real*8 tptot               ! pulse duration
      real*8 lskin
!     real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surfa

      real*8  phipon, epon_x, epon_y, epon_z ! pond. potential and field

      real*8  xf, yf, zf, Tpon, Xpon, gamma_pon, phi_m
      real*8 dxlab, xcold, nonc_eff, rhoup
      integer nover, icrit

!     linear rise

      t=dt*itime
      pha=w0*t
      Tpon = min(1.d0,t/trise)
!  Find critical surface
      dxlab = dx
      xcold = xcrit
      call findxc(xcold,xcrit,rhoi,nx,dxlab,rhotrak)

!  Determine mean ion density beyond critical surface (5 skin depths)
!  - move to findxc.f

      icrit = xcrit/dx
      nover = 5./dx
      rhoup = 0.

      do i=icrit+1,icrit+nover
         rhoup = rhoup+rhoi(i)
      end do

      nonc_eff = rhoup/nover   ! Effective n/nc seen by EM wave
      if (mod(itime,iout).eq.0) &
      write (*,'(a20,2f9.2)') 'n/nc_eff, xc: ',nonc_eff, xcrit

!     Standing wave in vacuum; evanescent in overdense plasma
!     Use standard solution of Helmholtz equation for step profile

      gamma_c = sqrt(1.+4*a0**2/nonc_eff) ! gamma at surface
      wp_r = sqrt(nonc_eff/gamma_c)  ! effective plasma frequency
      lskin = 1./wp_r   ! Rel. skin depth in EM units

!   Phase factor given by tan(phi) = -k0 * l_s = k0 * c/wp
      phi_m = atan(-w0/wp_r)

      do i=1,nx
         xc = (i-1)*dx - xcrit  ! adjust for fpond bunching

         if (xc.le.0) then
!     vacuum solution
            chi_vac = xc + phi_m      ! Vacuum phase
            xf = sin(2*chi_vac)       ! epond_x

            ezhelm = sin(chi_vac)     ! laser 'Ez'
            byhelm = cos(chi_vac)     ! laser 'Bz'

         else
!   evanescent wave
            xf = -2/lskin*sin(phi_m)**2*exp(-2*xc/lskin) ! epond_x

            ezhelm = sin(phi_m)*exp(-xc/lskin) ! laser Ez inside
            byhelm = -sin(phi_m)/lskin*exp(-xc/lskin) ! laser Ez inside

         endif

         phipon = (2*a0*ezhelm)**2
         gamma_pon= sqrt(1.+ phipon*Tpon*sin(pha)**2) ! relativistic phi
         Epon_x = 2*a0**2*Tpon*xf
         epond(i) = Epon_x*sin(pha)**2  ! pond. field, including time de

         ez(i) = 2*a0*ezhelm*cos(pha)*tpon  ! reconstruct EM fields (s-p
         by(i) = 2*a0*byhelm*sin(pha)*tpon
         az(i) = -2*a0*ezhelm*sin(pha)*tpon
         ey(i) = 0.
         bz(i) = 0.
      end do
      end









