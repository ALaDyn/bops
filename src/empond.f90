
!     ==================================
!
!     Electromagnetic fields
!
!     - ponderomotive laser model for step-profile
!
!
!     ==================================

!     fpond(t,tpuls,sigma0,vosc,omega,x,y,z,epon_x,epon_y,epon_z,phipon)

      subroutine empond
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

!     linear rise

      t=dt*itime
      pha=w0*t
      if (ilas.eq.1) then
       Tpon = min(1.d0,t/trise)

      else if (ilas.eq.5) then
!     convert from fwhm to total length
         tptot=tpulse*2

         if (t.le.tptot) then
            Tpon=a0*sin(pi*t/tptot)
         else
            Tpon=0.
         endif
      endif

!     Standing wave in vacuum; evanescent in overdense plasma
!     Use standard solution of Helmholtz equation for step profile

      gamma_c = sqrt(1.+4*a0**2/nonc) ! gamma at surface
!      gamma_c = 1.
      wp_r = wp/sqrt(gamma_c)
      lskin = 1./wp_r   ! Rel. skin depth in EM units

!   Phase factor given by tan(phi) = -k0 * l_s = k0 * c/wp
      phi_m = atan(-w0/wp_r)

      do i=1,nx
         xc = (i-1)*dx - xcrit

         if (xc.le.0) then
!     vacuum solution
            chi_vac = xc + phi_m      ! Vacuum phase
!            xf = sin(2*vchi)
            ezhelm = sin(chi_vac)     ! laser 'Ez'
            byhelm = cos(chi_vac)     ! laser 'Bz'

         else
!   evanescent wave
!            xf = -2/lskin*sin(phi_m)**2*exp(-2*xc/lskin)
            ezhelm = sin(phi_m)*exp(-xc/lskin) ! laser Ez inside
            byhelm = cos(phi_m)*exp(-xc/lskin) ! laser Ez inside

         endif

!         gamma_pon= sqrt(1.+abs(phipon*Tpon)/2.) ! relativistic phi_pond
!         Epon_x = 2*a0**2*Tpon*xf/gamma_pon
!         epond(i) = Epon_x  ! pond. field

         ez(i) = 2*a0*ezhelm*cos(pha)*tpon  ! recons. EM fields (s-pol)
         by(i) = 2*a0*byhelm*sin(pha)*tpon
         az(i) = -2*a0*ezhelm*sin(pha)*tpon
      end do
      end









