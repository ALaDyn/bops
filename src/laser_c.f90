
!     ===================================
!
!     Launch counterpropagating laser field from RH boundary
!
!     ===================================

      subroutine laser_c
      use bopsvars
      implicit none

      real*8 t, pha_c
      real*8 tp_c, u, tpend_c, tptot_c, tp1tail_c, tp1end_c
      real*8 tdum
      t=dt*itime
      pha_c=w0_c*t
    

!  Correct for chirp (if any)
!  TODO: Need to correct pulse length tp if bandwidth to be held constant

      if (bchirp.ne.0) then
         pha_c = pha_c - bchirp*t**2
      endif

      if (ilas.eq.0) then
         a0max_c=0.

!   constant antenna at x=0 with linear rise-time

      else if (ilas.eq.1) then
         if (t.le.trise_c) then
            a0max_c=a0_c*t/trise_c
         else
            a0max_c=a0_c
         endif

!  Gaussian

      else if (ilas.eq.2) then
!  convert tfwhm to 1/e (tp ~ 0.6 tfwhm)
         tp_c=2.*sqrt(log(2.))*tpulse_c
! revised by Qiao for cut of counterpropagating laser
         
         u=min(20.d0,(t-tdel_c)*(t-tdel_c)/2./tp_c**2)
         if (t>2*tdel_c) then
            a0max_c=0
         else  
            a0max_c=a0_c*exp(-u)
         endif

!     triangular with flat top (intensity)

      else if (ilas.eq.4) then
         tpend_c=trise_c+tpulse_c
         if (t.le.trise) then
            a0max_c=a0_c*(t/trise_c)**0.5
         else if (t.gt.trise_c.and.t.le.tpend_c) then
            a0max_c=a0_c
         else if (t.gt.tpend_c.and.t.le.tpend_c+tfall_c) then
            a0max_c=a0_c*(1.-(t-tpend_c)/tfall_c)**0.5
         else
            a0max_c=0.
         endif

!     sin squared single and double pulse

      else if (ilas.eq.5 .or. ilas.eq.6) then
!     convert from fwhm to total length
         tptot_c=tpulse_c*2.13


         if (t.le.tptot_c) then
            a0max_c=a0_c*sin(pi*t/tptot_c+pi/2)
         else
            a0max_c=0.
         endif


!     double pulse

      else if (ilas.eq.9) then
         tp1tail_c=trise_c+tpulse_c
         tp1end_c=tp1tail_c+tfall_c
    
!     Pulse 1
         if (t.le.trise_c) then
            a0max_c=a0_c*(t/trise_c)**0.5
            fb(nx)=fb(nx)+a0max_c*sin(w0_c*t+pi/2)
         else if (t.gt.trise_c .and. t.le.tp1tail_c) then
            a0max_c=a0_c
            fb(nx)=fb(nx)+a0max_c*sin(w0_c*t+pi/2)
         else if (t.gt.tp1tail_c .and. t.le.tp1end_c) then
            a0max_c=a0_c*(1.-(t-tp1tail_c)/tfall_c)**0.5
            fb(nx)=fb(nx)+a0max_c*sin(w0_c*t+pi/2)

!     Pulse 2
         else
            a0max_c=0.
         endif


      endif

!     phase and polarization of LH attenna

      if (ilas.le.6) then
         fb(nx) = ppol_c*( a0max_c*sin(pha_c+pi/2) )
         gb(nx) = -spol_c*( a0max_c*cos(pha_c+pi/2))

!        write (*,'(f12.5,4(1pe12.5))') pha,ppol,spol,ff(1),gf(1)

      else if (ilas.eq.7) then

!   read attenna amplitudes from file
         read (85,*) tdum,fb(nx)
         read (86,*) tdum,gb(nx)
      endif



      end

