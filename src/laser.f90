
!     ===================================
!
!     Launch laser field from LH boundary
!
!     ===================================

      subroutine laser
      use bopsvars
      implicit none

      real*8 t, pha, pha2, a0max2
      real*8 tp, u, tpend, tptot, tptot2, tp1tail, tp1end
      real*8 tp2tail, tp2end, t2, tdum, tp2beg
      real*8 :: hdot, hdash  ! pulse envelope time and space derivatives
      real*8 :: eys, ezs, bys, bzs  ! Field antennae

      t=dt*itime
      pha=w0*t*ttrans  ! keep same # launch cycles as for normal incidence
      pha2 = w02*t*ttrans

!  Correct for chirp (if any)
!  TODO: Need to correct pulse length tp if bandwidth to be held constant

      if (bchirp.ne.0) then
         pha = pha - bchirp*t**2
      endif

      a0max2=0. ! 2nd pulse off by default
     hdot=0.
     hdash=0.  ! Long pulse by default

     pulse_shape: select case(ilas)
      case(0)
         a0max=0.


      case(1,11)
!   constant antenna at x=0 with linear rise-time
         if (t.le.trise) then
            a0max=a0*t/trise
         else
            a0max=a0
         endif
	 

      case(2)
!  Gaussian
!  convert tfwhm to 1/e (tp ~ 0.6 tfwhm)
         tp=1./(2.*sqrt(log(2.)))*tpulse
         u=min(20.d0,(t-tdel)*(t-tdel)/2./tp**2)
         a0max=a0*exp(-u)


      case(4)
!     triangular with flat top (intensity)
         tpend=trise+tpulse
         if (t.le.trise) then
            a0max=a0*(t/trise)**0.5
         else if (t.gt.trise.and.t.le.tpend) then
            a0max=a0
         else if (t.gt.tpend.and.t.le.tpend+tfall) then
            a0max=a0*(1.-(t-tpend)/tfall)**0.5
         else
            a0max=0.
         endif


       case(5,6,15)
!     sin squared single and double pulse
!     convert from fwhm to total length
         tptot=tpulse*2
         tptot2=tpulse2*2

         if (t.le.tptot) then
            a0max=a0*sin(pi*t/tptot)
	    hdot= a0*pi/tptot*cos(pi*t/tptot)
	    hdash = -hdot
         else
            a0max=0.
	    hdot = 0.
 	    hdash = 0.
         endif
!     2nd pulse
         if (t.le.tptot2 .and. tptot2.gt.0) then
            a0max2=a02*sin(pi*t/tptot2)
         else
            a0max2=0.
         endif


      case(9)
!     double pulse
         tp1tail=trise+tpulse
         tp1end=tp1tail+tfall
         tp2beg=tp1end+tsep
         tp2tail=trise2+tpulse2
         tp2end=tp2tail+tfall2
         t2=t-tp2beg

!     Pulse 1
         if (t.le.trise) then
            a0max=a0*(t/trise)**0.5
            ff(1)=ff(1)+a0max*sin(w0*t)
         else if (t.gt.trise .and. t.le.tp1tail) then
            a0max=a0
            ff(1)=ff(1)+a0max*sin(w0*t)
         else if (t.gt.tp1tail .and. t.le.tp1end) then
            a0max=a0*(1.-(t-tp1tail)/tfall)**0.5
            ff(1)=ff(1)+a0max*sin(w0*t)

!     Pulse 2
         else if (t2.gt.0 .and. t2.le.trise2) then
            a0max=a02*(t2/trise2)**0.5
            ff(1)=ff(1)+a0max*sin(w02*t2)
         else if (t2.gt.trise2 .and. t2.le.tp2tail) then
            a0max=a02
            ff(1)=ff(1)+a0max*sin(w02*t2)
         else if (t2.gt.tp2tail .and. t2.le.tp2end) then
            a0max=a02*(1.-(t2-tp2tail)/tfall2)**0.5
            ff(1)=ff(1)+a0max*sin(w02*t2)
         else
            a0max=0.
         endif


      end select pulse_shape


!     phase and polarization of LH attenna

      if (ilas.le.6) then
	 eys = ppol*( a0max*sin(pha) - hdot*cos(pha) )
	 bzs = ppol*( a0max*sin(pha) - hdot*cos(pha) )
	 ezs = spol*(-a0max*cos(pha) - hdot*sin(pha) )
	 bys = spol*( a0max*cos(pha) + hdot*sin(pha) )
         ff(1) = 0.5*(eys + bzs)  ! F+ at LHB
         gf(1) = 0.5*(ezs - bys)  ! G- at LHB

!        write (*,'(f12.5,4(1pe12.5))') pha,ppol,spol,ff(1),gf(1)


      else if (ilas.eq.7) then

!   read attenna amplitudes from file
         read (85,*) tdum,ff(1)
         read (86,*) tdum,gf(1)
      endif

      Uemin=Uemin+(ff(1)**2+gf(1)**2)*dt  ! Incoming Poynting flux

      end

