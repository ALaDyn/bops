!     ======================================================
!
!     Derive initial simulation conditions from input data
!
!     ======================================================

      subroutine derdat
      use bopsvars
      implicit none

      real*8 ampp, xlol
      integer isp1, i, im1

      write (*,*) "uimax",uimax
!     boost velocity and gamma
      the0=theta0*pi/180
      vy0=sin(the0)
      gam0=1.d0/cos(the0)
      gamma0=gam0
!      uy0 = vy0*gam0
      uy0 = tan(the0)

!     units
      if (iunits.eq.2) then
!     time in fs, length in microns
         tconv = 1.88/xlambda
         xconv = 2*pi/xlambda

!     ES system: time normalised to 1/wp, lengths to c/wp
!     - only good for normal incidence

      else if (iunits.eq.1) then
         wp=sqrt(nonc)
         tconv = 1./wp
         xconv = 1./wp

      else
!     default: time in w_0^-1, length in c/w_0
         wp=sqrt(nonc)
         w0=1.
         tconv = 1.
         xconv = 1.

      endif

!   output  conversion factor from w0^-1 -> fs
         tcfs = 1.88/xlambda
!  from c/w0 to mu
         xmu = 2*pi/xlambda

!     time and length unit conversion
!     box lengths xl,xm1,xsol input in lab frame -> convert to sim
!     frame with boost factor

      xllab=xl*xconv
      xl=xl/gam0*xconv
      xm1=xm1/gam0*xconv
      xsol=xsol/gam0*xconv
      xsol2=xsol2/gam0*xconv
      xm2=xm2/gam0*xconv
      dfoil = dfoil/gam0*xconv
      x_layer=x_layer/gam0*xconv
      gap_foil = gap_foil/gam0*xconv
      xcur1 = xcur1/gam0*xconv
      xcur2 = xcur2/gam0*xconv
      dx=xl/nx
      dt=dx
      xload=xload*xl
      rho0lab=nonc
!  Sim frame time slowed by gam0**2

      if (ioboost==0) then
	ttrans=gam0**2
      else
	ttrans=1.
      endif

      trun = trun*tconv/ttrans
      tpulse = tpulse*tconv/ttrans
      trise = trise*tconv/ttrans
      tfall = tfall*tconv/ttrans

!  chirp in /ps**2 - convert to /w0**2
      w0r=1.88/xlambda*1.e3
      bchirp = bchirp/w0r**2

!     timesteps
      if (trun.gt.0) then
         nt = trun/dt
         if (tgr.gt.0) then
	   igr = tgr*tconv/dt/ttrans
	   write (*,*) 'tgr=',tgr
           write (*,*) 'igr=',igr
	 else 
           igr = igr*tconv/dt/ttrans
	   iplot2d = iplot2d*tconv/dt/ttrans
           write (*,*) 'igr=',igr
           write (*,*) 'iplot2d=',iplot2d
         endif
         iout = iout*tconv/dt/ttrans
         itc = 2*pi/dt/itc  ! Sample itc x per cycle
         nftem = nftem*tconv/dt/ttrans-1

      else
!     run time given as # timesteps
         trun = nt*dt

      endif


!     a0 i/p as Intensity if +ve, vosc/c if -ve
      if (a0.lt.0) then
         a0 = -a0
      else
         a0 = sqrt(a0*xlambda**2/1.38e18)
         a02 = sqrt(a02*xlambda**2/1.38e18)
      endif

!     evaluate polarization vectors

!   convert to upper case
      if (cpolzn.eq.'s') cpolzn='S'
      if (cpolzn.eq.'c') cpolzn='C'
      if (cpolzn.eq.'p') cpolzn='P'

!     s-pol
      if (cpolzn.eq.'S') then
         ppol=0.
         spol=1.
	 push_scheme=2

!     circular
      else if (cpolzn.eq.'C') then
         if (ppol.eq.0 .or. spol.eq.0) then
!  ensure pol. vectors defined if amplitudes not set
            ppol=1.0
            spol=1.0
         endif
         ampp=ppol**2+spol**2
         ppol=ppol/sqrt(ampp)
         spol=spol/sqrt(ampp)
	 push_scheme=3

!     p-pol by default
      else
         ppol=1.
         spol=0.
      endif

      if (nonc.gt.0) then
!     lab frame electron density from i/p
         nonc=nonc*gam0**3
         wp=nonc**0.5
         rho0=nonc
         rho_layer=rho_layer*gam0**3
      else
         rho0=wp*wp
      endif

!     adjust ramp scale-length
      if (inprof.le.3) then
!         xsol=xm1+xlolam*2*pi*rho0lab/gam0

      else if (inprof.eq.7) then
         xsol = xm1+dfoil
      endif

!     density scale length
      gnon=(xsol-xm1)/rho0
!     ginzburg parameter (k0l)**2/3*sin2(theta)
      xgin=gnon**(2./3.)*vy0*vy0
      if (xlolam.gt. 0.1) then
         theta_opt = 180/pi*asin(0.8/(xlolam*2*pi)**0.33)
      else
         theta_opt = 0.
      endif
!     mass ratio
      miome=miome*amass

!     Temperatures in keV
! TODO: check the conditions!!
!#### Anupam & Bin 2009/2010
      if (Te.gt.0) then
      vte=sqrt(Te/511)
      endif
      if (Ti.gt.0) then
      vti=sqrt(Ti/511./miome)
      endif
      if (Tproton.gt.0) then
      vtp=sqrt(Tproton/511./mpome)
      endif
!#### Anupam & Bin 2009/2010
!
      duin = duin/511.
!     Left and right particle walls
      wl=5*dx/2.
      wr=xl-dx/2.
      xdebye=vte/wp
      xdodx=xdebye/dx


      nf=max(min(nt,nftem)/ift,1)

      nxo2=nx/2
      nfo2=nf/2
      dto2=dt/2.
      dkx=2.*pi/xl
      dw=2*pi/nf/(ift*dt)
      tav=min(2*pi,igr*dt)
      itav=tav/dt
      ntav=0
      icstart=tcstart/dt
!     isube=max(itav/100,1)
      isube=1
      dte=dt*isube
      isubi=isubi*isube
      dti=dt*isubi
      if (ni.eq.0) isubi=nt+1
      a0max=a0

!     limit output switches
!      iout=max(iout,nt/50)
!      igr=max(igr,nt/5)
      itsk=max(itsk,nt/1000)

!     time-averages
      ncyc=min(ncyc,igr/itav)
      ncyc=max(ncyc,1)


!     movie over 4 laser cycles starting at isp(1)
      if (nsp.lt.0) then
         nsp=itav*2
!     start movie when laser reaches edge of plasma
         isp1 = xm1/dt
         do i=1,nsp
            isp(i) = isp1 + 8*(i-1)
         end do
      endif
      end












