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
      trun = trun*tconv
      tpulse = tpulse*tconv
      trise = trise*tconv
      tfall = tfall*tconv

!  chirp in /ps**2 - convert to /w0**2
      w0r=1.88/xlambda*1.e3
      bchirp = bchirp/w0r**2

!     timesteps
      if (trun.gt.0) then
         nt = trun/dt
         if (tgr.gt.0) then
	   igr = tgr*tconv/dt
	   write (*,*) 'tgr=',tgr
           write (*,*) 'igr=',igr
	 else 
           igr = igr*tconv/dt
	   iplot2d = iplot2d*tconv/dt
           write (*,*) 'igr=',igr
           write (*,*) 'iplot2d=',iplot2d
         endif
         iout = iout*tconv/dt
         itc = 2*pi/dt/itc  ! Sample itc x per cycle
         nftem = nftem*tconv/dt-1

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

! ---------------------------------------------------------
!  TODO: Move this section to 'denprof' (target_config)

!     electron cloud charge normalised to give ncrit=1 (rhoc=-1)
!     all space quantities now code frame

      if (inprof.eq.1) then
!     uniform
         if (ipbc.eq.2.or.ipbc.eq.3.or.ipbc.eq.4) then
            xm2=xl-dx*0.5
            xm1=dx*0.5
         else if (ipbc.eq.1) then
            xm1=0.
            xm2=xl
         endif
         wl=xm1
         wr=xm2
         xload=xm2-xm1
         xcrit =xl/2.
         qe=-rho0*xload/ne

      else if (inprof.eq.2.or.inprof.eq.3) then
!     linear/linear with flat top
         xm2=xl-1.5*dx
         wr = xm2
         xm1=max(xm1,wl)
         if (inprof.eq.2) xsol=xm2
         qe=-rho0/ne*(xm2-xsol/2-xm1/2)
         xcrit = xm1 + 2*pi*xlolam/gam0

      else if (inprof.eq.4) then
!     trapezoidal
         wr = xl-1.5*dx
         wrf = wr
         qe=-rho0/2/ne*(xm2+xsol2-xsol-xm1)
         xcrit = xm1 + 2*pi*xlolam/gam0

      else if (inprof.eq.5) then
!     exponential
         wr = xl-1.5*dx
         xm2 = wr
         xsol2 = wr
         xload = (xsol - xm1)
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xlol*(1.-exp(-xload/xlol)) + xm2 - xsol)
         xcrit = xsol-xlol*log(rho0lab)

      else if (inprof.eq.6) then
!     tanh
         xm2=xsol2
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xm2-xm1/2.-xsol/2.)
         xcrit = xsol-xlol*log(rho0lab)

      else if (inprof.eq.7) then
!     foil target in middle of grid
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

!#### Anupam & Bin 2009: for double/three layer targets
      else if (inprof.eq.9) then
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

      target: select case(target_config)

      case(3)        ! three layers with proton layer on both front and rear sides
         n_pp=nint(-rho_layer*2*x_layer/qe) ! add particle numbers for 
                                            ! both electrons and prtons in the proton layers   
         n_pp1 = nint(n_pp/2.0)		   ! the particle number (electrons and ions) for 1st proton layer 
           n_pp2=n_pp-n_pp1                ! the particle number (electrons and ions) for 2nd proton layer
         
         case(4)     ! double layers with proton layer on the rear side   
           n_pp=nint(-rho_layer*x_layer/qe)  ! add particle numbers for
                                            ! both electrons and prtons in the proton layers
           n_pp1 = 0.0d0                    ! the particle number (electrons and ions) for 1st proton layer
           n_pp2=n_pp-n_pp1                 ! the particle number (electrons and ions) for 2nd proton layer
      
         case(5)     ! double layers with proton layer on the front side
           n_pp=nint(-rho_layer*x_layer/qe) ! add particle numbers for
                                            ! both electrons and prtons in the proton layers
           n_pp1 = n_pp                     ! the particle number (electrons and ions) for 1st proton layer
           n_pp2=n_pp-n_pp1                 ! the particle number (electrons and ions) for 2nd proton layer

        case default
           ! do nothing
        end select target

!#### Anupam & Bin 2009/2010 : inprof == 10 : multispecies target 

      else if (inprof.eq.10) then
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

         n_pp=nint(-rho_layer*xload/qe)       !add particle numbers for
                                              !both electrons and prtons for proton species
         n_pp1 = n_pp
         n_pp2=n_pp-n_pp1

!#### Anupam & Bin 2009/2010

      else if (inprof.eq.17) then
!    picket-fence 'foam' target in middle of grid
!   nfoil foils of thickness dfoil, spacing gap_foil
         xm2=xm1+dfoil*nfoil + gap_foil*(nfoil-1)
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=dfoil
         qe=-rho0*xload/ne*nfoil
         xcrit = xm1

      else if (inprof.eq.8) then
!     layered
         xsol2=xl-dx/2.
         qe=-(rho_layer*(xm2-xm1) + rho0*(xsol2-xsol))/ne
         xcrit = xm1

      else if (inprof.eq.57) then
!     foil target in middle of grid with exponential front side
         wr = xl-1.5*dx
         wl=dx+0.5
         xm2=xsol+dfoil
         xsol2=xm2
         xload = (xsol - xm1)
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xlol*(1.-exp(-xload/xlol)) + xm2 - xsol)
         xcrit = xsol-xlol*log(rho0lab)
         inprof = 5

      else if (inprof.eq.575) then
!  'soft' foil target in middle of grid with exponential front and back side
         wr = xl-1.5*dx
         wl=dx+0.5
         xsol2=xsol+dfoil
         xload = (xsol - xm1)   ! ramp length
         xm2 = xsol2 + xload
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xlol*(2.-2.*exp(-xload/xlol)) + xsol2 - xsol)
         xcrit = xsol-xlol*log(rho0lab)
         inprof = 5
      endif

!     position of first ion in particle arrays
      ion1=ne+1

!     NB: qe(sim) = qe(lab)*gam0**2.  Same for me,mi etc

! Define particle macro charges and masses
! qome=-1     already set in init

!#### Anupam & Bin 2009/2010 : 
      qi = -z*qe	! ion macro-charge
!     to calculate the density rho0_i=ni*qi=z*ni*qe=z*ne*qe
!     this needs to be divided by z

      qp = -qe   	! proton macro-charge
      me = qe*qome	! electron macro-mass (qome=-1)
      mi = me*miome	! ion macro-mass
      mp = me*mpome	! proton macro-mass
      qomi = qi/mi	! Zm_e/Am_p
      qomp = qp/mp	! q_p/m_p

!#### Anupam & Bin 2009/2010 :Total # simulation particles
      ne=n_pp+ne        ! in total electrons add particles for proton layers
      np=n_pp           ! as in bopsvars.f90, np is the particle number for proton layers 
      npart = ne+ni+np  ! for ni, the charge is qi=-z*qe; for np, the charge is qp=-qe
      ni_tot=ni+np      ! total # ions

!#### Anupam & Bin 2009/2010: position of first ion in particle arrays
!     change calculated ion1 after adding the particles for proton layer 
      ion1=ne+1          ! It is very IMPORTANT  !!!!!!!!!

!#### Anupam & Bin 2009/2010

!     initial drift energies - thermal 'zero'
      uthe0=ne*(gam0-1.d0)*me  ! electrons
      uthi0=ni*(gam0-1.d0)*mi  ! heavy ions
      uthp0=np*(gam0-1.d0)*mp  ! protons

! -------------------------------------------

      nf=min(nt,nftem)/ift

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












