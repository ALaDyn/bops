c  *************************************************
c
c    Electrostatic PIC code
c
c     slab geometry
c
c  *************************************************

      program pices
      include 'picesx.h'
      common/fftarr/ rw1(nffm),rw2(nffm)
c
      call init
      call derdat
      call arinit
      call denprof(1,ne)
      call v1q(1,ne,vte)
      call ixq(1,ne)
      call xion
      call calcu(np)
      call eden
      call addneut
c      call iden
      call fieldx
      call laser
      call diagno
      call pouti
      call pout
      call gsnap
c
      do 100 itime=1,nt
c  electrons
	if (mod(itime,isube).eq.0) then
	  call espush(1,ne,qome,dte)
	  call move(1,ne,dte)
	  call pbcs(1,ne,vte,me)
	  call eden
	endif
c  ions
	if (mod(itime,isubi).eq.0) then
	  call espush(ne+1,ni,qomi,dti)
	  call move(ne+1,ni,dti)
	  call pbcs(ne+1,ni,vti,mi)
	  call iden
	endif
	call fieldx
        call laser
	call pshif
	call diagno
	call pout
	call gsnap
 100  continue
       call ghist
       call savear
      end
c
c     ==========
c
      subroutine diagno
      include 'picesx.h'
      character ct*40,cbl*16


      if (inorm.eq.2) then
c  time in wp**-1
        tlb=dt*itime*wp
      else if (inorm.eq.3) then
c  time in fs
        tlb=dt*itime/w0PHz
      else
c  time in w0**-1
        tlb=dt*itime
      endif

c  convert cpu time to char variable
      cbl='                '
      call chr(tlb,1,ct,lct)
      lctime=12-lct
      ctime='  t='//ct(1:lct)//cbl(1:lctime)

      if (mod(itime,itc).ne.0) return
      call energy(1,np)
c      call slice3d
      ntc=ntc+1
      end
c
c     ==========
c
      subroutine init
      include 'picesx.h'
c
c   default data
      nt=1
      nx=10
      ne=10
      ni=0
      ipbc=1
      ifbc=1
      ilas=1
      xl=1.
      xrzone=10.
      dt=nx/xl
      z=1
      qome=-1.
      miome=1836.
      vy0=0.
      theta0=0.
      vte=0.
      vti=0.
      vxm=0.2
      vym=0.1
      uxma=1.0
      nvx=100
      nvy=200
      nf=20
      wp=0.1
      w0=1.0
      tp=1.0
      tp1=1.0
      tp2=1.0
      trise=1.0
      xm1=0.
      xm2=1.0
      xsol=xm2
      erhb=0.
      rsv2=0.
      idc=0
      itc=1
      ntc=0
      idia=1
      igr=1
      igxs=2
      iout=1
      itime=0
      inorm=1
      inprof=1
      ipskip=1
      nesc=0
      isgx=1
      isgy=1
      jcon=0
c      abo=0.5
      abo=0.
c      asm=0.5
      asm=0.0
      ampl=0.
      yi=(0.,1.)
      imaxw=1
c
c   input data
      read(10,lwf)
      end

c     ==========

      subroutine derdat
      include 'picesx.h'
      neold=ne
      ne0=ne
      ncell=ne/nx

      if (inorm.eq.2) then
c  lwf/beat-wave input deck - quantities in c/wp
	w0=1.0
	xl=xl/wp
        dt=dt/wp
	trise=trise/wp
	tp=tp/wp
	tdel=tdel/wp
        sigma=sigma/wp
	tp1=tp1/wp
	tp2=tp2/wp
	xrzone=xrzone/wp
c      else if (inorm.eq.3) then
c  exptl i/p:  density in cm**-3
c              Intensity in Wcm**-2
c              Wavelength in microns

      endif

      rwp=1./wp
      dx=xl/nx
      xdebye=vte/wp
      xdodx=xdebye/dx
      n0=wp*wp
      wl = 0.
      wr = xl
c  laser frequency in 10^15 s**-1 (PetaHz)
      w0PHz = 2*pi*c/(xlam*1.e-6)/1.e15
c
c   cloud charge normalised to give ncrit=1 (rhoc=-1)
      if (inprof.eq.1) then
	xload=xl
	qe=-n0*xload/ne
	ncell=nint(ne*dx/xload)
        xm1=0.
        xm2=xl
      else if (inprof.eq.2.or.inprof.eq.3) then
	xm2=xl
	if (inprof.eq.2) xsol=xm2
	qe=-n0/ne*(xm2-xsol/2-xm1/2)
      else if (inprof.eq.4) then
	qe=-n0/2/ne*(xm2+xsol2-xsol-xm1)
      endif
      qomi=-z*qome/miome
      qi=-z*qe
      me=qe*qome
      mi=me*miome
      np=ne+ni
      nxo2=nx/2
      nfo2=nf/2
      dto2=dt/2.
      rdx=1./dx
      dkx=2.*pi/xl
      dw=pi/nkm/2/dt
      dvx=2*vxm/nvx 
      tav=amin1(2*pi,igr*dt)
      itav=tav/dt
      ntav=0
      isube=1
      dte=dt*isube
      isubi=isubi*isube
      dti=dt*isubi
      if (ni.eq.0) isubi=nt+1
      isgx=nx/ncxm
      isgy=max(isgy,itav/ncym)
      if (mod(itav,ncym).ne.0) isgy=isgy+1
      end

c     ==========

      subroutine arinit
      include 'picesx.h'
      do  i=0,nx
	xx(i)=i*dx
	xk(i)=i*dkx
	ex(i)=0.
      end do

      do 170 i=1,nkm
 170    xw(i)=(i-1)*dw
c
c   smoothing and boost factor for fft
      do 200 i=1,nxo2
	dxko2=i*pi/nx
	smarg=amax1(abo*sin(dxko2)**2-asm*tan(dxko2)**4,-30.0)
	sm(i)=exp(smarg)
 200    rk2(i)=1./((2.0*sin(dxko2)/dx)**2)*sm(i)**2
      end

c  =======================================
c
c    Particle loading routines
c
c  =======================================

      subroutine ostart
      include 'picesx.h'
      dpx = xload/ne
      do i=1,ne
        xo(i) = dpx*(i-0.5) + xm1
      end do
      end


c     ==========

      subroutine xq(i1,n,xlo)
      include 'picesx.h'
      if (n.eq.0) return
c
c   binary reversal for positions
      dpx=xlo/n
      rs=0.
      do 100 l=1,n
	ip=l+i1-1
	xo(ip)=dpx/2.+rs*xlo+xm1
	rsi=1.0
  10      rsi=rsi/2.
	  rs=rs-rsi
	  if (rs.ge.0.) goto 10
	rs=rs+2*rsi
 100  continue
      end

c     ==========

      subroutine ixq(i1,n)
      include 'picesx.h'
      if (n.eq.0) return

c   binary reversal for indices
      rs=0.
      do 100 l=1,n
	ip=l+i1-1
	ipn=rs*n+i1
	xn(ipn)=xo(ip)
	rsi=1.0
  10      rsi=rsi/2.
	  rs=rs-rsi
	  if (rs.ge.0.) goto 10
	rs=rs+2*rsi
 100  continue

      do 200 l=1,n
	ip=l+i1-1
 200    xo(ip)=xn(ip)
      end

c     ==========

      subroutine xion
      include 'picesx.h'
      if (ni.eq.0) return
      do 100 l=1,ne
	ipi=l+ne
	xn(ipi)=xn(l)
	xo(ipi)=xn(ipi)
 100  continue
      end

c     ==========

      subroutine v1q(i1,n,vt)
      include 'picesx.h'
      real f0,df,rs,v,dv,vmax,finf
      if (n.eq.0) return
      nv=20*n
      vmax=3.0d0
      dv=vmax/nv
      f0=0.
      finf=sqrt(pi/2)
      l=1
      ip1=n/2
      ip2=n/2+1
      do 100 i=1,nv
	v=vmax-(i-0.5d0)*dv
	df=dexp(dmax1(-30.0d0,-0.5d0*v*v))*dv/finf*n/2.
	f0=f0+df
	if(f0.ge.l) then
c          vip=vt*(v-dv*(f0-l)/df)
	  vip=vt*v
	  if (imaxw.eq.1) then
	    vx(ip1)=vip
	    vx(ip2)=-vip
	  else if (imaxw.eq.2) then
	    vx(ip1)=0.
	    vx(ip2)=0.
	  endif
	  l=l+1
	  ip1=ip1-1
	  ip2=ip2+1
	endif
 100  continue
      end

c     ==========

      subroutine v2q(i1,n,vt)
      include 'picesx.h'
      double precision rs,rsi
      if (n.eq.0) return

c   trinary reversal for 2v velocity space
      rs=0.d0
      do 200 l=1,n
	ip=i1+l-1
	vm=vt*(-2.*alog((l-0.5)/n))**0.5
c        vm=amin1(vm,0.99)
	theta=2*pi*rs
	vx(ip)=vm*sin(theta)
	rsi=1.0d0
  20      rsi=rsi/3.
	  rs=rs-2*rsi
	  if (rs.ge.0.d0) goto 20
 200  rs=rs+3*rsi
      end


c  ==================================
c
c   Initialise density profile
c
c  ==================================

      subroutine denprof(i1,n)
      include 'picesx.h'
      real den,dxs,dmax
      if (n.eq.0) return

c   uniform
      if (inprof.eq.1) then
	call ostart
      endif

c   linear
      if (inprof.eq.2) then
	n1=n
	xsol2=xm2
	xsol=xm2

c   linear with flat top
      else if (inprof.eq.3) then
	n1=0.5*n0/(-qe)*(xsol-xm1)
	n2=n-n1
	xsol2=xm2

c  trapezoidal
      else if (inprof.eq.4) then
	n1=0.5*n0/(-qe)*(xsol-xm1)
	n3=0.5*n0/(-qe)*(xm2-xsol2)
	n2=n-n1-n3
      endif

c
c   linear with moat
      if (inprof.ge.2.and.inprof.le.4) then
	den=0.
	xplas=xsol-xm1
	dmax=2.*n1/xplas
	l=1
	nstrip=n1*5
	dxs=xplas/nstrip
	do 200 i=1,nstrip+1
	  x=xm1+i*dxs-dxs
	  den=den+dxs*(x-xm1)/xplas*dmax
	  if (den.ge.l) then
	    xo(l+i1-1)=x
	    l=l+1
	  endif
 200    continue
      endif

      end

c  =========================================
c
c    Prescribed laser field for
c     ponderomotive force
c
c  =========================================


      subroutine laser
      include 'picesx.h'
      if (ilas.eq.0) return
c
      t=dt*itime

      if (ilas.eq.1) then
c  Gaussian in x and t
        t0=tdel
        tau=tp
        do i=0,nx+1
          r = dx*i
          argr = amin1( r**2/2/sigma**2,20.)
          argt = amin1( (t-t0)**2/2/tau**2,20.)
          a(i) = a0*exp(-argr-argt)
        end do

      else if (ilas.eq.2) then
c  Gaussian x, constant with linear rise-time
        do i=0,nx+1
          r = dx*i
          argr = amin1( r**2/sigma**2,20.)
          if (t.lt.tp) then
            a(i) = a0*sqrt(t/tp)*exp(-argr)
          else 
            a(i) = a0*exp(-argr)
          endif
        end do

      endif

c  ponderomotive force (field) 
c   - excluding gamma factor because this contains particle momentum


      do i=1,nx
        epond(i) = 0.25*(a(i+1)**2 - a(i-1)**2)/2/dx
      end do
      epond(0) = 0.



      end

c  ==========================================
c
c    Find momenta from loaded velocities
c
c  ==========================================

      subroutine calcu(n)
      include 'picesx.h'
      do 100 l=1,n
	xn(l)=xo(l)
	g=(1-vx(l)**2)**(-0.5)
	ux(l)=g*vx(l)
 100  continue
      end

c  ===================================
c
c    Update particle momenta
c
c  ===================================  

      subroutine espush(ip1,n,qom,dts)
      include 'picesx.h'

      do l=1,n
	ip=ip1+l-1
 
c   interpolate field, fpond and a**2 (for gamma)
	xa = xo(ip)*rdx
	i1 = xa
	i2 = i1+1
	b2 = xa-i1
	b1 = 1.-b2

	exi = b1*ex(i1) + b2*ex(i2)
        epondi = b1*epond(i1) + b2*epond(i2)
        a2 = b1*a(i1)**2 + b2*a(i2)**2

c  net gamma including py=a
        gam1 = sqrt(1. + ux(ip)**2 + 0.5*a2 )
        
	ux(ip) = ux(ip) + qom*dts*( exi + epondi/gam1)

c  guessed value for new a2 = old a2 
c  - changes slowly compared to ux
	gam2 = sqrt( 1.+ ux(ip)**2 + 0.5*a2 )
	vx(ip) = ux(ip)/gam2
      end do

      end
 

c  ===================================
c
c    Update particle positions
c
c  ===================================  

 
      subroutine move(ip1,n,dts)
      include 'picesx.h'

      do l=1,n
	ip=ip1+l-1
	xn(ip)=xo(ip)+dts*vx(ip)
      end do
      end

c  ===================================
c
c    Register shift for positions
c
c  ===================================  


      subroutine pshif
      include 'picesx.h'
      do 100 l=1,np
 100    xo(l)=xn(l)
      end

c  ===================================
c
c    Particle boundary conditions
c
c  ===================================  

      subroutine pbcs(i1,n,vt,mass)
      include 'picesx.h'
      real mass

      do 100 l=1,n
	ip=i1+l-1

c   periodic
	if (ipbc.eq.1) then
	  if (xn(ip).lt.0.) xn(ip)=xn(ip)+xl
	  if (xn(ip).ge.xl) xn(ip)=xn(ip)-xl

c   reflective at x=wl and x=wr

	else if (ipbc.eq.2) then
	  if (xn(ip).le.wl) then
	    xn(ip) = 2*wl - xn(ip)
	    vx(ip) = -vx(ip)
	    ux(ip) = -ux(ip)
	  endif
	  if (xn(ip).ge.wr) then
	    xn(ip) = 2*wr - xn(ip)
	    vx(ip) = -vx(ip)
	    ux(ip) = -ux(ip)
	  endif 

	else if (ipbc.eq.3) then
c   re-inject with thermal velocity at rh boundary
	  if (xn(ip).le.wl) then
	    xn(ip) = 2*wl - xn(ip) 
	    vx(ip) = -vx(ip)
	    ux(ip) = -ux(ip)
	  endif
	  if (xn(ip).ge.wr) then
	    xn(ip) = 2*wr - xn(ip)
	    eth1=((ux(ip)*ux(ip)+1)**(0.5)-1)*mass
	    nesc=nesc+1
	    if (nesc.ge.npm) nesc=1
c  find new random position in cold maxwellian
	    iinj=n*rano(1.)+1
	    vm=vt*(-2.*alog((iinj-0.5)/n))**0.5
	    theta=2*pi*rano(1.)
	    vx(ip)=-vm*abs(cos(theta))
            g = sqrt(1.-vx(ip)**2)
	    eth2=(g-1)*mass
c  add net escape energy
	    erhb=erhb+eth1-eth2
            uesc(nesc) = eth1-eth2
	    ux(ip)=vx(ip)*g
	  endif
	endif
c  max and min momenta
	uxmax=max(uxmax,ux(ip))
	uxmin=max(uxmin,ux(ip))
 100  continue
      end


c  ======================================
c
c    Electron density
c
c  ======================================


      subroutine eden
      include 'picesx.h'
      re=qe/dx
      do 10 i=0,nx
  10  rhoe(i)=0.

      do 100 l=1,ne
	xa=xn(l)*rdx
	i1=xa
	i2=i1+1
	f2=xa-i1
	f1=1.-f2
	rhoe(i1)=rhoe(i1)+re*f1
	rhoe(i2)=rhoe(i2)+re*f2
 100  continue

      if (ipbc.eq.1) then
c   periodic 
	rhoe(0)=rhoe(0)+rhoe(nx)
	rhoe(nx)=rhoe(0)

      else if ((ipbc.eq.2).or.(ipbc.eq.3)) then
c   reflective
	rhoe(0) = 2*rhoe(0)
        rhoe(nx) = 2*rhoe(nx)
      endif

      rhodc=0.
      do 200 i=0,nx
 200    rhodc=rhodc+rhoe(i)+rhoi(i)
      end


c  ======================================
c
c    Fixed ions:  n = ne + no
c
c  ======================================


      subroutine addneut
      include 'picesx.h'
c
c   add neutralising background
      do 100 i=0,nx
 100    rhoi(i)=-rhoe(i)
      end


c  ======================================
c
c    Ion density
c
c  ======================================


      subroutine iden
      include 'picesx.h'
      if (ni.eq.0) return
      ri=qi/dx

      do 10 i=0,nx
  10  rhoi(i)=0.

      do 100 l=1,ni
	ip=ne+l
	xa=xn(ip)*rdx
	i1=xa
	i2=i1+1
	f2=xa-i1
	f1=1.-f2
	rhoi(i1)=rhoi(i1)+ri*f1
	rhoi(i2)=rhoi(i2)+ri*f2
 100  continue

      if (ipbc.eq.1) then
c   periodic 
	rhoi(0)=rhoi(0)+rhoi(nx)
	rhoi(nx)=rhoi(0)

      else if ((ipbc.eq.2).or.(ipbc.eq.3)) then
c   reflective
	rhoi(0) = 2*rhoi(0)
	rhoi(nx) = 2*rhoi(nx)
      endif

      do 200 i=0,nx
 200    rhot(i)=rhoe(i)+rhoi(i)
      end

c  ===============================
c
c    Poisson solver using FFTs
c
c  ===============================

      subroutine esf1
      include 'picesx.h'
      common/fftarr/ rw1(nffm),rw2(nffm)
      do 100 i=1,nx
	rw1(i)=rhoe(i)+rhoi(i)
 100    rw2(i)=0.0
      ifail=0
      call c06ecf(rw1,rw2,nx,ifail)
c   zero dc part
      rw1(1)=0.0
      rw2(1)=0.0
      yrhok(1)=(0.,0.)
      yphik(1)=(0.,0.)
      do 200 ipl=2,nxo2
	imi=nx+2-ipl
	yrhok(ipl)=rw1(ipl)+yi*rw2(ipl)
	yrhok(imi)=rw1(imi)+yi*rw2(imi)
c        rksq=1.0d0/((ipl-1)*(ipl-1)*dkx*dkx)
	rw1(ipl)=rw1(ipl)*rk2(ipl-1)
	rw2(ipl)=rw2(ipl)*rk2(ipl-1)
	rw1(imi)=rw1(imi)*rk2(ipl-1)
	rw2(imi)=rw2(imi)*rk2(ipl-1)
	yphik(ipl)=rw1(ipl)+yi*rw2(ipl)
 200    yphik(imi)=rw1(imi)+yi*rw2(imi)
      yrhok(nxo2+1)=rw1(nxo2+1)+yi*rw2(nxo2+1)
c      rksq=1.0d0/(nxo2*nxo2*dkx*dkx)
      rw1(nxo2+1)=rw1(nxo2+1)*rk2(nxo2)
      rw2(nxo2+1)=rw2(nxo2+1)*rk2(nxo2)
      yphik(nxo2+1)=rw1(nxo2+1)+yi*rw2(nxo2+1)
c
c   ift for potential
      call c06gcf(rw2,nx,ifail)
      call c06ecf(rw1,rw2,nx,ifail)
      call c06gcf(rw2,nx,ifail)
      do 300 i=1,nx
 300    phi(i)=rw1(i)
c
c   periodic bcs
      if (ifbc.eq.1) then
	phi(nx+1)=phi(1)
	phi(0)=phi(nx)
c   bounded
      else if (ifbc.eq.2) then
	phi(nx+1)=phi(1)
	phi(0)=phi(nx)
	b=0.5*(phi(nx)-phi(2))*rdx
	do 350 i=0,nx+1
 350      phi(i)=phi(i)+b*(i-1)*dx
      endif
c
c   es field ex
      do 400 i=1,nx
	ipl=i+1
	imi=i-1
 400    ex(i)=0.5*(phi(imi)-phi(ipl))*rdx
c   bcs
      if (ifbc.eq.1) then
	ex(0)=ex(nx)
	ex(nx+1)=ex(1)
      else if (ifbc.eq.2) then
	ex(nx+1)=0.
      endif
      end


c  ============================================
c
c    Direct field integration - slab geometry
c
c  ============================================

      subroutine fieldx
      include 'picesx.h'
      common/ftarr/ rw1(nffm),rw2(nffm)
      do i=0,nx
	rhot(i)=rhoe(i)+rhoi(i)
      end do

c   integrate div.e=rho directly (trapezium approx)
      if (ifbc.eq.2) then
c   end point - ex=0 mirror at wl
	ex(0)=0.

	do i=1,nx
          ex(i) = ex(i-1) + 0.5*dx*(rhot(i)+rhot(i-1))
        end do

c   potential - same again
	phi(0)=0.
	do i=1,nx
          phi(i) = phi(i-1) - 0.5*(ex(i)+ex(i-1))*dx
        end do 
      endif
      end

c  =====================================
c
c   Field FT  E(k)
c
c  =====================================

      subroutine esft
      include 'picesx.h'
      common /fftarr/rw1(nffm),rw2(nffm)
      real uesk(0:nxm),wk1(nkm),wk2(nkm)
      complex cx(nkm)
      integer iwk(20)
c
c real work array for CRAY, double prec. for IBM
c
      do 50 i=1,nx-2
  50    wk1(i)=ex(i)
      do 55 i=nx-1,nkm
  55    wk1(i)=0.
      ifail=0
c     call c06faf(rw1,nx,rw2,ifail)
c  IMSL fft
      call fftrc(wk1,nkm,cx,iwk,wk2)
      uesk(0)=0.
      do 60 i=2,nx/2
	j=i-1
   60   uesk(j)=abs(cx(i))
      uekm=0.
      uesk(0)=0.
      xmax=w0/wp*4
      xmin=0.
      dax=xmax/5.
      nk=(nx-2)/2
      igk=1+nk/2000
      do 100 i=1,nk
 100    uekm=amax1(uekm,uesk(i))
      ukmax=5*10**(ze(uekm)+1)
      ukmin=ukmax/1.e6
      do 300 j=1,nk
	xk(j)=j*dkx
 300    uesk(j)=amax1(uesk(j),ukmin)
c  plot ft Ek**2 on (manual) log-lin
      call grxy(xk,uesk,nk,730+idc,igk,-2
     :,'      k        ','     Ues(k)    ','uesk'//ctime(1:12)
     :,1./wp,1.,xmin,xmax,dax,ukmin,ukmax,1.)
      end
c
c     ==========
c
      subroutine vdist
      include 'picesx.h'
      do 50 i=1,nvx
	work2(i)=0.
  50    work1(i)=0.
      do 100 l=1,ne
	vax=abs(vx(l)+vxm)/dvx
	i0=vax+1
	j0=vay+1
	if (i0.le.nvx) work1(i0)=work1(i0)+1
 100  continue
c
      do 150 i=1,nvx
 150    xv(i)=i*dvx-dvx/2.-vxm
      call grxy(xv,work1,nvx,600+idc,igxs,1
     :,'      vx       ','     f(vx)     ','fvxe'//ctime(1:12),1.,1.)
      end
c
c     ==========
c
      subroutine udist
      include 'picesx.h'
      do 50 i=1,nvx
  50    work1(i)=0.
      umax=(1.-vxm*vxm)**(-0.5)-1.
      du=umax/nvx
c
c  relativistic thermal energy distn
      do 100 l=1,ne
	u=(sqrt(ux(l)*ux(l)+1.)-1)/du
	i0=u+1
	if (i0.le.nvx) work1(i0)=work1(i0)+1
 100  continue
      fmin=work1(1)
      do 200 i=1,nvx
	if (work1(i).ne.0) fmin=amin1(fmin,work1(i))
 200  continue
      if (fmin.le.0) fmin=1
      do 210 i=1,nvx
	if (work1(i).ne.0) then
	  work2(i)=alog(work1(i))
	else
	  work2(i)=log(fmin)
	endif
	xv(i)=i*du-du/2
 210  continue
      call grxy(xv,work2,nvx,670+idc,igxs,1,2
     :,'      U(keV)   ','     fe(U)     ','fuep'//ctime(1:12),1.,1.)
      end

c     ==========
c
      subroutine energy(i1,n)
      include 'picesx.h'
      real mass

      complex yes
      if (i1.eq.1) then
	mass=me
	uth(ntc)=0.
	udr(ntc)=0.
      else
	mass=mi
      endif

c   potential energy

      ese=0.
      do 100 i=0,nx
	ese = ese + 0.5*ex(i)**2*dx
 100  continue
      ues(ntc)=ese

c   drift energy

      vs=0.
      do 150 l=1,n
	ip=i1+l-1
150     vs=vs+vx(ip)
      vs=vs/n
      udr(ntc)=udr(ntc)+0.5*vs*vs*n*mass

c   thermal energy - relativistic

      do 200 l=1,n
	ip=i1+l-1
	uth(ntc) = uth(ntc) + ((ux(ip)**2+1.)**(0.5)-1.)
     :*mass
 200  continue

c   energy lost to rh boundary
      erh(ntc)=erhb

c   total
      utot(ntc)=ues(ntc)+uth(ntc)

c   ratio field to thermal
      r1=0.
      if (uth(ntc).ne.0) r1=ues(ntc)/uth(ntc)
      esoth(ntc)=r1

c  Max electron density
      demax=0.
      do i=0,nx
        demax=amax1(demax,abs(rhoe(i)))
      end do
      denem(ntc) = demax
c  Density on axis - ave. of 1st 3 cells
      dene0(ntc) = abs(rhoe(0)+rhoe(1)+rhoe(2))/3.

      end

c     ==========

      subroutine savear
      include 'picesx.h'
      write (8,101) np,ne,ni
  101 format(3i8)
      write (8,102) (xo(i),xn(i),ux(i),vx(i),i=1,np)
  102 format(2(1pe16.8)/4(1pe16.8))
      write (8,103) nx,(ex(i),i=1,nx)
  103 format(i8/(1pe16.8))
      write (8,104) itime
  104 format(i8)
      end
c
c     ==========
c
      subroutine pouti
      include 'picesx.h'
      character*40 cw

c     write(11,picohd)
      call i0prnt('itime ',itime)
      call r0prnt('  tmax',nt*dt*wp,'f12.2')
      call i0prnt('    nt',nt)
      call i0prnt('    ne',ne)
      call i0prnt('    ni',ni)
      call i0prnt('    nx',nx)
      call i0prnt(' ncell',ncell)
      call i0prnt('  ipbc',ipbc)
      call i0prnt('  ifbc',ifbc)
      call i0prnt('inprof',inprof)
      call i0prnt('   igr',igr)
      call i0prnt('  idia',idia)
      call i0prnt('  iout',iout)
      call r0prnt('    dt',dt*wp,'f14.4')
      call r0prnt('    xl',xl*wp,'f12.2')
      call r0prnt('    dx',dx*wp,'f14.4')
      call r0prnt('    a0',a0,'f12.2')
      call r0prnt(' sigma',sigma*wp,'f12.2')
      call r0prnt('tpulse',tp*wp,'f12.2')
      call r0prnt('  tdel',tdel*wp,'f12.2')
      call r0prnt('w0/PHz',w0PHz,'f12.2')
      call r0prnt('lambda',xlam,'f12.2')
      call r0prnt('     z',z,'f11.1')
      call r0prnt('    wp',wp,'f12.2')
      call r0prnt('   vte',vte,'f13.3')
      call r0prnt('   vti',vti,'f13.3')
      call r0prnt(' miome',miome,'f11.1')
      call r0prnt('  qome',qome,'f13.3')
      call r0prnt('  qomi',qomi,'f13.3')
      call r0prnt('    qe',qe,'f13.3')
      call i0prnt('  isgx',isgx)
      call i0prnt('  isgy',isgy)

      if (nx.ge.nxm.or.np.gt.npm) then
	cw=' not enough array space'
	call warn(cw)
	stop
      endif
      if (xdodx.lt.1.) then
	cw=' dx less than debye length '
	call warn(cw)
      endif
      end
c
c     ==========
c
      subroutine pout
      include 'picesx.h'
      if (mod(itime,iout).ne.0) return
      call blank
      write (6,'(a)') ctime
      call i0prnt(' itime',itime)
      call r0prnt('  time',itime*dt*wp,'f12.2')
c     call r1prnt('rhot  ',rhot,nx)
c     call r1prnt('ex    ',ex,nx)
      end

c  =========================================
c
c   Graphical snapshots
c
c  =========================================

      subroutine gsnap
      include 'picesx.h'
      if (mod(itime,igr).ne.0) return
      xlno=xl*wp/w0
      exno=1.0

c  max and min momenta
      uxmax=0.
      uxmin=0.
      do 150 l=1,ne
	uxmax=max(uxmax,ux(l))
150     uxmin=min(uxmin,ux(l))

c  electron px-x phase space
      call grps(xn,ux,ne,0.,xlno,uxmin,uxmax,100+idc,ipskip
     :,'       x       ','      px       ','pxxe'//ctime(1:12) 
     :,wp,1.)

c  ion density
      call grxy(xx,rhoi,nx+1,220+idc,igxs,-1
     :,'       x       ','        ni     ','nimv'//ctime(1:12) 
     :,wp,1.,0.,xlno,0.2*xlno,0,2*n0,0.5*n0)
c      call vdist
c      call udist
c      call uescd
      if (itime.eq.0) then
	idc=idc+1
	return
      endif

c  electron density
      demax=10.
      demin=0.
      call grxy(xx,rhoe,nx+1,200+idc,igxs,1
     :,'       x       ','      ne       ','eden'//ctime(1:12) 
     :,wp,-rwp**2,0.,xlno,0.2*xlno,demin,demax,demax/2.)

c  pulse amplitude
      call grxy(xx,a,nx+1,500+idc,igxs,-1
     :,'       x       ','      vosc/c   ','alas'//ctime(1:12) 
     :,wp,1.,0.,xlno,0.2*xlno,0,a0,0.5*a0)

c  fpond
      call grxy(xx,epond,nx+1,600+idc,igxs,1
     :,'       x       ','      epond    ','epon'//ctime(1:12) 
     :,wp,1.)

c  e.s. potential
      call grxy(xx,phi,nx+1,250+idc,igxs,1
     :,'       x       ','        phi    ','phip'//ctime(1:12) 
     :,wp,1.)

c  e.s. field
      exmx=0.16
      exmn=0.
      call grxy(xx,ex,nx+1,300+idc,igxs,1
     :,'       x       ','        Ex     ','expw'//ctime(1:12) 
     :,wp,rwp,0.,xlno,0.2*xlno,exmn,exmx,exmx/2.)

      call esft
      idc=idc+1
      end
c
c     ==========
c
      subroutine slice3d
      include 'picesx.h'
      write (70,'(1pe13.4/)') (ex(i),i=1,nx+1)
      end

c  =========================================
c
c   Time histories
c
c  =========================================

      subroutine ghist
      include 'picesx.h'
      if (ntc.eq.0) return
      ntc=ntc-1
      do 100 i=1,ntc
 100    xt(i)=dt*(i-1)*itc

c  Electrostatic energy
      call grxy(xt,ues,ntc,1510,igxs,1
     :,'       t       ','      Ues      ','uese            '
     :,wp,1.)

c  Thermal energy
       call grxy(xt,uth,ntc,1520,igxs,1
     :,'       t       ','      Uth      ','ueth            '
     :,wp,1.)

c  Total plasma energy
      call grxy(xt,utot,ntc,1530,igxs,1
     :,'       t       ','      Utot     ','utot            '
     :,wp,1.)

C  Max electron density
      call grxy(xt,denem,ntc,1600,igxs,1
     :,'       t       ','      ne(max)   ','nemx            '
     :,wp,rwp**2)

C  electron density at x=0
      call grxy(xt,dene0,ntc,1610,igxs,1
     :,'       t       ','      ne(x=0)   ','nem0            '
     :,wp,rwp**2)
      end

c  =========================================
c
c   Ouput graphical data in ODPLOT format
c
c  =========================================

      subroutine grxy(xx,yy,nd,id,igx,iman,chx,chy,ctitle
     :,rnx,rny,xmin,xmax,dax,ymin,ymax,day)

      common/pplo/iplot,nppage
      real x(2002),y(2002),store(15),xx(0:nd),yy(0:nd)
      real xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,dax,day,rnx,rny
      integer iman
      character*8 chars(1),chx*15,chy*15,ctitle*16,cfile*12,cid*1

c  open data file
      call chr(mod(id,10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xydata'
      write (15,'(a)') cfile
      open (51,file=cfile)

      if (igx.eq.0) igx=1
      n=nd/igx
      if (mod(nd,igx).ne.0) n=n+1
      do i=1,nd,igx
	ip=(i+igx-1)/igx
	x(ip)=xx(i-1)*rnx
        y(ip)=yy(i-1)*rny
      end do

      idp=id
      idash=1
      ispage=1
      if (iman.lt.0) then
c  manual axes: ityp=-iman; id=-id
	ityp=-iman
	idp=-id
      else if (iman.eq.0) then
c  plot using old set of axes
	ityp=1
	ispage=2
	idash=2
      else if (iman.gt.0) then
c  auto axes: ityp=iman
	ityp=iman
      endif
c
c  write out to 'ODPLOT' datafile

      write (50,101) idp,n,ityp,idash,0,ispage,chx,chy,ctitle
c  manual axes
      if (iman.lt.0) then
	write(50,102) xmin,xmax,dax,ymin,ymax,day
      endif
      write (51,103) (x(i),y(i),i=1,n)
      close(51)

  101 format (2i6,4i4,2a15,1a16)
  102 format (6(1pe12.3))
  103 format (2(1pe11.3))

      end

c  =========================================
c
c   Ouput phase-space data in ODPLOT format
c
c  =========================================

      subroutine grps(x,y,n,xmin,xmax,ymin,ymax,id,iskip,chx,chy
     :,ctitle,rnx,rny)
      common/pplo/iplot,nppage
      real x(0:n),y(0:n),store(15),xp(10000),yp(10000)
      real xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,rnx,rny
      character*8 chars(1),chx*15,chy*15,ctitle*16,cfile*12,cid*1

c  open data file
      call chr(mod(id,10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xydata'
      open (61,file=cfile)

      if (ymin.eq.ymax) ymax=ymax+1.
      iskip=iskip+n/6000/iskip
      np=n/iskip
      if (mod(n,iskip).ne.0) np=np+1
      do 100 i=1,n,iskip
	iadd=iskip*rano(1.)
	ip=(i+iskip-1)/iskip
	xp(ip)=x(i+iadd)*rnx
 100    yp(ip)=y(i+iadd)*rny

      idp=-id
      day=(ymax-ymin)/5
      dax=(xmax-xmin)/5

      write (50,101) idp,np,1,0,9,1,chx,chy,ctitle
c  manual axes
      write(50,102) xmin,xmax,dax,ymin,ymax,day
      write (61,103) (xp(i),yp(i),i=1,np)
      close(61)

  101 format (2i6,4i4,2a15,1a16)
  102 format (6(1pe12.3))
  103 format (2(1pe15.5))


      end
c
c  ===================
c  header page
c
      subroutine header
      include 'picesx.h'
      npars=18
      write (50,101) 0,npars
  101 format(2i6)
      call hop('    nt',nt*1.,0)
      call hop('  tmax',nt*dt*wp,1)
      call hop('delay ',tdel*wp,1)
      call hop('rezone',xrzone*wp,1)
      call hop('#elecs',ne*1.,0)
      call hop('# ions',ni*1.,0)
      call hop('vosc/c',a0,4)
      call hop('  I 16',a0*a0*140,4)
      call hop('om/omp',1./wp,2)
      call hop('n0    ',n0,1)
      call hop('    Te',vte*vte*511,3)
      call hop('    Ti',vti*vti*511,3)
      call hop(' mi/me',miome,1)
      call hop('     Z',z,0)
      call hop('    xl',xl*wp,1)
      call hop('    nx',nx*1.,0)
      call hop('    dx',dx*wp,4)
      call hop('ldebye',vte,4)
      end
c
c  ====================
c
      subroutine hop(cname,z,ndp)
      character cname*6
      write (50,101) cname,z,ndp
  101 format(a6,1pe15.6,i6)
      end
c
c  ====================
c
      subroutine minmax(f,fmin,fmax,n)
      real f(n)
      real fmin,fmax
      xmin=f(1)
      xmax=f(n)
      do 100 i=1,n
      xmin=amin1(xmin,f(i))
      xmax=amax1(xmax,f(i))
 100  continue
      fmin=xmin
      fmax=xmax
      end
      subroutine dum
      entry c06gcf
      entry c06ecf
      end
