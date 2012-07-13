
!     ================================
!
!     Graphical snapshots: frequency
!     determined by igr
!
!     ================================


         subroutine gsnap
         use bopsvars
         implicit none

         real*8 xln, uzma, uzmi, uymi, uxmi, xmin,nmin
         real*8 vzma, vzmi, vymi, vyma, gammin, gammax
         real*8 vxma, vxmi !#### Anupam & Bin 2009/2010: vxmi & vxma added!
         real*8 nmax, xstep, nstep
         real*8, dimension(0:nx+1) :: work1, work2, work3
         ! plot work arrays for back transformation
         real*8, dimension(0:nx+1) :: rhoi_p, rhoe_p, rhot_p, rhop_p
         real*8, dimension(0:nx+1) :: jyi_p, jzi_p, jye_p, jze_p, jyt_p, jzp_p
         real*8, dimension(0:nx+1) :: ex_p, ey_p, ez_p, bx_p, by_p, bz_p, ay_p
         integer i

         character chx*15
         integer ii

!#### Anupam & Bin 2009/2010: for outputing reflected Ez
!      real*8, dimension(0:nx+1) :: ezreflect
!#### Anupam & Bin 2009/2010

         if (iunits.eq.2) then
            chx = '   x/microns   '
            xln = xllab/xmu
         else if (iunits.eq.1) then
            chx = '     k_px      '
            xln = xllab/xconv
         else
            chx = '    k_0x       '
            xln = xllab
         endif
         xmin=0.
         xstep=xln/5.
         nmin=0.
         nmax=2*rho0lab
         nstep = rho0lab/2.

!     find max/min of momenta in lab/boost frame
         uxma=ux(1)
         uxmi=ux(1)
         uzma=uz(1)
         uzmi=uz(1)

         if (ioboost.eq.1) then
            uyma=uy(1)
         else
            uyma=gam0*( uy(1) - vy0*gamma(1) )
         endif

         uymi=uyma

         do i=1,ne
            uux(i) = ux(i)
            uuz(i) = uz(i)

            if (ioboost.eq.1) then
!     boost frame
               uuy(i) = uy(i)
               xnlab(i) = xn(i)/xconv
               vy(i) = uuy(i)/gamma(i)   ! velocities
               vz(i) = uuz(i)/gamma(i)
               vx(i) = uux(i)/gamma(i)   ! anupam 2010:added vx
               gamlab(i) = gamma(i)

            else
               !     lab frame
               uuy(i) = gam0*( uy(i) - vy0*gamma(i) )
               xnlab(i) = gam0*xn(i)/xconv
               vy(i) = uuy(i)/gamma(i)   ! velocities
               vz(i) = uuz(i)/gamma(i)
               vx(i) = uux(i)/gamma(i)   ! anupam 2010:added vx
               gamlab(i) = gamma(i)

            endif

            uxma=max(uxma,uux(i))
            uxmi=min(uxmi,uux(i))
            uyma=max(uyma,uuy(i))
            uymi=min(uymi,uuy(i))
            uzma=max(uzma,uuz(i))
            uzmi=min(uzmi,uuz(i))
         end do

!#### Anupam & Bin 2009/2010
!TODO: change uyma, uzma, uxma
!         uyma = 3*a0
         uymi = -uyma
!         uzma = 3*a0
         uzmi = -uzma
!         uxma = 2*a0
         uxmi = -uxma
         vyma = 1.5
         vymi = -vyma
         vzma = 1.5
         vzmi = -vyma
         gammin=0.5
         gammax = sqrt(1+2*uyma**2)

!#### Anupam & Bin 2009/2010

	call grps(xnlab(1:ne),uux(1:ne),ne,xmin,xln &
	,uxmi,uxma,1000+idc,ipskip &
	,chx,'      px       ','pxxe'//ctime(1:12) )

!#### Anupam & Bin 2009 : add output of "vxxe"

!	call grps(xnlab(1:ne),vx(1:ne),ne,xmin,xln &
!        ,vxmi,vxma,11600+idc,ipskip &            
!        ,chx,'      vx       ','vxxe'//ctime(1:12) )
!##### Anupam & Bin 2009/2010

!     if (ppol.ne.0) then
	call grps(xnlab(1:ne),uuy(1:ne),ne,xmin,xln &
	,uymi,uyma,1500+idc,ipskip &
	,chx,'      py       ','pyxe'//ctime(1:12) )

	call grps(uux(1:ne),uuy(1:ne),ne,uxmi,uxma &
	,uymi,uyma,1600+idc,ipskip &
	,'      py       ','      px       ','pxpy'//ctime(1:12) )
	call grps(xnlab(1:ne),vy(1:ne),ne,xmin,xln &
	,vymi,vyma,11500+idc,ipskip &
	,chx,'      vy       ','vyxe'//ctime(1:12) )

	call grps(xnlab(1:ne),gamlab(1:ne),ne,xmin,xln &
	,gammin,gammax,12500+idc,ipskip &
	,chx,'      gamma    ','gaxe'//ctime(1:12) )
!     x, ux, uy plot in rayshade format
!     call gr3d(xnlab,uux,uuy,0,ne,ipskip,idc,'elec_x2v'//ctime(1:8))

!     endif

!#### Anupam & Bin 2009/2010 : output for either TM or TE mode

!if (spol.ne.0) then
	call grps(xnlab(1:ne),uuz(1:ne),ne,xmin,xln &
	,uxmi,uxma,1700+idc,ipskip &
	,chx,'      pz       ','pzxe'//ctime(1:12) )

	call grps(uux(1:ne),uuz(1:ne),ne,uxmi,uxma &
	,uxmi,uxma,1800+idc,ipskip &
	,'      pz       ','      px       ','pxpz'//ctime(1:12) )

	call grps(xnlab(1:ne),vz(1:ne),ne,xmin,xln &
	,vzmi,vzma,11700+idc,ipskip &
	,chx,'      vz       ','vzxe'//ctime(1:12) )

!     x, ux, uy plot in rayshade format
!     call gr3d(xnlab,uux,uuz,0,ne,ipskip,idc,'elec_x2v'//ctime(1:8))

!endif		!#### Anupam & Bin 2009/2010


!  ions
! TODO: need to separate heavy ions and protons

        if (ni_tot.gt.0) then
	uxma=ux(ne+1)
	uxmi=uxma
!#### Anupam & Bin 2009/2010: added uzma and uzmi
         uzma=uz(ne+1)
         uzmi=uzma
!#### Anupam & Bin 2009/2010

	uyma=gam0*( uy(ne+1) - vy0*gamma(ne+1) )
	uymi=uyma

	do i=1,ni_tot
!     ion index
	ii = i+ion1-1
	uux(i)=ux(ii)

	if (ioboost.eq.1) then
!     boost frame
	uuy(i) = uy(ii)
	xnlab(i) = xn(ii)

!#### Anupam & Bin 2009/2010: output for vx,vy,vz
        vy(i)=uuy(i)/gamma(ii)
        vz(i)=uuz(i)/gamma(ii)
        vx(i)=uux(i)/gamma(ii)
        gamlab(i)=gamma(ii)
!#### Anupam & Bin 2009/2010

	else
!     lab frame
	uuy(i)=gam0*(uy(ii) - vy0*gamma(ii) )
	xnlab(i) = gam0*xn(ii)/xconv
	endif

	uxma=max(uxma,uux(i))
	uxmi=min(uxmi,uux(i))
	uyma=max(uyma,uuy(i))
	uymi=min(uymi,uuy(i))
!#### Anupam & Bin 2009/2010
        uzma=max(uzma,uuz(i))
        uzmi=min(uzmi,uuz(i))
!#### Anupam & Bin 2009/2010
	end do

!#### Anupam & Bin 2009/2010: setup for velocity max. and min.
         vxmi=-1.5
         vxma=1.5
         vymi=-1.5
         vyma=1.5
         vzmi=-1.5
         vzma=1.5

         gammin=0.5
         gammax=sqrt(1+uxma**2+uyma**2+uzma**2)
!#### Anupam & Bin 2009/2010

	uxma=2*sqrt(uimax/.511/mpome)
!#### Anupam & Bin 2009/2010: seperate outputing information of ions and protons
        if (np.gt.0) call grps(xnlab(1:np),uux(1:np),np,xmin,xln &
              		,-uxma,uxma,9200+idc,ipskip & 
              		,chx,'      px       ','pxxp'//ctime(1:12) ) 
!#### Anupam & Bin 2009/2010
	if ((ni_tot-np).gt.0) call grps(xnlab(1:ni_tot),uux(1:ni_tot),ni_tot,xmin,xln &
             ,-uxma,uxma,1200+idc,ipskip &
             ,chx,'      px       ','pxxi'//ctime(1:12) )

!#### Anupam & Bin 2009/2010: add others output for ions 

         if (np.gt.0) call grps(xnlab(1:np),uuy(1:np),np,xmin,xln &
              ,uymi,uyma,9300+idc,ipskip &
              ,chx,'      py       ','pyxp'//ctime(1:12) )
         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),uuy(np+1:ni_tot),ni_tot-np,xmin,xln &
              ,uymi,uyma,1300+idc,ipskip &            
              ,chx,'      py       ','pyxi'//ctime(1:12) )
        
         if (np.gt.0) call grps(xnlab(1:np),uuz(1:np),np,xmin,xln &
              ,uzmi,uzma,9400+idc,ipskip &
              ,chx,'      pz       ','pzxp'//ctime(1:12) )
         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),uuz(np+1:ni_tot),ni_tot-np,xmin,xln &
              ,uzmi,uzma,1400+idc,ipskip &            
              ,chx,'      pz       ','pzxi'//ctime(1:12) )

!         call grps(uux(1:ni_tot),uuy(1:ni_tot),ni_tot,uxmi,uxma &
!              ,uymi,uyma,3400+idc,ipskip &            
!              ,'      py       ','      px       ','pypi'//ctime(1:12) )
!         call grps(uux(1:ni_tot),uuz(1:ni_tot),ni_tot,uxmi,uxma &
!              ,uzmi,uzma,3300+idc,ipskip &            
!              ,'      pz       ','      px       ','pzpi'//ctime(1:12) )

         if (np.gt.0) call grps(xnlab(1:np),vx(1:np),np,xmin,xln &
           ,vxmi,vxma,9500+idc,ipskip &
           ,chx,'      vx       ','vxxp'//ctime(1:12) )
!         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),vx(np+1:ni_tot),ni_tot-np,xmin,xln &
!           ,vxmi,vxma,3200+idc,ipskip &            
!           ,chx,'      vx       ','vxxi'//ctime(1:12) )
!
         if (np.gt.0) call grps(xnlab(1:np),vy(1:np),np,xmin,xln &
           ,vymi,vyma,9600+idc,ipskip &
           ,chx,'      vy       ','vyxp'//ctime(1:12) )
         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),vy(np+1:ni_tot),ni_tot-np,xmin,xln &
           ,vymi,vyma,4200+idc,ipskip &
           ,chx,'      vy       ','vyxi'//ctime(1:12) )

         if (np.gt.0) call grps(xnlab(1:np),vz(1:np),np,xmin,xln &
           ,vzmi,vzma,9700+idc,ipskip &
           ,chx,'      vz       ','vzxp'//ctime(1:12) )
         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),vz(np+1:ni_tot),ni_tot-np,xmin,xln &
           ,vzmi,vzma,5200+idc,ipskip &
           ,chx,'      vz       ','vzxi'//ctime(1:12) )

         if(np.gt.0) call grps(xnlab(1:np),gamlab(1:np),np,xmin,xln &
           ,gammin,gammax,9800+idc,ipskip &
           ,chx,'      gamma    ','gaxp'//ctime(1:12) )
         if ((ni_tot-np).gt.0) call grps(xnlab(np+1:ni_tot),gamlab(np+1:ni_tot),ni_tot-np,xmin,xln &
           ,gammin,gammax,3100+idc,ipskip &        
           ,chx,'      gamma    ','gaxi'//ctime(1:12) )

!#### Anupam & Bin 2009/2010

!     x, ux, uy plot in rayshade format
!     call gr3d(xnlab,uux,uuy,ion1,ni,ipskip
!     :,idc,'ions_x2v'//ctime(1:8),gam0*xsol,2*vti)

      endif


!!  FIELD PLOTS

!  Do back-transformations if lab-frame I/O requested

        do i=1,nx+1
	if (ioboost.eq.1) then
!  Boost frame fields
		jye_p(i) = jye(i)
		jyi_p(i) = jyi(i)
		jyt_p(i) = jp(i)    ! TODO: clean up current notation

		ay_p(i) = ay(i)
		rhoe_p(i) = rhoe(i)
		rhoi_p(i) = rhoi(i)
		rhop_p(i) = rhop(i)
		rhot_p(i) = rhot(i)
                bz_p(i) = bz(i)
		bx_p(i) = bx(i)
		by_p(i) = by(i)
		ey_p(i) = ey(i)
		ex_p(i) = ex(i)
		ez_p(i) = ez(i)
           else
!  Lab frame
		jye_p(i) = (jye(i) - vy0*rhoe(i))/gam0
		jyi_p(i) = (jyi(i) - vy0*rhoi(i))/gam0
		jyt_p(i) = (jp(i) - vy0*rhot(i))/gam0
		jze_p(i) = jze(i)/gam0**2
		jzi_p(i) = jzi(i)/gam0**2
!		jzp_p(i) = jzp(i)/gam0**2  ! TODO: need proton current

		ay_p(i)=gam0*(ay(i) - vy0*phi(i))
		rhoe_p(i) = (rhoe(i) - vy0*jye(i))/gam0
		rhoi_p(i) = (rhoi(i) - vy0*jyi(i))/gam0
		rhot_p(i) =  (rhot(i) - vy0*jp(i))/gam0
!		rhop_p = (rhop(i) - vy0*jyp(i))/gam0 
                bz_p(i) = bz(i)+vy0*ex(i)
		bx_p(i) = -vy0*ez(i)
		by_p(i) = by(i)/gam0
		ey_p(i) = ey(i)/gam0
		ex_p(i) = ex(i)+vy0*bz(i)
		ez_p(i) = ez(i)

	endif
       end do

if (itime.eq.0) then
!   initial ion density - fixed y-axis

	call grxy(xx,rhoi_p,nx,2100+idc,igxs,-1 &
	,chx,'       n_i/n_c ','ninc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,1.5*nonc,nstep)
!   initial proton density
	call grxy(xx,rhop_p,nx,2200+idc,igxs,-1 &
	,chx,'       n_p/n_c ','npnc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,1.5*nonc,nstep)
!   initial electron density
	call grxy(xx,-rhoe_p,nx,2300+idc,igxs,-1 &
	,chx,'       n_e/n_c ','nenc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,1.5*nonc,nstep)
!  Net charge density rhoe+rhoi+rhop
	call grxy(xx,rhot_p,nx+1,2000+idc,igxs,-1 &
	,chx,' \sum_j n_j    ','rhot'//ctime(1:12) &
	,xmin,xln,xstep,-nonc/1000.,nonc/1000.,nonc/2000.)
	idc=idc+1

	return
endif

!     ion density
call grxy(xx,rhoi_p,nx,2100+idc,igxs,-1 &
	,chx,'         ni/nc ','ninc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!  proton density
call grxy(xx,rhop_p,nx,2200+idc,igxs,-1 &
	,chx,'       n_p/n_c ','npnc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!  electron density
call grxy(xx,-rhoe_p,nx,2300+idc,igxs,-1 &
	,chx,'       n_e/n_c ','nenc'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

call grxy(xx,rhot_p,nx+1,2000+idc,igxs,1 &
	,chx,'       rho     ','rhot'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)


call grxy(xx,phi,nx+1,2500+idc,igxs,1 &
	,chx,'       phi     ','phsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

call grxy(xx,ex_p,nx+1,3000+idc,igxs,1 &
	,chx,'        E_x    ','exsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     if (ppol.ne.0) then

!     TE fields and sources

call grxy(xx,ey_p,nx+1,4000+idc,igxs,1 &
	,chx,'        E_y    ','eysi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

call grxy(xx,bz_p,nx+1,5000+idc,igxs,1 &
	,chx,'        B_z    ','bzsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

call grxy(xx,ff,nx+1,40000+idc,igxs,1 &
	,chx,'        f+     ','tefo'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

call grxy(xx,fb,nx+1,40500+idc,igxs,1 &
	,chx,'        f-     ','teba'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)


!     electron current

call grxy(xx,jye_p,nx+1,41000+idc,igxs,1 &
	,chx,'        jye    ','jyel'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     ion current

call grxy(xx,jyi_p,nx+1,43000+idc,igxs,1 &
	,chx,'        jyi    ','jyio'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     net current

call grxy(xx,jyt_p,nx+1,41500+idc,igxs,1 &
	,chx,'        jp     ','jtot'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     vector potential

call grxy(xx,ay_p,nx+1,42500+idc,igxs,1 &
	,chx,'        ay     ','ayem'//ctime(1:12) &
	,xmin,xllab,0.25*xllab,-2*a0,2*a0,a0)

!     endif
!  pond force  Epond from pond. laser model
do i=1,nx+1
	work1(i)=epond(i)/sqrt(1+az(i)**2)  ! correct for fluid gamma
end do
	call grxy(xx,work1,nx+1,28100+idc,igxs,1 &
,chx,'      fp       ','fpon'//ctime(1:12) & 
,xmin,xllab,0.25*xllab,-2*a0,2*a0,a0)

!#### Anupam & Bin 2009/2010: snapshoting for either TM or TE mode
!if (spol.ne.0) then ! changed

!     TM fields and sources

	call grxy(xx,ez_p,nx+1,3500+idc,igxs,1 &
	,chx,'        Ez     ','ezsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

	call grxy(xx,bx_p,nx+1,4500+idc,igxs,1 &
	,chx,'        Bx     ','bxsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

	call grxy(xx,by_p,nx+1,5500+idc,igxs,1 &
	,chx,'        By     ','bysi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

	call grxy(xx,gf,nx+1,44000+idc,igxs,1 &
	,chx,'        g+     ','tmfo'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

	call grxy(xx,az,nx+1,42000+idc,igxs,1 &
	,chx,'        az     ','azsi'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

	call grxy(xx,gb,nx+1,44500+idc,igxs,1 &
	,chx,'        g-     ','tmba'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     electron current

	call grxy(xx,jze_p,nx+1,45000+idc,igxs,1 &
	,chx,'        jze    ','jzel'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)

!     ion current

	call grxy(xx,jzi_p,nx+1,45500+idc,igxs,1 &
	,chx,'        jzi    ','jzio'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)


!     net current
	
	call grxy(xx,jzp_p,nx+1,46000+idc,igxs,1 &
	,chx,'        jzp    ','jzto'//ctime(1:12) &
	,xmin,xln,xstep,nmin,nmax,nstep)


!     vector potential


!     call grxy(xx,work1,nx+1,4200+idc,igxs,-1
!     :,chx,'        ay     ','ayem'//ctime(1:12)
!     :,xmin,xllab,0.25*xllab,-2*a0,2*a0,a0)

!endif !#### Anupam & Bin 2009/2010



end











