

!  ======================================================
!
!                      AVLAS
!
!    Cycle-average diagnostics
!
!  ======================================================

subroutine avlas
  use bopsvars
  implicit none

  integer icycdel, i, ip, l
  real*8 xinc, dxlab, ese, xmelab, xmilab, xmplab, ute, u2, uti, utp, xclab,  xc1,xc2
  real*8 emte, emtm, rhom, dmax, bzl, bxl, byl, exl, eyl, ezl, ud

  real*8 uxl,uyl,uzl,gl
  real*8 rhlab(0:nx+1)

  !  absorbed laser energy = poynting flux at lh boundary
  !  averaged over 1 cycle
  !  - corrected for time taken for em wave to reflect.

  if (mod(itime,itav).eq.1) then

     ntav=ntav+1
     eta(ntav)=0.
     uinc(ntav)=0.
     utrans(ntav)=0.
     uback(ntav)=0.
     uesa(ntav)=0.
     utha(ntav)=0.
     uthea(ntav)=0.d0
     uthia(ntav)=0.d0
     uema(ntav)=0.
     utma(ntav)=0.
     urha(ntav)=0.
     ulha(ntav)=0.
     usysa(ntav)=0.
     utota(ntav)=0.
     xi0(ntav)=0.
     rhoim(ntav)=0.

  else

     !  absorption computed from
     !  poynting flux through lh boundary
     ud=ff(1)*ff(1)+fb(1)*fb(1)
     u1a(ntav)=u1a(ntav)+ud*dt/tav
     if (ud.eq.0) ud=a0*a0/2

     !  sim frame forward and backward em flux (invariant)
     !  TE mode

     if (ppol.ne.0) then
        uinc(ntav)=uinc(ntav)+ff(1)**2*dt/tav
        uback(ntav)=uback(ntav)+fb(1)**2*dt/tav
        utrans(ntav)=utrans(ntav)+ff(nx)**2*dt/tav
     endif

     !  TM mode
     if (spol.ne.0) then
        uinc(ntav)=uinc(ntav)+gf(1)**2*dt/tav
        uback(ntav)=uback(ntav)+gb(1)**2*dt/tav
        utrans(ntav)=utrans(ntav)+gf(nx)**2*dt/tav
     endif


     xi0(ntav)=xi0(ntav)+a0max**2*dt/tav

     !   potential energy (includes dc) - lab frame
     dxlab=dx*gam0
     ese=0.
     do i=1,nx+1
        ese = ese + 0.5*dxlab*(ex(i)+vy0*bz(i))**2
     end do
     uesa(ntav)=uesa(ntav)+ese*dt/tav


     xmelab=me/gam0**2
     xmilab=mi/gam0**2
     xmplab=mp/gam0**2

     !  electron thermal energy - lab frame
     ute=0.d0
     do ip=1,ne
        call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyl,gl)
        u2=ux(ip)**2+uz(ip)**2+uyl**2
        ute=ute+xmelab*u2/(gl+1.)
     end do

     !  ion thermal energies - lab frame
     uti=0.d0
     utp=0.d0

     do l=1,ni_tot
        ip=ne+l
        call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyl,gl)
        u2=ux(ip)**2+uz(ip)**2+uyl**2
 	if (species(ip).eq.2) then
          uti=uti+xmilab*u2/(gl+1.)  ! heavy ions
  	else if (species(ip).eq.3) then
          utp=utp+xmplab*u2/(gl+1.)  ! protons
	endif
     end do


     uthea(ntav)=uthea(ntav)+ute*dt/tav
     uthia(ntav)=uthia(ntav)+uti*dt/tav
     uthpa(ntav)=uthpa(ntav)+utp*dt/tav
! total thermal
     utha(ntav)=uthea(ntav)+uthia(ntav)+uthpa(ntav)

     !  escaped particles
     urha(ntav)=urha(ntav)+erhb*dt/tav

     !  escaped particles
     ulha(ntav)=ulha(ntav)+elhb*dt/tav

     !  total field energy (ES+EM) - lab frame
     emte=0.
     emtm=0.
     rhom=0.
     dmax=0.
     do i=1,nx
        bzl = bz(i)+vy0*ex(i)
        bxl = -vy0*ez(i)
        byl = by(i)/gam0
        eyl = ey(i)/gam0
        exl = ex(i)+vy0*bz(i)
        ezl = ez(i)
        emte = emte + 0.5*( eyl**2 + exl**2 + bzl**2)*dxlab
        emtm = emtm + 0.5*( ezl**2 + bxl**2 + byl**2)*dxlab
        !  lab frame ion density
        rhlab(i)=(rhoi(i)-vy0*jyi(i))/gam0
        dmax=max(rhlab(i),dmax)
     end do
     uema(ntav) = uema(ntav) + emte*dt/tav
     utma(ntav) = utma(ntav) + emtm*dt/tav

     !  total average
     usysa(ntav)=utha(ntav)+uema(ntav)+utma(ntav)
     utota(ntav)=usysa(ntav)+urha(ntav)+ulha(ntav)

     rhoim(ntav) = rhoim(ntav)+dmax/itav

     !  position of critical density (lab frame)
     !  xx(i) is x-axis array

     xclab=xcrit*gamma0
     dxlab=dx*gamma0
     xc1 = xclab/2.
     call findxc(xc1,xc2,rhlab,nx,dxlab,rhotrak)
     xcrni(ntav) = xc2/xconv
     xcrit = xc2/gamma0

! Estimated absorption fraction allowing for reflection path
     icycdel=(2*xm1/tav+0.5)
     xinc=uinc(ntav-icycdel)
     if(xinc.ne.0) then
        eta(ntav)=(uinc(ntav-icycdel)-uback(ntav)-utrans(ntav))/xinc
     else
        eta(ntav)=0.
     endif

  endif


end subroutine avlas
