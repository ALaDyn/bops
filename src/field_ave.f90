
!     ================================================================
!
!     Field averaging:
!
!     lab/boost frame selected with ioboost
!
!     ================================================================

subroutine fieldave
  use bopsvars
  implicit none

  integer icycav, i, l, i1, i2
  real*8 tcycav, sth, cth, dxl, exl1, exl2, re, xan,f1,f2,p0,rhel, phil
  real*8 pys, gams, pyl, gaml, poynt, gam1
  real*8 exl, eyl, ezl, bxl, byl, bzl, c2, ft1, rdx
  real*8  jxl,jyl,jzl,vxl,vyl,vzl,uxi,uyi,uzi,gami

  real*8 curx(0:nx+1),cury(0:nx+1),vye(0:nx+1),vye2(0:nx+1),vze(0:nx+1)
  real*8 wk1(0:nx+1),wk2(0:nx+1)

  rdx=1./dx
  icycav=ncyc*itav
  tcycav=ncyc*tav
  sth=vy0
  cth=(1.-vy0*vy0)**(0.5)
  if (mod(itime,igr).eq.0) then
     if (itime.gt.0) call grav
     !     rezero
     do i=1,nx+1
        avni(i)=0.
        avex(i)=0.
        avey(i)=0.
        avez(i)=0.
        avbx(i)=0.
        avby(i)=0.
        avbz(i)=0.
        avjy(i)=0.
        avjz(i)=0.
        avfp(i) = 0.
        exlab(i)=0.
        bzlab(i)=0.
        avext(i)=0.
        avphi(i)=0.
        dcphi(i)=0.
        dcbx(i)=0.
        dcby(i)=0.
        dcbz(i)=0.
        dcex(i)=0.
!#### Anupam & Bin 2009/2010 : for diagnose of all three directions
        dcey(i)=0.
!#### Anupam & Bin 2009/2010
        dcez(i)=0.
        edotj(i)=0.
        vxb(i)=0.
        vxb2(i)=0.
        dcvy(i)=0.
        dcjy(i)=0.
        dcjz(i)=0.
     end do

  else if (mod(itime,igr).ge.igr-icycav) then

     !     Evaluate E.J in lab frame: transform ux, uy, ex, qe, dx
     !     qe = qe'/gam0**2; dx = gam0*dx'
     do i=0,nx+1
        curx(i)=0.0
        cury(i)=0.
        wk1(i)=0.0
     end do

     !     electrons
     re=qe*rdx/gam0**3

     do  l=1,ne
        xan=xn(l)*rdx
        i1=xan+1
        i2=i1+1
        f1=i1-xan
        f2=1.-f1

        !     boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)
        vxl=uxi/gami
        vyl=uyi/gami

        jxl=re*vxl
        jyl=re*vyl
        cury(i1)=cury(i1)+jyl*f1
        cury(i2)=cury(i2)+jyl*f2
        curx(i1)=curx(i1)+jxl*f1
        curx(i2)=curx(i2)+jxl*f2
     end do

     p0=vy0*gam0

     !     electron 'fluid' vy
     !     use vy = jy/ne

     do i=1,nx+1
        rhel = (rhoe(i) - vy0*jye(i))/gam0
        jyl = cury(i)

        !     use py'=ay'+p0
        pys = ay(i)+p0
        gams = sqrt(1.+pys**2)
        pyl = gam0*(pys-vy0*gams)
        gaml = sqrt(1.+pyl**2+az(i)**2)  ! lab frame gamma
        vye(i) = pyl/gaml
        vze(i) = az(i)/gaml

     end do


     !     add ion current
     re=qi*rdx/gam0**3
     do  l=ne+1,ne+ni
        xan=xn(l)*rdx
        i1=xan+1
        i2=i1+1
        f1=i1-xan
        f2=1.-f1

        !     boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)
        vxl=uxi/gami
        vyl=uyi/gami

        jxl=re*vxl
        jyl=re*vyl
        cury(i1)=cury(i1)+jyl*f1
        cury(i2)=cury(i2)+jyl*f2
        curx(i1)=curx(i1)+jxl*f1
        curx(i2)=curx(i2)+jxl*f2
     end do

     call filter1(cury,nx+1)
     call filter1(curx,nx+1)

     !     cumulative e.j integral
     wk2(0)=0.0

     do i=1,nx+1
        wk1(i)=curx(i)*(ex(i)+vy0*Bz(i))+cury(i)*ey(i)/gam0
        wk2(i)=wk2(i-1)+wk1(i)*dxl
     end do

     !     time-av normalised to lab-frame Poynting vector = 0.5*vosc**2 cos(theta)
     poynt = 0.5*a0**2/gam0
     do i=1,nx+1
        edotj(i)=edotj(i)+wk2(i)/poynt*dt/tcycav
     end do

     c2=(cos(the0))**2
     ft1 = dt/tcycav
     do i=1,nx+1

        if (ioboost.eq.1) then
           !     boost frame fields
           bxl = 0.
           byl = by(i)
           bzl = bz(i)
           exl = ex(i)
           eyl = ey(i)
           ezl = ez(i)
           jyl = jp(i)
           rhel = rhot(i)
        else
           !     lab frame fields bz and ex to compare with 2D
           bzl = bz(i)+vy0*ex(i)
           bxl = -vy0*ez(i)
           byl = by(i)/gam0
           eyl = ey(i)/gam0
           exl = ex(i)+vy0*bz(i)
           ezl = ez(i)
           jyl = (jp(i) - vy0*rhot(i))/gam0
           jzl = jzp(i)/gam0**2
	   phil = gam0*(phi(i) - vy0*ay(i))
        endif

        gam1 = sqrt(1.+ay(i)**2 + az(i)**2)
        avex(i) = avex(i)+ft1*exl**2
        avey(i) = avey(i)+ft1*eyl**2
        avez(i) = avez(i)+ft1*ezl**2
        avbx(i) = avbx(i)+ft1*bxl**2
        avby(i) = avby(i)+ft1*byl**2
        avbz(i) = avbz(i)+ft1*bzl**2
        dcbx(i) = dcbx(i)+ft1*bxl
        dcby(i) = dcby(i)+ft1*byl
        dcbz(i) = dcbz(i)+ft1*bzl
        dcex(i) = dcex(i)+ft1*exl
!##### Anupam & Bin 2009/2010: added dcey dcez for 3D diagnose
        dcey(i) = dcey(i)+ft1*eyl     
        dcez(i) = dcez(i)+ft1*ezl     
!##### Anupam & Bin 2009/2010
        avphi(i) = avphi(i)+ft1*phil**2
        dcphi(i) = dcphi(i)+ft1*phil
        vxb(i) = vxb(i) + ft1*vye(i)*bzl
        vxb2(i) = vxb2(i) - ft1*vze(i)*byl
        avfp(i) = avfp(i) + ft1*epond(i)/gam1 ! need to correct f
        dcvy(i) = dcvy(i) + ft1*vye(i)
        dcjx(i) = dcjx(i) + ft1*curx(i)
        dcjy(i) = dcjy(i) + ft1*jyl
        dcjz(i) = dcjz(i) + ft1*jzl
        avjy(i) = avjy(i) + ft1*cury(i)**2
        avjz(i) = avjz(i) + ft1*jzl**2
     end do

  endif
end subroutine fieldave















