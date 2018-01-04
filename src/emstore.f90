!
!      EMSTORE
!
!     $Revision: 31 $
!
!      Store EM fields for spectra
!
!     ========================================================

      subroutine emstore
      use bopsvars
      implicit none
      integer nemw, iprobe, ipw, ipcurr
      real*8 fflab, fblab, gflab, gblab
      save nemw
      data nemw/1/

      if (ilas.eq.0) return
!  Time-integrated Poynting flux at LHB -SIM frame
      fflab = ff(1)
      fblab = fb(1)
      gflab = gf(1)
      gblab = gb(1)

      if (mod(itime,ift).eq.0) then

!  TE mode
        if (ppol.ne.0) then
!          Uemin = Uemin + dt*ift*fflab**2
!          Uemout = Uemout + dt*ift*fblab**2
!          Uemtr = Uemtr + dt*ift*ff(nx)**2
          Upoy_in = Upoy_in + dt*ift*Ey(1)*Bz(1)
          Upoy_out = Upoy_out + dt*ift*Ey(nx)*Bz(nx)
        endif

!  TM mode
        if (spol.ne.0) then
!          Uemin = Uemin + dt*ift*gflab**2
!          Uemout = Uemout + dt*ift*gblab**2
!          Uemtr = Uemtr + dt*ift*gf(nx)**2
          Upoy_in = Upoy_in - dt*ift*Ez(1)*By(1)
          Upoy_out = Upoy_out - dt*ift*Ez(nx)*By(nx)
        endif

!  probe antenna for forward wave
        if (inprof.eq.6) then
          iprobe=(xm2+xsol)/dx/2.
        else
        iprobe=nx
        endif


!  forward em
          ffl(nemw) = ff(iprobe)
          gfl(nemw) = gf(iprobe)
!  backward em
          fb0(nemw) = fb(1)
          gb0(nemw) = gb(1)


          ipw=min(xl/dx-dx,xcrit/dx)
          fex(nemw)=ex(ipw)

!  lab frame electron current at solid
          ipcurr=xsol/dx
          fjy(nemw) = (jye(ipcurr) - vy0*rhoe(ipcurr))/gam0



          if (mod(nemw,nf).eq.0 .and. mod(em_scheme,10).eq.1) then
!  do spectral snapshot
             call teft
             call tmft
             nemw=0
	  else if (nemw.eq.nf) then
             nemw=0
          endif

          nemw=nemw+1


      endif

!  dump time series of backscatter
! back transformations for fields:
!           bzl = bz(i)+vy0*ex(i)
!           bxl = -vy0*ez(i)
!           byl = by(i)/gam0
!           eyl = ey(i)/gam0
!           exl = ex(i)+vy0*bz(i)
!           ezl = ez(i)
      write (70,'(2(1pe16.8))') itime*dt/tconv,bz(1)
      write (71,'(2(1pe16.8))') itime*dt/tconv,ez(1)

      end
