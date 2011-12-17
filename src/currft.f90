!  =======================================================
!
!       2D FT of current source
!
!  =======================================================

      subroutine currft
      use bopsvars
      implicit none

      integer jcou, ncont, ncy, icend, i1, iskip, ix, i, isol, j
      integer i3, i4, i2, ncx, nfom, k
      real*8 tspec, dts, dwc, wcmax, pha1, pha, pha2, pha3, pha4

      real*8 cur2d(nx,ncycle*itav)
      real :: rw1(nf), rw2(6*nf+150)  ! FT work arrays
      integer :: iwk(6*nf+150)    !

      complex cx(nf/2+1)
      real*8 j1amp(0:nx+1),j2amp(0:nx+1),j3amp(0:nx+1),j4amp(0:nx+1) &
          ,j1pha(0:nx+1),j2pha(0:nx+1),j3pha(0:nx+1),j4pha(0:nx+1)
      real*8 :: work1(0:nf+1), work2(0:nf+1)
      integer :: ncym, ncxm
      save jcou,ncont
      data jcou,ncont/1,1/


      if (itime.eq.0) return
!  ncy should be power of 2
!  itav = isgy/ncycle*2**a

      ncym = itav
      ncy=ncym
      ncxm = nx
      icend = icstart+ncycle*itav


!  accumulate 2D lab frame y-current (x,y)
!  sampled over total time  ncycle*itav*dt
      if (itime.ge.icstart .and. jcou.le.ncym) then
!        call filter1(jye,nx)
        i1 = xcur1/dx
        i2 = xcur2/dx
        if ((i2-i1).lt.ncxm) then
          ncx=i2-i1
          iskip=1
        else
          iskip = (i2-i1)/ncxm+1
          ncx=(i2-i1)/iskip
        endif

        ix = 0
        do i=i1,i2,iskip
          ix=ix+1
        cur2d(ix,jcou)=(jye(i)-vy0*rhoe(i))/gam0
        end do
        jcou=jcou+1
      endif

!  fourier transform current in y
      if (jcou.gt.ncym) then
!  reset sample start time
        icstart=icstart+ncycle*itav

        tspec = icft*ncym*dt
        dts = dt*icft
        dwc = 2*pi/tspec
        wcmax = pi/dts
        if (ncont.eq.1) then
          call blank
          write(15,*) 'Current ft'
          write(15,*) 'dw = ',dwc
          write(15,*) 'wmax = ',wcmax
          call blank
!  check current ft at xsol
          isol=(xsol-xcur1)/dx/iskip
          do j=1,ncy
            rw1(j)=cur2d(isol,j)
          end do
          call fftrc(rw1,ncy,cx,iwk,rw2)
          do j=1,ncy/2
            work2(j)=abs(cx(j+1))/ncy
            work1(j)=dwc*j
          end do
          call grxy(work1,work2,ncy/2,80000+ncont,1,1 &
      ,'     w          ',' Jy(w) at xsol ','jyxs'//ctime(1:12))
!  check current ft at xsol+c/w0
          isol=(xm2-xcur1)/dx/iskip
          do j=1,ncy
            rw1(j)=cur2d(isol,j)
          end do
          call fftrc(rw1,ncy,cx,iwk,rw2)
          do j=1,ncy/2
            work2(j)=abs(cx(j+1))/ncy
            work1(j)=dwc*j
          end do
          call grxy(work1,work2,ncy/2,80500+ncont,1,1 &
      ,'     w          ',' Jy(w) at xsol2','jyxp'//ctime(1:12))
        endif

        i1 = 1./dwc+1
        i2 = 2./dwc+1
        i3 = 3./dwc+1
        i4 = 4./dwc+1
        write (15,*) i1,i2,i3,i4

        open(51,file='currft.dat')

        do i=1,ncx
          do j=1,ncy
            rw1(j) = cur2d(i,j)
          end do
          call fftrc(rw1,ncy,cx,iwk,rw2)
          nfo2=ncy/2
!  Average over 3*dwc
          j1amp(i) = (abs(cx(i1-1))+abs(cx(i1))+abs(cx(i1+1)) )/3./ncy
          j2amp(i) = (abs(cx(i2-1))+abs(cx(i2))+abs(cx(i2+1)) )/3./ncy
          j3amp(i) = (abs(cx(i3-1))+abs(cx(i3))+abs(cx(i3+1)) )/3./ncy
          j4amp(i) = (abs(cx(i4-1))+abs(cx(i4))+abs(cx(i4+1)) )/3./ncy
          pha1 = ( pha(cx(i1-1))+pha(cx(i1))+pha(cx(i1+1)) )/3.
          pha2 = ( pha(cx(i2-1))+pha(cx(i2))+pha(cx(i2+1)) )/3.
          pha3 = ( pha(cx(i3-1))+pha(cx(i3))+pha(cx(i3+1)) )/3.
          pha4 = ( pha(cx(i4-1))+pha(cx(i4))+pha(cx(i4+1)) )/3.
          j1pha(i) = mod(pha1+pi,2*pi)-pi
          j2pha(i) = mod(pha2+pi,2*pi)-pi
          j3pha(i) = mod(pha3+pi,2*pi)-pi
          j4pha(i) = mod(pha4+pi,2*pi)-pi
!  space coord
          work2(i) = xcur1*gam0+(i-1)*iskip*dx*gam0
!  o/p whole 2d array: omega-x space
          nfom = omegm/dwc
          do k=1,nfom
            work1(k) = abs(cx(k+1))/ncy
          end do
!  smooth in w-space
          call filter1(work1,nfom)
          write (51,'(3(1pe15.4))') (work2(i),k*dwc,work1(k),k=1,nfom)
          write (51,'(a)')
        end do

!  plot currents
       call grxy(work2,j1amp,ncx,81000+ncont,1,1 &
      ,'     x          ',' J_w           ','jy1w'//ctime(1:12))
       call grxy(work2,j1pha,ncx,81500+ncont,1,1 &
      ,'     x          ',' Jpha_w        ','jp1w'//ctime(1:12))
       call grxy(work2,j2amp,ncx,82000+ncont,1,1 &
      ,'     x          ',' J_2w          ','jy2w'//ctime(1:12))
!     :,xcur1,xcur2,dom,uwmin,uwmax,10.)
       call grxy(work2,j2pha,ncx,82500+ncont,1,1 &
      ,'     x          ',' Jpha_2w       ','jp2w'//ctime(1:12))
       call grxy(work2,j3amp,ncx,83000+ncont,1,1 &
      ,'     x          ',' J_3w          ','jy3w'//ctime(1:12))
       call grxy(work2,j3pha,ncx,83500+ncont,1,1 &
      ,'     x          ',' Jpha_3w       ','jp3w'//ctime(1:12))
       call grxy(work2,j4amp,ncx,84000+ncont,1,1 &
      ,'     x          ',' J_4w          ','jy4w'//ctime(1:12))
        ncont=ncont+1
        jcou=1
        close(51)
      endif
      end
