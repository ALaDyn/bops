
!  =============================================
!
!                      GRPS
!
!   Write out scatter graph
!     in odplot format
!
!  =============================================

      subroutine gr3d(x,y,z,i1,n,iskip,id,ctitle &
      ,xmin,xmax,ymin,ymax,zmin,zmax,xcut,vcut)

      implicit none
      integer i1, n, iskip, id

      real*8 x(0:n),y(0:n),z(0:n),store(15),xp(40000),yp(40000) &
      ,zp(40000)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,x1,x2,y1,y2,xl,yl &
      ,xcut,vcut, v2
      character chars(1)*8,ctitle*16,cfile*17,cid*3
      integer i, ip, idum, np, idp, lc
      data idum/-11/
      save idum

!  open data file
      call chr(id*1.0,0,cid,lc)
      if (lc.eq.1) then
        cfile=ctitle(1:10)//"_00"//cid
      else if (lc.eq.2) then
        cfile=ctitle(1:10)//"_0"//cid
      else
        cfile=ctitle(1:10)//"_"//cid
      endif
      write(*,'(2a)') 'writing movie file ',cfile
      open (51,file=cfile)
!  header
      write (51,105) xmin,xmax,ymin,ymax,zmin,zmax
  105 format('#define xmin',f12.1/ &
            '#define xmax',f12.1/ &
            '#define ymin',f12.1/ &
        '#define ymax',f12.1/ &
        '#define zmin',f12.1/ &
        '#define zmax',f12.1)

      ip=0
      do i=i1,i1+n-1,iskip
      v2 = y(i)**2+z(i)**2
      if (x(i).lt.xcut.or.v2.gt.vcut**2) then
          ip=ip+1
        xp(ip)=x(i)
          yp(ip)=y(i)
          zp(ip)=z(i)
      endif
      end do
      np=ip
      idp=-id

!  write out in 'rayshade' format
      do i=1,np
      if (i1.eq.1) then
        write (51,103) xp(i),yp(i),zp(i)
      else
        write (51,104) xp(i),yp(i),zp(i)
      endif
      end do
      close(51)

  103 format ('elec(',1pe9.3,',',1pe10.3,',',1pe10.3,')')
  104 format ('ion(',1pe9.3,',',1pe10.3,',',1pe10.3,')')

      end
