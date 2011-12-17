
!  =================================
!
!   Write xy particle tracking data in cplot format
!
!  =================================

      subroutine gtrak(x,y,n,id,igx,iman,chx,chy,ctitle &
                     ,xmin,xmax,dax,ymin,ymax,day)

      implicit none


      integer iman, id, igx, n
      real*8 x(0:n),y(0:n)
      real*8 xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,dax,day
      integer lc, i
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*3

!  open data file
      call chr(id*1.0,0,cid,lc)
      cfile=ctitle(1:4)//'/'//cid(1:lc)
      write (15,'(i6,1x,a)') id,cfile

      open (51,file=cfile)

      write (51,103) (x(i),y(i),i=1,n)
      close(51)

  103 format (2(1pe12.4))

      end
