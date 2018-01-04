
!  =================================
!
!   Write xy graph in odplot format
!
!  =================================

      subroutine grxy(xx,yy,nd,id,igx,iman,chx,chy,ctitle &
                     ,xmin,xmax,dax,ymin,ymax,day)

      implicit none

      integer, intent(in) :: nd
      integer id, igx, iman
      real*8, intent(in) :: xx(0:nd),yy(0:nd)
      real*4 :: x(nd),y(nd)
      real*8 :: store(15), xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,dax,day
      integer lc, i, n, ip, idp, ispage, idash, ityp, ibox
      character chars(1)*8,chx*15,chy*15,cfile*12,cid*2
      character(len=*), intent(in) :: ctitle
!  open data file
!  2-digit timestamp
     cid(2:2) = achar(mod(id,10) + 48)  
     cid(1:1) = achar(mod(id/10,10) + 48)
!      call chr(mod(id,10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xy'
      write (15,'(i6,1x,a)') id,cfile
      open (51,file=cfile)

      if (igx.eq.0) igx=1
      n=nd/igx
      if (mod(nd,igx).ne.0) n=n+1

      do i=1,nd,igx
      ip=(i+igx-1)/igx
      x(ip)=xx(i)
        y(ip)=yy(i)
      end do

      idp=id
      ispage=1
      idash=1

      if (iman.lt.0) then
!  manual axes: ityp=-iman; id=-id
      ityp=-iman
      idp=-id

      else if (iman.eq.0) then
!  plot using old set of axes
      ityp=1
      ispage=2
      idash=2

      else if (iman.gt.0) then
!  auto axes: ityp=iman
      ityp=iman
      endif

!  write out to 'ODPLOT' datafile

      write (50,101)idp,n,ityp,idash,0,ispage,chx,chy,ctitle
!  manual axes
      if (iman.lt.0) then
      write(50,102) xmin,xmax,dax,ymin,ymax,day
      endif
      write (51,103) (x(i),y(i),i=1,n)
      close(51)

  101 format (2i6,4i4,2a15,1a16)
  102 format (6(1pe12.3))
  103 format (2(1pe15.5))

      end
