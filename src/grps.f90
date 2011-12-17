
!  =============================================
!
!                      GRPS
!
!   Write out scatter graph
!     in odplot format
!
!  =============================================

      subroutine grps(x,y,n,xmin,xmax,ymin,ymax,id,iskip &
                     ,chx,chy,ctitle)
      implicit none
      integer n,id,iskip
      integer idum, lc, is, np, i, iadd, ipick, ip
      real*8 x(n),y(n),store(15)
      real*4 xp(n),yp(n)
      real*8 xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl, day, dax
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*2
      real*8 rano

!      data idum/-11/
!      save idum
	idum=-111
	
!  open data file
!  2-digit timestamp
     cid(2:2) = achar(mod(id,10) + 48)  
     cid(1:1) = achar(mod(id/10,10) + 48)
!  call chr(mod(id,10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xy'
      open (51,file=cfile)

      if (ymin.eq.ymax) ymax=ymax+1.
      is=iskip
      np=n/is
      if (mod(n,iskip).ne.0) np=np+1

       do i=1,n,is
       iadd=is*rano(idum)+1
        ipick = max(1,min(iadd+i,n))
        ip=(i+iskip-1)/iskip
        xp(ip)=x(ipick)
        yp(ip)=y(ipick)
       end do

      day=(ymax-ymin)/5
      dax=(xmax-xmin)/5

!  write out to 'ODPLOT' datafile

      write (50,101) -id,np,1,0,9,1,chx,chy,ctitle
!  manual axes
      write(50,102) xmin,xmax,dax,ymin,ymax,day
      write (51,103) (xp(i),yp(i),i=1,np)
      close(51)

  101 format (2i6,4i4,2a15,1a16)
  102 format (6(1pe12.3))
  103 format (2(1pe14.4))

      end











