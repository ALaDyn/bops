
!  =========================
!
!  Run data output formatter
!
!  =========================

      subroutine hop(cname,z,ndp)

      implicit none
      real*8 z
      integer ndp, n

      character cname*6,cform*12,cnum*9
      data cnum/'123456789'/
!  printed header
      if (ndp.eq.0) then
      cform='(a6,i13)    '
      write (20,cform) cname,int(z)
      else if (ndp.lt.0) then
       n=-ndp
      cform='(a6,1pe13.'//cnum(n:n)//')'
      write (20,cform) cname,z
      else
      cform='(a6,f13.'//cnum(ndp:ndp)//')  '
      write (20,cform) cname,z
      endif

      end

