
!  =========================
!
!     Time stamp
!
!  =========================

      subroutine stamp(istream)
      character*24 fdate
      character cdate*8, ctime*10, czone*5
      integer vals(4)

      call DATE_AND_TIME(cdate,ctime,czone,vals)

      write(istream,'(a6,a12/a6,a12/a6,a12)') 'BOPS  ' &
      ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
      ,'Time: ',ctime(1:2)//':'//ctime(3:4),'Zone:',czone

!  =========================
!
!     SUN date and time
!
!  =========================

!      write(istream,'(/a24/)') fdate()

      end


