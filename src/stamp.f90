
!     =========================
!
!     Time stamp
!
!     =========================

      subroutine stamp(istream,ibegin)
      implicit none

      character fdate*24
      character cdate*8, ctime*10, czone*5
      integer vals(4)
      integer ibegin, istream

!      call DATE_AND_TIME(cdate,ctime,czone,vals)
       call DATE_AND_TIME(cdate,ctime,czone)

      if (ibegin.eq.1) then

         write(istream,'(//a6,a12/a6,a12/a6,a12)') 'BOPS  ' &
              ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
              ,'Time:',ctime(1:2)//':'//ctime(3:4),'Zone:',czone

      else
         write(istream,'(a,a6)') 'Finished run at time: ' &
              ,ctime(1:2)//':'//ctime(3:4)
      endif
!     =========================
!
!     SUN date and time
!
!     =========================

!     write(istream,'(/a24/)') fdate()

      end


