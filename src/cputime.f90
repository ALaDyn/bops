
      subroutine cputime(t,i1)
      implicit none

      real t, tar(2), timef, tuser
      real etime
      integer i1

!  IBM AIX
       call cpu_time(t)
!  timef returns elapsed wall-clock time in ms
!      t=1.e-3*timef()
!  Linux i386, GNU f77
!      tuser=etime(tar)
!      t=tar(1)

!  Cray t3e
!       call second(t1)
!       t=t1
      end
