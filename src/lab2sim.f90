
!  =============================
!
!    Boost lab frame -> sim frame
!
!  =============================

      subroutine lab2sim(v0,g0,uylab,glab,uysim,gsim)
      implicit none
      real*8 v0,g0,uylab,glab,uysim,gsim

      uysim = g0*( uylab + v0*glab )
      gsim = g0*( glab + v0*uylab )

      end
