
!  =============================
!
!    Boost sim frame -> lab frame
!
!  =============================

      subroutine sim2lab(v0,g0,uysim,gsim,uylab,glab)
      implicit none
      real*8 v0,g0,uylab,glab,uysim,gsim

      uylab = g0*( uysim - v0*gsim )
      glab = g0*( gsim - v0*uysim )

      end
