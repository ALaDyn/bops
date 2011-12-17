

!  ==========================================
!
!    Fill remaining octants in velocity space
!
!  ==========================================

subroutine fill3v(i1,n)
  use bopsvars
  implicit none

  integer l, n, ip, i1, iquad

  integer ioc(8),joc(8),koc(8)
  data ioc/1,1,1,1,-1,-1,-1,-1/
  data joc/1,-1,1,-1,1,-1,1,-1/
  data koc/1,1,-1,-1,1,1,-1,-1/

  do l=1,n
     ip=i1+l-1
     iquad=mod(l-1,8)+1
     ux(ip) = ioc(iquad)*ux(ip)
     uy(ip) = joc(iquad)*uy(ip)
     uz(ip) = koc(iquad)*uz(ip)
  end do

end subroutine fill3v
