
!  ==========================================================
!
!        Plot exit phase of tracked particles
!
!  ==========================================================

subroutine scatrk
  use bopsvars
  implicit none

  real*8 uxm
  integer i

  uxm=0.

  do i=1,ntrack
     uxm=max(uxm,uest(i))
  end do

  call grps(phat(1:ntrack),uest(1:ntrack),ntrack,0,2*pi,0,uxm,6000,1 &
       ,'     phi       ','       U       ','escp            ')
end subroutine scatrk
