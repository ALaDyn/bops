
!  ==================================
!
!    Sum electron and ion current currents
!    for em field solver
!
!   $Revision: 6 $
!  ==================================

subroutine addcurr
  use bopsvars
  implicit none
  integer i

  do i=1,nx+1

     !  TE mode

     !   j-
     jm(i) = jyem(i) + jyim(i)
     !   j+
     jp(i) = jye(i) + jyi(i)

     !  TM mode

     !   j-
     jzm(i) = jzem(i) + jzim(i)
     !   j+
     jzp(i) = jze(i) + jzi(i)
  end do

end subroutine addcurr
