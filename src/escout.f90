
!  =======================================================
!
!     Pick particles to be tracked u**2>uhot
!
!  =======================================================

subroutine escout
  use bopsvars
  implicit none

  integer i
  real*8 du

  do i=1,5
     write(90+i,*) nesc
  end do
  du=a0**2/2.

  !  do approx sort into energies
  do i=1,nesc
     if (uesc(i).ge.uhot .and. uesc(i).lt.uhot+du) then
        write(90,101) iesc(i),uesc(i)
     else if (uesc(i).ge.uhot+du &
          .and. uesc(i).lt.uhot+2*du) then
        write(91,101) iesc(i),uesc(i)
     else if (uesc(i).ge.uhot+2*du &
          .and. uesc(i).lt.uhot+3*du) then
        write(92,101) iesc(i),uesc(i)
     else if (uesc(i).ge.uhot+3*du &
          .and. uesc(i).lt.uhot+4*du) then
        write(93,101) iesc(i),uesc(i)
     else if (uesc(i).ge.uhot+4*du &
          .and. uesc(i).lt.uhot+5*du) then
        write(94,101) iesc(i),uesc(i)
     else if (uesc(i).ge.uhot+5*du) then
        write(95,101) iesc(i),uesc(i)

     endif
101  format (i8,1pe15.4)

  end do

end subroutine escout
