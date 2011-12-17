!  =======================
!  SLAB - set up multi(currently three) flat slab
!  ======================

subroutine multislab(n,np1,np2,x,x_layer1,x_start,x_finish,x_layer2)
                     ! n : total particle numbers
                     ! np1: the particle for the first proton layer 
                     ! np2: the particle for the second proton layer
                     ! n-np1-np2 : the particle for the ion layer 
  implicit none
  real*8 :: x(n)  ! Array returning particle positions
  real*8, intent(in) :: x_start, x_layer1, x_layer2, x_finish
  integer, intent(in) :: n  ! # total particles
  integer, intent(in) :: np1, np2  ! # particles for proton layers 1 and 2 
                                   ! # front and rear side of the foil   
  real*8 :: dpx
  integer :: i

  if (n.eq.0) return

! set up the 1st proton layer for np1 particles 
  if (np1.gt.0) then
     dpx = (x_start - x_layer1)/np1
     do i=1,np1
        x(i) = x_layer1 + i*dpx - dpx/2.
     end do
  end if

!  set up the 2nd proton layer for np2 particles
  if (np2.gt.0) then 
     dpx = (x_layer2 - x_finish)/np2
     do i=np1+1,np1+np2
        x(i) = x_finish + (i-np1)*dpx - dpx/2.
     end do
  end if

! set up the 3rd ion layer for n-np1-np2 particles
  if ((n-np1-np2).gt.0) then
    dpx = (x_finish - x_start)/(n-np1-np2)
    do i=np1+np2+1,n
       x(i) = x_start + (i-np1-np2)*dpx - dpx/2.
    end do
  end if

end subroutine multislab
