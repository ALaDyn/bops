
!  ==================================
!
!    Electrostatic field
!      - direct integration
!
!  ==================================

subroutine esfield
  use bopsvars
  implicit none

  integer i

  !  net charge density - rho

  do i=1,nx+1
!#### Anupam & Bin 2009/2010
!    for carbon ions: should multiply by z as charge in the heavy ion is z
     rhot(i)=rhoe(i)+rhoi(i)*z+rhop(i)    
!#### Anupam & Bin 2009/2010
  end do

  ! filter net charge density
!  call filter1(rhot,nx+1)

  !   integrate div.e=rho directly (trapezium approx)


  !   end point - ex=0 mirror at particle boundary wr

  !        iwr = wr/dx + 1
  !        ex(iwr+1) = 0.5*dx*rhot(iwr)

  !        do i=iwr,1,-1
  !          ex(i)=ex(i+1)-0.5*dx*(rhot(i)+rhot(i+1))
  !        end do

  !        do i=iwr+2,nx
  !          ex(i)=0.
  !        end do

  ex(nx+1) = 0.
  do i=nx,1,-1
     ex(i)=ex(i+1)-0.5*dx*(rhot(i)+rhot(i+1))
  end do


  !   potential - same again
  ! 21.4.2011 fixed sign consistent with E=-d(phi)/dx

  phi(nx+1)=0.

  do i=nx,1,-1
     phi(i)=phi(i+1)+0.5*(ex(i)+ex(i+1))*dx
  end do



end subroutine esfield
