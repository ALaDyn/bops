
!     ========================================================
!
!     Particle tracking initialisation
!
!     $Revision: 26 $
!     ========================================================

subroutine trkini
  use bopsvars
  implicit none

  real*8 x1, xdum
  integer i,l, ntrk

  x1=xpstart
  call r0form(' xpint',xpint,'f12.2')
  call r0form('    x1',x1,'f12.2')



  if (itropt.eq.1) then
     do i=1,ntrack
	itrack(i)=i
     end do

  else if (itropt.eq.2) then
!     track particles with x>=xpstart
     do i=1,ntrack
        l=0
10      l=l+1
        if ((xn(l).gt.x1).and.(xn(l).le.x1+xpint)) then
           itrack(i)=l
           x1=x1+xpint
           call i0prnt('itrack',i)
           call r0form(' posn ',x1,'f12.2')
           goto 100
        endif
        if (l.lt.ne) goto 10
        x1=x1+xpint
100  continue
     end do


     !     read particle labels from file

  else
     call warn('Starting particle tracking   ')
     rewind(11)
     read (11,*) ntrk
     if (ntrk.lt.ntrack) ntrack=ntrk

     do i=1,ntrack
        read(11,'(i8,1pe15.4)') itrack(i),xdum
        ytrk(i,0) = 0.
     end do

  endif
end subroutine trkini









