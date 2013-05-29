


!  =================================
!
!    Electron density gather
!
!  =================================

subroutine eden
  use bopsvars
  implicit none

  integer i, i1, l, i2, iwl, iwr
  real*8 re, xa, f2, f1, rdx, ared
  real*8 rhopriv(0:nx+1)
  rdx = 1./dx

  do i=1,nx+1
     rhopriv(i)=0.
     rhoe(i)=0.
  end do

  !$omp do
  do l=1,ne
     xa=xn(l)*rdx
     re=q(l)/dx
     i1=xa
     i2=i1+1
     f2=xa-i1
     f1=1.-f2
     rhoe(i1+1)=rhoe(i1+1)+re*f1
     rhoe(i2+1)=rhoe(i2+1)+re*f2

 if (i1.lt.0 .or. i1.gt.nx) then
	write (*,*) 'i1:',i1,'l=',l
	write (*,*) 'i2:',i2,'l=',l
 endif
! OMP version
!     rhopriv(i1+1)=rhopriv(i1+1)+re*f1
!     rhopriv(i2+1)=rhopriv(i2+1)+re*f2
  end do


  !  Now do OMP reduce on grid: thread-local copies of rhopriv
  !  contain particle density sums

!  do i=0,nx+1
!     rhoe(i) = rhoe(i) + rhopriv(i)
!  end do

  !      write (*,'(30f8.2)') (rhopriv(i),i=1,nx)

  !      write (*,'(30f8.2)') (rhoe(i),i=1,nx)
  !   periodic bc
  if (ipbc.eq.1) then
     rhoe(1)=rhoe(1)+rhoe(nx+1)
     rhoe(nx+1)=rhoe(1)

     !   reflective
  else if (ipbc.ge.2) then
     iwl=wl/dx+1
     rhoe(iwl+1)=rhoe(iwl+1)+rhoe(iwl)
     rhoe(iwl)=0.

     iwr=wr/dx+1
     rhoe(iwr)=rhoe(iwr)+rhoe(iwr+1)
     rhoe(iwr+1)=0.
  endif
  Qgrid = 0.
  do i=1,nx+1
     Qgrid=Qgrid+ rhoe(i)*dx
     rhot(i)=rhoe(i)+rhoi(i)
  end do
end subroutine eden
