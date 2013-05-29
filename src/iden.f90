
!  =================================
!
!    Heavy ion and proton density gather
!
!  =================================

subroutine iden(ip1)
  use bopsvars
  implicit none

  integer i, i1, l, i2, iwl, iwr, ip, ip1
  real*8 ri, xa, f2, f1, rdx

  rdx = 1./dx

  if (ni_tot.eq.0) return


  do i=1,nx+1
     rhoi(i)=0.
     rhop(i)=0.
  end do

  do l=1,ni_tot
     ip=ip1+l-1
     xa=xn(ip)*rdx
     i1=xa
     i2=i1+1  
     f2=xa-i1  ! i1<xa<i2; f2 prop to distance between particle and furthest gp
     f1=1.-f2
     ri=q(ip)/dx  ! particle charge density

     if (species(ip).eq.2) then
       rhoi(i1+1)=rhoi(i1+1)+ri*f1   ! Index offset by 1, since rho(1) at x=0
       rhoi(i2+1)=rhoi(i2+1)+ri*f2
     else
! proton layer
       rhop(i1+1)=rhop(i1+1)+ri*f1
       rhop(i2+1)=rhop(i2+1)+ri*f2
     endif
  end do

  if (ipbc.eq.1) then
     !   periodic bc
     rhoi(1)=rhoi(1)+rhoi(nx+1)
     rhoi(nx+1)=rhoi(1)
     rhop(1)=rhop(1)+rhop(nx+1)
     rhop(nx+1)=rhop(1)
  else if (ipbc.ge.2) then
     !   reflective at particle walls wl and wr
     iwl=wl/dx+1
     rhoi(iwl+1)=rhoi(iwl+1)+rhoi(iwl)
     rhop(iwl+1)=rhop(iwl+1)+rhop(iwl)
    !      rhoi(iwl)=rhoi(iwl+1)
     rhoi(iwl)=0.
     rhop(iwl)=0.
     iwr=wr/dx+1
     rhoi(iwr)=rhoi(iwr)+rhoi(iwr+1)
     rhoi(iwr+1)=0.
     iwr=wr/dx+1
     rhop(iwr)=rhop(iwr)+rhop(iwr+1)
     rhop(iwr+1)=0.
     !      rhoi(iwr+1)=rhoi(iwr)
  endif

  do i=1,nx+1
!#### Anupam & Bin 2009/2010
     rhoi(i) = rhoi(i)/z      ! for ions, its density should divide by its charge z 
                              ! and initially rho_i=rho_e/z
!#### Anupam & Bin 2009/2010

     rhot(i) = rhoe(i)+rhoi(i)+rhop(i)
!     if (abs(rhot(i)).lt.1.e-8) rhot(i)=0.
  end do

end subroutine iden
