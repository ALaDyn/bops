
!  ======================================
!
!    Electron density
!
!  ======================================


subroutine density

include 'es.h'

 re=qe/dx    !  charge weighting factor
     
 rhoe = 0.   ! zero density array

!  gather density onto grid from particles

 do i=1,ne
      xa = x(i)/dx
      j1 = xa
      j2 = j1+1
      f2 = xa-j1
      f1 = 1.-f2
      rhoe(j1) = rhoe(j1) + re*f1
      rhoe(j2) = rhoe(j2) + re*f2
 end do

 if (bc_field.eq.1) then
!   periodic boundaries 
      rhoe(0) = rhoe(0) + rhoe(nx)
      rhoe(nx) = rhoe(0)

 else if (bc_field.eq.2) then
!   reflective - 1st and last (ghost) cells folded back onto physical grid
	iwl=wall_left/dx
	rhoe(iwl+1)=rhoe(iwl+1)+rhoe(iwl)
        rhoe(iwl)=0.	

	iwr=wall_right/dx
	rhoe(iwr)=rhoe(iwr)+rhoe(iwr+1)
        rhoe(iwr+1)=rhoe(iwr)     

 endif

end
