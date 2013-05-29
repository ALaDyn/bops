
!  ============================================
!
!    Direct field integration - slab geometry
!
!  ============================================

subroutine field

include 'es.h'

real :: rhot(0:nxm)     !  net density

!  Add neutral background to get net density
 
 rhot(1:nx) = rhoe(1:nx) + rho0
      

!   integrate div.E=rho directly (trapezium approx)


!   end point - ex=0 mirror at right wall
      Ex(nx+1)=0.
      edc = 0.

!  integrate from right to left

      do j=nx,1,-1
         Ex(j) = Ex(j+1) - 0.5*(rhot(j) + rhot(j+1))*dx
         edc = edc + Ex(j)
      end do


 if (bc_field.eq.1) then

!  periodic fields - remove DC component 
!  -- need this for consistency with charge conservation
      
      do i=1,nx
         ex(i)=ex(i) - edc/nx
      end do
      
!  end points periodic
      ex(0) = ex(nx)

 endif

!   potential - same again
!   - integration const. is arbitrary - phi not used for anything at present

      phi(nx+1)=0.

      do j=nx,1,-1
         phi(j) = phi(j+1) + 0.5*(Ex(j) + Ex(j+1))*dx
      end do 

 
end
