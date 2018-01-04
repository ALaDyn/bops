!
!   FINDXC
!
!   routine to return position of density surface (den=rhotrack)
!  ==========================================

 subroutine findxc(xc1,xc2,den,nx,dx,rhotrak)

        implicit none
        real*8 xc1, xc2, dx, rhotrak
        integer nx
        real*8 :: den(0:nx+1)
        real*8 xx,denmin
        integer icm, icp, i
        logical found

!  start indices
      icm = xc1/dx
      icp = icm+1
      found = .false.

    denmin = minval(den(1:nx))
    if (denmin.gt.1) then
	write(*,*) 'Profile has no critical density'
	return
    endif

!  sweep up
      i = min(icm,nx-1)
      do while (i.lt.nx .and. .not.found)
         if (den(i).le.rhotrak .and. den(i+1).gt.rhotrak) then
            found=.true.
                xx = (i-1)*dx
            xc2 = xx + dx*(rhotrak-den(i))/(den(i+1)-den(i))
         endif
         i = i+1
      end do

!  sweep down
      i = min(icm,nx-1)
      do while (i.gt.0 .and. .not.found)
         if (den(i).le.rhotrak .and. den(i+1).gt.rhotrak) then
            found=.true.
                xx = (i-1)*dx
            xc2 = xx + dx*(rhotrak-den(i))/(den(i+1)-den(i))
         endif
         i = i-1
      end do


      if (.not.found) xc2=xc1




      end
