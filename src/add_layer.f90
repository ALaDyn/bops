!     ======================================================
!
!     Create additional ion layer with different charge, density
!
!     ======================================================

      subroutine add_layer
      use bopsvars
      implicit none

      integer :: p,i, nprot, ni_layer, ratio_i

! Initial profile has ni ions, ne=Z.ni electrons

      target: select case(target_config)

      case(2)
! Create proton layer on rear of foil target at same total density
        nprot = 0
	qp = -qe  ! Proton and electron charge equal 
	np = rho_layer*x_layer/qp  ! # protons
        ni_layer = nonc*x_layer/qi  ! # heavy ions in layer
	ratio_i = ni_layer/np  ! ratio of ions:protons in layer (>1)

	do i=1,ni,ratio_i
	   p=ne+i
	   if (q(p)>0 .and. xo(p).ge.(xm2-x_layer) .and. xo(p).le.xm2) then
	      nprot=nprot+1
	      q(p) = qp   ! convert ion to proton
	      m(p) = mp
	      species(p) = 3  ! proton label
	   endif
        end do

      case(3)
! Create proton layer on front of foil target at same total density 
        nprot = 0
	qp = -qe  ! Proton and electron charge equal 
	np = rho_layer*x_layer/qp  ! # protons
        ni_layer = nonc*x_layer/qi  ! # heavy ions in layer
	ratio_i = ni_layer/np  ! ratio of ions:protons in layer (>1)

	do i=1,ni,ratio_i
	   p=ne+i
	   if (q(p)>0 .and. xo(p).ge.xm1 .and. xo(p).le.xm1+x_layer) then
	      nprot=nprot+1
	      q(p) = qp   ! convert ion to proton
	      m(p) = mp
	      species(p) = 3  ! proton label
	   endif
        end do

      case default
       ! do nothing 
      end select target

     write (*,*) 'Adding proton layer ..'
     write (*,*) 'Check: nprot= ',nprot,'np= ',np,' ni_layer=',ni_layer
     write (*,*) 'Ion/proton ratio',ratio_i
      if (target_config==2) then
        write (*,*) 'Limits',xm2-x_layer,'-',xm2
      else 
        write (*,*) 'Limits',xm1,'-',xm1+x_layer
      endif
     write (*,*) 'Density: ',rho_layer
     write (*,*) 'Total charge: ',rho_layer*x_layer

      write (*,*) 'Charges: ',qe, qi, qp
      write (*,*) 'Masses: ',me, mi, mp

      np = nprot
      ni = ni_tot - np      ! # heavy ions
      write (*,*) 'Heavy ions ni= ',ni
      write (*,*) 'Total ions = ',ni_tot

      end subroutine add_layer











