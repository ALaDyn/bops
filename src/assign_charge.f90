! Compute electron macro charge according to density profile
!

subroutine assign_charge
use bopsvars
implicit none

real*8 :: xlol, rhomin

!     electron cloud charge normalised to give ncrit=1 (rhoc=-1)
!     all space quantities have been transformed to code frame

      profile: select case(inprof)
 
      case(1)   !     uniform
         if (ipbc.eq.2.or.ipbc.eq.3.or.ipbc.eq.4) then
            xm2=xl-dx*0.5
            xm1=dx*0.5
         else if (ipbc.eq.1) then
            xm1=0.
            xm2=xl
         endif
         wl=xm1
         wr=xm2
         xload=xm2-xm1
         xcrit =xl/2.
         qe=-rho0*xload/ne

      case(2,3)  !     linear/linear with flat top
         xm2=xl-1.5*dx
         wr = xm2
         xm1=max(xm1,wl)
         if (inprof.eq.2) xsol=xm2
         qe=-rho0/ne*(xm2-xsol/2-xm1/2)
         xcrit = xm1 + 2*pi*xlolam/gam0

      case(4)  !     trapezoidal
         wr = xl-1.5*dx
         wrf = wr
         qe=-rho0/2/ne*(xm2+xsol2-xsol-xm1)
         xcrit = xm1 + 2*pi*xlolam/gam0

      case(5)  !     exponential
         wr = xl-1.5*dx
         xm2 = wr
         xsol2 = wr
	 rhomin=0.1  ! 1/10 nc - TODO: probably want to scale this by gam0**3
         xlol = 2*pi*xlolam/gam0
!  Adjust start of exp profile to begin at min. density above
	 xm1 = xsol - xlol*log(rho0lab/rhomin)
         xload = (xsol - xm1)
         qe = -rho0/ne*(xlol*(1.-exp(-xload/xlol)) + xm2 - xsol)
         xcrit = xsol-xlol*log(rho0lab)

      case(6)  !     tanh
         xm2=xsol2
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xm2-xm1/2.-xsol/2.)
         xcrit = xsol-xlol*log(rho0lab)

      case(7)  !     foil target in middle of grid
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*1.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

      case(8)  !     layered
         xsol2=xl-dx/2.
         qe=-(rho_layer*(xm2-xm1) + rho0*(xsol2-xsol))/ne
         xcrit = xm1

      case(9)  ! 2 or 3-layer targets
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

        target: select case(target_config)

         case(3)        ! three layers with proton layer on both front and rear sides
           n_pp=nint(-rho_layer*2*x_layer/qe) ! add particle numbers for 
                                            ! both electrons and prtons in the proton layers   
           n_pp1 = nint(n_pp/2.0)		   ! the particle number (electrons and ions) for 1st proton layer 
           n_pp2=n_pp-n_pp1                ! the particle number (electrons and ions) for 2nd proton layer
         
         case(4)     ! double layers with proton layer on the rear side   
           n_pp=nint(-rho_layer*x_layer/qe)  ! add particle numbers for
                                            ! both electrons and prtons in the proton layers
           n_pp1 = 0                        ! the particle number (electrons and ions) for 1st proton layer
           n_pp2=n_pp                       ! the particle number (electrons and ions) for 2nd proton layer
      
         case(5)     ! double layers with proton layer on the front side
           n_pp=nint(-rho_layer*x_layer/qe) ! add particle numbers for
                                            ! both electrons and prtons in the proton layers
           n_pp1 = n_pp                     ! the particle number (electrons and ions) for 1st proton layer
           n_pp2=n_pp-n_pp1                 ! the particle number (electrons and ions) for 2nd proton layer

         case default
            ! do nothing
         end select target

      case(10) ! multi-species target
         xm2=xm1+dfoil
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=xm2-xm1
         qe=-rho0*xload/ne
         xcrit = xm1

         n_pp=nint(-rho_layer*xload/qe)       !add particle numbers for
                                              !both electrons and prtons for proton species
         n_pp1 = n_pp
         n_pp2=n_pp-n_pp1

      case (17)  !    picket-fence 'foam' target in middle of grid
!   nfoil foils of thickness dfoil, spacing gap_foil
         xm2=xm1+dfoil*nfoil + gap_foil*(nfoil-1)
         wl=dx+0.5
         wr=xl-dx*0.5
         xload=dfoil
         qe=-rho0*xload/ne*nfoil
         xcrit = xm1

      case(57)  !     foil target in middle of grid with exponential front side
         wr = xl-1.5*dx
         wl=dx+0.5
         xm2=xsol+dfoil
         xsol2=xm2
 	 wr=xm2
	 rhomin=0.1*gam0**3  ! 1/10 nc
         xlol = 2*pi*xlolam/gam0
!  Adjust start of exp profile to begin at min. density above
!	 xm1 = xsol - xlol*log(rho_layer/rhomin)
         xload = (xsol - xm1)
         qe = -rho0/ne*dfoil - rho_layer/ne*xlol*(1.-exp(-xload/xlol))
         xcrit = xsol-xlol*log(rho0lab) ! TODO: needs adjusting for rho_layer<> rho0
         inprof = 5  ! reset profile flag for denprof routine

      case(575)  !  'soft' foil target in middle of grid with exponential front and back side
         wr = xl-1.5*dx
         wl=dx+0.5
         xsol2=xsol+dfoil
         xload = (xsol - xm1)   ! ramp length
         xm2 = xsol2 + xload
         xlol = 2*pi*xlolam/gam0
         qe = -rho0/ne*(xlol*(2.-2.*exp(-xload/xlol)) + xsol2 - xsol)
         xcrit = xsol-xlol*log(rho0lab)
         inprof = 5  ! reset profile flag for denprof routine

      case default
         wl=xm1
         wr=xm2
         xload=xm2-xm1
         xcrit =xl/2.
         qe=-rho0*xload/ne
      end select profile

!     position of first ion in particle arrays
      ion1=ne+1

!     NB: qe(sim) = qe(lab)*gam0**2.  Same for me,mi etc

! Define particle macro charges and masses
! qome=-1     already set in init

!#### Anupam & Bin 2009/2010 : 
      qi = -z*qe	! ion macro-charge
!     to calculate the density rho0_i=ni*qi=z*ni*qe=z*ne*qe
!     this needs to be divided by z

      qp = -qe   	! proton macro-charge
      me = qe*qome	! electron macro-mass (qome=-1)
      mi = me*miome	! ion macro-mass
      mp = me*mpome	! proton macro-mass
      qomi = qi/mi	! Zm_e/Am_p
      qomp = qp/mp	! q_p/m_p

!#### Anupam & Bin 2009/2010 :Total # simulation particles
      ne=n_pp+ne        ! in total electrons add particles for proton layers
      np=n_pp           ! as in bopsvars.f90, np is the particle number for proton layers 
      npart = ne+ni+np  ! for ni, the charge is qi=-z*qe; for np, the charge is qp=-qe
      ni_tot=ni+np      ! total # ions

!#### Anupam & Bin 2009/2010: position of first ion in particle arrays
!     change calculated ion1 after adding the particles for proton layer 
      ion1=ne+1          ! It is very IMPORTANT  !!!!!!!!!

!#### Anupam & Bin 2009/2010

!     initial drift energies - thermal 'zero'
      uthe0=ne*(gam0-1.d0)*me  ! electrons
      uthi0=ni*(gam0-1.d0)*mi  ! heavy ions
      uthp0=np*(gam0-1.d0)*mp  ! protons

! -------------------------------------------

end
