
!  ===========================================
!
!  Load particles according to density profile
!
!  ===========================================

subroutine denprof(i1,n,charge)
  use bopsvars
  implicit none
  integer, intent(in) :: n  ! # parts
  integer, intent(in) :: i1 ! # start index
  real*8, intent(in) :: charge ! particle macro-charge
  integer :: l, nstrip, n1, n2, n3, i, ip, ncheck
  integer :: npart_foil, istart, j
  real*8 :: nmin
  real*8 :: xplas, xlol, x, xc, a, xprime, xload1, xload2
  real*8 :: den,dxs,dmax, xstart

  if (n.eq.0) return

!   uniform profiles

!  whole grid
  if (inprof.eq.1) then
     call xq(i1,n,xm1,xload)

! foil
  else if (inprof.eq.7) then
!     call xq(i1,n,xm1,dfoil)

!#### Anupam & Bin 2009/2010
     call slab(n,xo(i1:i1+n-1),xm1,xm1+dfoil)   ! changed i1:n+1 to i1:n

!#### Anupam & Bin 2009/2010: profile for three layered target

  else if (inprof.eq.9) then
     call multislab(n,n_pp1,n_pp2,xo(i1:i1+n-1),xm1-x_layer,xm1,xm1+dfoil,xm1+dfoil+x_layer)

!#### Anupam & Bin 2010 : profile for multi-species target

  else if (inprof.eq.10) then
     	if (np .gt. 0) call slab(np,xo(i1:i1+np-1),xm1,xm1+dfoil)
     	if (n-np .gt. 0) call slab(n-np,xo(i1+np:i1+n-1),xm1,xm1+dfoil)

!#### Anupam & Bin 2009/2010
   
! picket fence
  else if (inprof.eq.17) then
    npart_foil = n/nfoil ! # particles per foil
    do j=1,nfoil 
	xstart = xm1 + (j-1)*(dfoil+gap_foil)
	istart = i1 + (j-1)*npart_foil  ! net particle offset
    	call xq(istart,npart_foil,xstart,dfoil)
!    	write (*,'(i8,f15.3)') (i,xn(i),i=istart,istart+npart_foil-1)
    end do	
 endif

  xlol = 2*pi*xlolam/gam0

  !  First find number of particles in each part of profile

  !   linear

  if (inprof.eq.2) then
     n1=n
     xsol2=xm2
     xsol=xm2

     !  linear with flat top

  else if (inprof.eq.3) then
     n1=0.5*rho0/(abs(charge))*(xsol-xm1)
     n2=n-n1
     xsol2=xm2

     !  trapezoidal

  else if (inprof.eq.4) then
     n1=0.5*rho0/(abs(charge))*(xsol-xm1)
     n3=0.5*rho0/(abs(charge))*(xm2-xsol2)
     n2=n-n1-n3

     !  exponential

  else if (inprof.eq.5) then
     xload1=xsol-xm1
     n1 = xlol*rho_layer/abs(charge)*(1.-exp(-xload1/xlol))
     nmin = rho_layer*exp(-xload1/xlol)
     if (n1.lt.0) then
	write(*,*) 'Something wrong with density profile setup: check xm1, xsol, xlolam'
     else
	write(*,*) 'Min density at xm1: ',nmin
     endif
     xload2=xm2-xsol2
     if (xload2 .gt. 0) then
        n3 = xlol*rho_layer/abs(charge)*(1.-exp(-xload2/xlol))
     else
        n3=0.
     endif

     ncheck = rho0/abs(charge)*(xsol2-xsol)
     n2 = n-(n1+n3)
     call i0prnt('n1    ',n1)
     call i0prnt('n3    ',n3)
     call i0prnt('n2    ',n2)
     call i0prnt('ncheck',ncheck)
     call i0prnt('n1+2+3',n1+n2+n3)


     !  tanh

  else if (inprof.eq.6) then
     n2 = rho0/abs(charge)*(xm2-xsol)
     n1 = n-n2
     xsol2 = xm2


  endif



  if (inprof.ge.2.and.inprof.le.4.and.n1.gt.0) then

     !   linear with moat: add first ramp

     den=0.d0
     xplas=xsol-xm1
     dmax=2.*n1/xplas**2
     l=1
     nstrip=n1*8
     dxs=xplas/nstrip

     do i=1,nstrip+1
        x=xm1+i*dxs-dxs/2.
        den=den+dxs*(x-xm1)*dmax
        if (den.ge.l) then
           xo(i1-1+l)=x
           l=l+1
        endif
     end do

  else if (inprof.eq.5 .and. n1.gt.0) then

     !  exponential ramp

     xplas=xsol-xm1

     do i=1,n1

        x = -xlol*log(1.-1.*i/n1*(1.-exp(-xplas/xlol)))
        xo(i1-1+i) = xsol-x

     end do

  else if (inprof.eq.6 .and. n1.gt.0) then

     !  tanh ramp

     xplas=xsol-xm1
     xc = xm1+xplas/2.
     xlol = 2*pi*xlolam/gam0
     a = 2/xlol/rho0*gam0**3
     nstrip=n1*8
     dmax = n1/xplas
     dxs=xplas/nstrip
     den=0.d0
     l=1

     do i=1,nstrip+1
        x=xm1+i*dxs-dxs/2.
        xprime=x-xc
        den=den+dxs*dmax*(1.+tanh(a*xprime))
        if (den.ge.l) then
           xo(i1-1+l)=x
           l=l+1
        endif
     end do


  endif

  !   add flat slab if not all particles loaded yet

  if (inprof.ge.3.and.inprof.le.6.and.n2.gt.0) then
     dxs=(xsol2-xsol)/n2
     do l=1,n2
        ip=i1-1+l+n1
        xo(ip)=xsol+l*dxs-dxs/2.
     end do
  endif

  !  add trailing ramp

  if (inprof.eq.4.and.n3.gt.0) then
     !  linear
     xplas=xm2-xsol2
     den=0.
     dmax=2.*n3/xplas**2
     l=1
     nstrip=n3*2
     dxs=xplas/nstrip

     do i=1,nstrip+1
        x=xm2-i*dxs+dxs/2.
        den=den+dxs*(xm2-x)*dmax
        if (den.ge.l) then
           xo(n+i1-l)=x
           l=l+1
        endif
     end do

  else if (inprof.eq.5.and.n3.gt.0) then
     ! exponential
     xplas=xm2-xsol2
     xlol = 2*pi*xlolam/gam0

     do i=1,n3

        x = -xlol*log(1.-1.*i/n3*(1.-exp(-xplas/xlol)))
        xo(i1-1+i+n1+n2) = xsol2+x

     end do
  endif

  !   layered target

  if (inprof.eq.8) then
     n1=rho_layer/abs(charge)*(xm2-xm1)
     n2=n-n1
     write (6,*) n1,n2,xm1,xm2,xsol,xsol2
     write (15,*) n1,n2,xm1,xm2,xsol,xsol2
     dxs=(xm2-xm1)/n1
     !  1st slab
     do l=1,n1
        ip=i1-1+l
        xo(ip)=xm1+l*dxs-dxs/2.d0
     end do
     !  High density slab
     dxs=(xsol2-xsol)/n2
     do l=1,n2
        ip=i1-1+l+n1
        xo(ip)=xsol+l*dxs-dxs/2.d0
     end do
  endif


end subroutine denprof
