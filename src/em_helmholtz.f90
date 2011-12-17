
!     ==================================
!     
!     Helmholtz solver for electromagnetic fields 
!     
!     
!     ==================================

!     call em_helmholtz(itime,itav,nx,dx,dt,theta0,ilas,a0,w0,trise,tpulse,dcrhe,Az,Ez,By,Epond)

subroutine em_helmholtz(itime,itav,n,dx,dt,theta,ilas,a0,w0,trise,tpulse,rhoe,Azr,Ezr,Byr,Eyr,Bzr,epond)

  implicit none
  integer, intent(in) :: n, itime, ilas, itav
  real*8, intent(in) :: theta, a0, w0, dx, dt, trise, tpulse
  real*8, intent(in) :: rhoe(0:n+1)  ! cycle-averaged electron density
  real*8, dimension(0:n+1), intent(out) :: Ezr, Byr, Eyr, Bzr, Azr, epond ! fields including time-dep phase factor
  complex, dimension(n) :: alpha,beta,gamma,y
  complex, dimension(0:n+1) :: Az,Ao,Ay	! Vector potential
  complex, dimension(n) :: Ey, Ez, By, Bx, Bz ! Fields derived from Az
  real, dimension(n) :: eps  ! permittivity
  ! real :: rho(0:n+1)  ! density
  real :: rgam(0:n+1)  !  relativistic gamma
  real :: err(0:n+1)
  complex :: yi,  aave, carg, cphase
  integer :: i,j,n_iter
  real*8 :: pi, s2th, cth, errmax
  real*4 :: pha
  real*8 :: t                    ! time
  real*8 :: tptot               ! pulse duration
  real*8 :: lskin

  real*8 :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real*8 :: amplitude, Tpon, gamma_pon, phi_m, g0, nu_eff
  integer :: iplas
  !     linear rise

  t=dt*itime

  pi = asin(1.0)*2
  yi = (0.,1.)
  nu_eff = 0.2
  n_iter = 3

! Pulse envelope
  if (mod(ilas,10).eq.1) then
! linear rise-time
     Tpon = min(1.d0,t/trise)

  else if (mod(ilas,10).eq.5) then
     !     convert from fwhm to total length
     tptot=tpulse*2

     if (t.le.tptot) then
	phi_m = 1.d0*pi*t/tptot
        Tpon=sin(phi_m)
     else
        Tpon=0.
     endif 

  endif

! Laser phase (C or P-pol)
  if (ilas.le.10) then
   pha=w0*t
  else
    pha=0.  ! remove quick phase: emulate C-pol light
  endif

! write (*,*) ilas,itime,pha
  amplitude=a0*tpon
  cphase = cexp(yi*pha)

! write(*,*) 'a0, t, tptot, tpon',a0, t, tptot, phi_m, tpon, amplitude


  err=0.
  !  s2th=sin(pi/180*theta)**2
  s2th=0.
  cth = cos(pi/180*theta)

  Az=(0.,0.)
  ao=(0.,0.)
  rgam(0:n+1)=1.

  do j=1,n_iter


     do i=1,n
        !  coefficients
        ! rhoe normalized to nc
	! include absorption fraction via effective collision freq.
        eps(i) = 1.-rhoe(i)/rgam(i)/(1+yi*nu_eff)
        y(i)=(0.,0.)
        alpha(i)=1
        beta(i)=-2 + dx**2*(eps(i)-s2th)
        gamma(i)=1
     end do

     !  BCs
     y(1) = 2*yi*amplitude*sin(dx*cth)
     carg = yi*dx*cth
     beta(1) = beta(1) + cexp(carg)

     call trisolve(alpha,beta,gamma,y,Az(1:n),n-1,n)


     ! relativistic factor - average old and new potentials
     errmax=0.
     do i=1,n
        rgam(i) = sqrt(1 + 0.5*abs(az(i))**2) 
        err(i) = sqrt(abs(az(i)**2-ao(i)**2)/4/amplitude**2)
     end do

     ! BCs for next iterate (Laplacian of gamma)
     rgam(0) = 2*rgam(1) - rgam(2)
     rgam(n+1) = rgam(n)

     errmax = maxval(err(1:n))
     Ao = Az  ! Store previous iterate

  end do
iplas = 9./dx
!write(*,'(i6,2f12.3)') itime,rhoe(iplas),abs(az(iplas))
if (itime .eq. itav) then
  write (*,'(a20,i2,a10,f12.5)') 'Iterate ',j,' error=',errmax
  g0 = sqrt(1+amplitude**2.2)
  open (40,file='a_error.dat')
  write(40,'(6(1pe12.3))') (dx*i,rhoe(i),eps(i),abs(az(i))/amplitude,rgam(i)/g0,err(i),i=1,n)
  close(40)
  
endif

  ! Derive fields from Az
  ! Bcs
  Az(0) = 2*Az(1) - Az(2)
  Az(n+1) = Az(n)

  Ay = yi*Az

  do i=1,n
     Ey(i) = yi*Ay(i)
     Ez(i) = yi*Az(i)
     By(i) = -(Az(i+1)-Az(i-1))/2/dx  ! By=iEz'
     Bz(i) = (Ay(i+1)-Ay(i-1))/2/dx  ! Bz=-iEy'
!     Bx(i) = sin(theta)*Ez(i) 
  end do

  Eyr(1:n) = Real(Ey(1:n)*cphase)*tpon  ! reconstruct real EM fields including time-variation (s-pol)
  Ezr(1:n) = Real(Ez(1:n)*cphase)*tpon  ! reconstruct real EM fields including time-variation (s-pol)
  Byr(1:n) = Real(By(1:n)*cphase)*tpon
  Bzr(1:n) = Real(Bz(1:n)*cphase)*tpon
  Azr(1:n) = Real(Az(1:n)*cphase)*tpon
!  Ezr(1:n) = Real(Ez(1:n))*tpon  ! reconstruct real EM fields excluding time-variation (s-pol)
!  Byr(1:n) = Real(By(1:n))*tpon
!  Azr(1:n) = Real(Az(1:n))*tpon
  ! pond force - without gamma factor, as in emfield
  do i=1,n
     epond(i) = .25/dx*( azr(i+1)**2-azr(i-1)**2 )
  end do

end subroutine em_helmholtz









