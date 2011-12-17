
!  ===========================
!
!   Time-averaged ion density
!
!  ===========================

subroutine rhoiav
  use bopsvars
  implicit none

  integer isnap, i
  real*8 xln

  character chx*15
  real*8 denav(0:nx+1)
  real*8 :: work1(0:nx)
  real*8 xmin, nmin, xmax, nmax, xstep, nstep
  data isnap/0/
  save isnap

  xmin=0.
  xln = xllab/xconv
  xstep = 0.2*xln
  nmin=0.
  nmax = 2.*rho0lab
  nstep = nmax/2.
  if (iunits.eq.2) then
     chx = '     x/mu      '
  else if (iunits.eq.1) then
     chx = '     k_px      '
  else
     chx = '     k_0x      '
  endif

  if (mod(itime,igr).eq.0) then
     do i=1,nx
        work1(i) = denav(i)
     end do
     if (itime.gt.0 .and. ni.gt.0 ) then
        call grxy(xx,work1,nx,2200+isnap,igxs,-1 &
             ,chx,'       <ni>/nc ','niav'//ctime(1:12) &
             ,xmin,xln,xstep,nmin,nmax,nstep)
        isnap=isnap+1

     else
        !  initial profile
        do i=1,nx
           work1(i)=rhoi(i)/gam0**3
        end do
        call grxy(xx,work1,nx,2200+isnap,igxs,-1 &
             ,chx,'       ni/nc   ','nit0'//ctime(1:12) &
             ,xmin,xln,xstep,nmin,nmax,nstep)
        isnap=isnap+1
     endif

     do i=1,nx+1
        denav(i)=0.
     end do

     !  running average of density
  else if (mod(itime,igr).ge.igr-ndenav) then
     !  lab frame ion density = rhoi/gam0**3 if ions fixed.

     do i=1,nx
        if (ioboost.eq.1) then
           !  boost frame
           denav(i)=denav(i) + rhoi(i)/ndenav
        else
           !  lab frame
           denav(i)=denav(i) + (rhoi(i)-vy0*jyi(i))/gam0/ndenav
        endif
     end do

  endif
end subroutine rhoiav










