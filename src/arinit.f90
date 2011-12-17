
! ============= ==========
!
!   Array initialization
!
!  $ Revision: $
!
! =======================

subroutine arinit
  use bopsvars
  implicit none
  integer i,ncyc_delay
  real*8 smarg, dxko2, tdelay

  do i=0,nx+1
     if (ioboost.eq.1) then
        !  plots in boost frame
        xx(i) = (i-1)*dx/xconv
     else
        !  lab frame x-axis stretched relative to code frame
        xx(i)=(i-1)*dx*gam0/xconv
     endif
     xk(i)=(i)*dkx
     ff(i)=0.
     fb(i)=0.
     jm(i)=0.
     jp(i)=0.
     jye(i)=0.
     jyem(i)=0.
     jyi(i)=0.
     jyim(i)=0.
     gf(i)=0.
     gb(i)=0.
     jzm(i)=0.
     jzp(i)=0.
     jze(i)=0.
     jzem(i)=0.
     jzi(i)=0.
     jzim(i)=0.
     ex(i)=0.
     ey(i)=0.
     ez(i)=0.
     bx(i)=0.
     by(i)=0.
     bz(i)=0.
     epond(i)=0.
  end do

  do i=1,nf
     xw(i)=i*dw
  end do

  do i=1,nptrk
     uest(i)=0.
     phat(i)=0.
  end do


  !   smoothing and boost factor for fft
  do i=1,nxo2
     dxko2=i*pi/nx
     smarg=max(abo*sin(dxko2)**2-asm*tan(dxko2)**4,-30.d0)
     sm(i)=exp(smarg)
     rk2(i)=1./((2.0*sin(dxko2)/dx)**2)*sm(i)**2
  end do

  tdelay=2*xm1
  ncyc_delay=tdelay/tav+0.5
  write (6,*) '# delay cycles =',ncyc_delay

  do i=0,ncyc_delay
     xi0(-i)=0.
     Uinc(-i)=0.
  end do

end subroutine arinit
