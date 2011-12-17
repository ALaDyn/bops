
!  ====================
!
!  Run-data header page
!
!  ====================

subroutine header
  use bopsvars
  implicit none

  integer npars
  real*8 xnclab, vboost, qelab, cq, xline, tspec
  npars=26
  !     write (20,101) 0,npars
101 format(2i6)
  call stamp(20,1)
  call hop('#elecs',1.d0*ne,0)
  call hop('# ions',1.d0*ni,0)
  call hop('vosc/c',a0*1.,3)
  call hop('I_{18}',a0*a0*1.40,5)
  !      call hop('lambda',xlambda*1.,3)
  !      call hop(' sigma',sigma*1.,1)
  !      call hop('wp/w0 ',1.*wp/gamma0**1.5,2)
  xnclab=wp*wp/gam0**3
  call hop('rho0/nc ',1.*xnclab,1)
  call hop('L/lamd',1.*xlolam,4)
  call hop('theta0',1.*theta0,1)
  call hop('Pol: '//cpolzn,1.*ppol+spol,2)
  vboost=vy0
  !      call hop('   vy0',1.*vboost,3)
  call hop('  gam0',1.*gamma0,2)
  call hop('    Te',vte*vte*511,3)
  call hop('    Ti',vti*vti*511*miome,3)
  call hop(' mi/me',1.*miome,0)
  !      call hop('Pos/Pe',a0**2/2./xnclab/vte**2,-2)
  !      call hop('     Z',1.*z,0)
  !      call hop('niprof',inprof*1.,0)
  call hop('    nt',nt*1.d0,0)
  !      call hop(' itim0',itim0*1.,0)
  call hop('  Trun',1.*nt*dt/tconv,1)
  call hop('   *w0',1.*nt*dt,1)
  call hop('   /fs',1.*nt*dt/tcfs,1)
  call hop('Tpulse',1.*tpulse/tconv,1)
  call hop('   /fs',1.*tpulse/tcfs,1)
  !      call hop(' Trise',1.*trise/tconv,1)
  !      call hop(' Tfall',1.*tfall/tconv,1)
  call hop('    nx',nx*1.d0,0)
  call hop('sim dx',1.*dx/xconv,4)
  call hop('    dt',1.*dt/tconv,4)
  !      call hop('Ldebye',1.*vte/wp,4)
  call hop('Ld/dx ',1.*vte/wp/dx,2)
  call hop('lab xl',1.*xllab/xconv,1)
  call hop('vacuum',1.*xm1*gamma0/xconv,2)
  call hop(' solid',1.*xsol*gamma0/xconv,2)
  call hop(' dfoil',1.*dfoil/xconv,2)
  call hop('  d/mu',1.*dfoil/xmu,2)
  !  line density
  !   qe(lab)= qe(sim)/gam0**2
  qelab=qe/gam0**2
  !  particles per cell
  ncell = -nonc*dx/qe
  !      call hop('   qe ',1.*qelab,-3)
  call hop(' ncell',ncell*1.,0)
  cq = 5.e8/xlambda*sigma**2*abs(qelab)
  !      call hop('Ne/Nsm',1.*cq,-5)
  xline = -xnclab/qe*gam0**2
  call hop(' ldens',1.*xline,-3)
  call hop(' Tsnap',1.*igr*dt/tconv,1)
  call hop('Tdelay',2.*xm1/tav,1)
  tspec = dt*nf*ift
  call hop(' Tspec',1.*tspec/tconv,1)
  call hop('  wmax',1.*pi/dt/ift,1)
  call hop(' dwbin',2.*pi/tspec*ifbin,4)
  call hop(' units',iunits*1.d0,0)
end subroutine header



