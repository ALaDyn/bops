
!  ================================
!
!   Echo input data to output file
!     and check array space
!
!  ================================

subroutine pouti
  use bopsvars
  implicit none

  integer nspec, io, i
  real*8 vboost, dterr, tsp0, tsp1, dtsp,cq, qelab

  vboost=vy0
  qelab=qe/gam0**2

  do io=6,15,9
     write (io,'(/a/a)') 'Laser parameters' &
          ,                      '----------------'
     write (io,'(a30,1pe13.4,a15)') 'Intensity' &
          ,a0**2*1.38e18/xlambda**2,'W cm^-2'
     write (io,'(a30,1pe13.4,a15)') 'I-lambda^2' &
          ,a0**2*1.38e18,'W cm^-2 mu^2'
     write (io,'(a30,f13.4)') 'vosc/c',a0
     write (io,'(a30,f13.3,a15)') 'Wavelength',xlambda,'microns'
     write (io,'(a30,1pe13.3,a15)') 'Frequency',1.88/xlambda*1.e15 &
          ,'/s'
     write (io,'(a30,f13.3,a15)') 'Pulse duration',tpulse/tcfs &
          ,'femtoseconds'
     write (io,'(a30,f13.3,a15)') 'Pulse delay',tdel/tcfs &
          ,'femtoseconds'
     write (io,'(a30,f13.2,a15)') 'Incidence angle',theta0 &
          ,'degrees'
     write (io,'(a30,a13)') 'Polarization',cpolzn
     write (io,'(a30,f13.4,a15)') 'Spot size',sigma,'microns'
     write (io,'(a30,1pe13.3,1pe13.3,a8)') 'Chirp',bchirp &
          ,bchirp*w0r**2,'/ps^2'
     write (io,'(/a/a)') 'Target parameters' &
          ,                      '-----------------'
     write (io,'(a30,f13.4)') 'Density scale-length L/lambda ' &
          ,xlolam
     write (io,'(a30,f13.4)') 'k0 L',xlolam*2*pi
     write (io,'(a30,f13.4)') 'Max density n/n_c ',rho0lab
     write (io,'(a30,2f13.4)') 'Grid length (norm/microns) ' &
          ,xllab,xllab/xmu
     write (io,'(a30,2f13.4)') 'LH vacuum (norm/microns) ' &
          ,xm1*gamma0,xm1*gamma0/xmu
     write (io,'(a30,2f13.4)') 'RH vacuum (norm/microns) ' &
          ,xm2*gamma0,xm2*gamma0/xmu
     write (io,'(a30,2f13.4)') 'Shelf position (norm/microns) ' &
          ,xsol*gamma0,xsol*gamma0/xmu
     write (io,'(a30,2f13.4)') 'Foil width (norm/microns) ' &
          ,dfoil*gamma0,dfoil*gamma0/xmu
     write (io,'(a30,2f13.4)') 'Destruction threshold  n/nc d/lambda, a0 ' &
          ,rho0lab*dfoil*gamma0/xlambda,a0
     if (inprof.eq.5) then
        denmin = rho0*exp(-(xsol-xm1)/xlolam)
        write (io,'(a30,f13.4)') 'Exp. profile nmin ',denmin
     else
        write(io,'(a30,i8)') 'Density profile:',inprof
     endif
     if (target_config.gt.1) then
	write (io,'(/a/a)') 'Additional layer(s)' &
          ,                      '-----------------'
     	write (io,'(a30,i8)') 'Config:',target_config
     	write (io,'(a30,2f13.4)') 'Layer width (norm/microns) ' &
          ,x_layer*gamma0,x_layer*gamma0/xmu
     	write (io,'(a30,f13.4)') 'Layer density n_p/n_c ',rho_layer/gam0**3

     endif
     write (io,'(/a/a)') 'Plasma characteristics' &
          ,                      '-----------------'
     write (io,'(a30,f13.4)') 'Boost velocity ',vboost
     write (io,'(a30,f13.4)') 'Boost gamma ',gamma0
     write (io,'(a30,f13.4)') 'Plasma frequency wp/w0 ',wp
     write (io,'(a30,f13.4)') 'Electron temperature (keV) ',Te
     write (io,'(a30,f13.4)') 'Thermal velocity (vte/c) ',vte
     write (io,'(a30,2f13.4)') 'Debye length (k0 l_D; l_D/dx)' &
          ,xdodx*dx,xdodx
     write (io,'(a30,f13.4)') 'Ion temperature (keV) ',Ti
     write (io,'(a30,f13.4)') 'Ginzburg angle ',theta_opt
     cq = 5.e8/xlambda*sigma**2*abs(qelab)
     write (io,'(a30,1pe13.4)') 'Charge conversion factor N/Nsim',cq
     write (io,'(a30,2f13.4)') 'Charge on grid/particles ' &
          ,Qgrid/gamma0**2,ne*qe/gamma0**2

     write (io,'(/a/a)') 'Control parameters' &
          ,                      '------------------'
     write (io,'(a30,i13)') '# electrons ',ne
     write (io,'(a30,i13)') '# ions ',ni
     write (io,'(a30,i13)') '# protons ',np
     write (io,'(a30,i13)') '# grid points ',nx
     write (io,'(a30,f13.2)') '# grid points across foil ',(xm2-xm1)/dx
     write (io,'(a30,f13.2)') '# particles/cell (foil) ',ne/((xm2-xm1)/dx)
     write (io,'(a30,i13)') '# timesteps ',nt
     write (io,'(a30,i13)') 'Pusher scheme: ',push_scheme
     write (io,'(a30,i13)') 'history frequency itc ',itc
!     write (io,'(a30,i13)') '# steps in time-ave plots',ntc
    write (io,'(a30,2f13.4)') 'Run time (norm/fs)',nt*dt,nt*dt/tcfs
     write (io,'(a30,2f13.4)') 'Mesh size (norm/microns)', &
          dx,dx/xconv
     write (io,'(a30,2f13.4)') 'Time step (norm,fs)',dt,dt/tcfs
     write (io,'(/a/)') '-----------------------------------'
  end do


  call i0prnt('itime ',itime)
  call i0prnt('itim0 ',itim0)


  call i0prnt(' isube',isube)
  call i0prnt(' isubi',isubi)
  call i0prnt('itropt',itropt)
  call i0prnt('  ipbc',ipbc)
  call i0prnt('  ifbc',ifbc)
  call i0prnt('  ilas',ilas)
  call i0prnt('   igr',igr)
  call i0prnt('  idia',idia)
  call i0prnt('  iout',iout)
  call i0prnt('   itc',itc)
  call i0prnt('  igxs',igxs)
  call i0prnt(' units',iunits)
  call i0prnt('igmovi',igmovie)
  call i0prnt('  itav',itav)
  call r0form('   tav',tav,'f12.2')
  call i0prnt('ncycav',ncyc)
  call r0form('    dt',dt,'f14.4')
  call r0form('    dx',dx,'f14.4')

  call r0form(' tconv',tconv,'f12.2')
  call r0form(' xconv',xconv,'f12.2')
  call r0form('     z',z,'f11.1')
  call r0form('    w0',w0,'f12.2')
  call r0form('   vti',vti,'f13.3')

  call r0form(' xsol2',xsol2*gamma0/xconv,'f12.2')
  call r0form(' xcrit',xcrit*gamma0/xconv,'f12.2')
  call r0form('l wall',wl*gamma0/xconv,'f12.2')
  call r0form('r wall',wr*gamma0/xconv,'f12.2')

  call i0prnt('ntrack',ntrack)
  call i0prnt('istart',itstart)
  call i0prnt('iend  ',itend)
  call i0prnt('   nf ',nf)
  call i0prnt('  ift ',ift)
  call r0form(' miome',miome,'f11.1')
  call r0form('  qome',qome,'f12.2')
  call r0form('  qomi',qomi,'f12.2')
  call r0form('    qe',qe/gamma0**2,'f16.6')
  call r0form('    qi',qi/gamma0**2,'f16.6')
  call r0form('    qp',qi/gamma0**2,'f16.6')
  !      call r0form('  xgin',xgin,'f13.3')
  call i0prnt('  isgx',igx2d)
  call i0prnt('  isgy',isgy)
  call r0form('  uhot',uhot,'f12.2')
  dterr=abs(dt-dx)/dx

  if (dterr.gt.0.001.and.ilas.gt.0) then
     call warn(' timestep not equal to dx/c    ')
     stop
  endif

  if (xdodx.lt.0.1) then
     write(6,*)
     write(6,*) 'Debye length badly under-resolved!'
     write (*,'(a,f12.3)') ' LD/dx =',xdodx
     write(6,*) ' - run will suffer numerical heating'
     write(6,*)
  else if (xdodx.lt.0.5) then
     call warn('Debye length not resolved     ')
     write (*,'(a,f12.3)') 'LD/dx=',xdodx
  endif
  if (xsol.gt.xl) then
     call warn(' Not enough space for profile  ')
     call warn(' Aborting run                  ')
     stop
  endif

  tsp0 = 2*xm1
  !      tsp1 = 2*xm1+nf*ift*dt
  !      tsp0=0.
  tsp1 = nf*ift*dt
  dtsp = nf*ift*dt
  nspec = (nt-tsp0/dt)/nf/ift
  do io=6,15,9
     write (io,*) 'Spectrum sampling start:',tsp0/tconv,' itime=', &
          tsp0/dt/tconv
     do i=0,nspec-1
        write (io,*) 'Outputs:',(tsp1+i*dtsp)/tconv,' itime=', &
             (tsp1+i*dtsp)/dt
     end do
  end do
end subroutine pouti





