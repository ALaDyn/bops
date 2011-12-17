!      *********************************************
!      *                                           *
!      * 1D3V EM PIC CODE FOR FS-SOLID SIMULATION  *
!      *                                           *
!      *********************************************
!
!      Author:        Paul Gibbon
!		      Institute for Advanced Simulation
!		      Juelich Supercomputing Centre
! 		      Research Centre Juelich
!




program bops
  use bopsvars
  implicit none
  real tleft, tst, tcyc, tnodia, tdiag, tim1, tr, td1, td2, tim2
  real tpush, tdens, tfield, tcur, ts1, ts2
  integer i1, iflag, pdi

  !      open O/P files
  call fileop
  call stamp(15,1)
  call stamp(6,1)
  call init
  !      input data
  call derdat
  !  allocate array space
  call setup_arrays
  call arinit
  !  write out header page
  !  call header

  !  load particles

  if (lrstrt) then
     !      read particle variables and fields from stream 9
     call restrt
     !      snapshot counter
     idc=(itim0-1)/igr
     call eden
  else
     itim0=1
     if (itropt.eq.3.or.itropt.eq.2) call trkini
     call parload
  endif

  call esfield
  call laser

!enable counterpropagating laser pulse from the right boundary (Qiao 2009 additon!)

!  call laser_c

!#############################

  call pouti
  call diagnostics
  call pout
  !      ibm timer
  itime=itim0
  !      call tremain(tleft)
  tleft=1000.
  tst=tleft

  !      MAIN LOOP

  tdiag=0.
  tpush =0.
  tdens = 0.
  tcur = 0.
  tfield = 0.

  call cputime(tim1,i1)

 do itime=1,nt

  tr=itime*dt
!      push particles and gather density
!      start pushing when laser reached surface
  if (ifreeze.eq.0 .or. tr.ge.xm1-5*dx) then

  call cputime(ts1,i1)

!## Bin & Anupam 09/10 
!## choose 3Dv push 'pushc' for all cases. Particle pusher selection abandoned.

!   pusher: select case(push_scheme)  ! Select particle pusher according to EM scheme
! 
!   case(1)  ! Full EM, p-pol
!      call pushp(1,ne,dte)
!   case(2) !      s-pol
!      call pushs(1,ne,dte)
!   case(3) !      c-pol
!      call pushc(1,ne,dte)
!   case(4)     !      fpond ES pusher
!      call push_pond(1,ne,dte)
!   case(5)  !      ES pusher, fpond computed in emfield
!      call push_pondes(1,ne,dte)
!   case default
!      call pushp(1,ne,dte)
!   end select pusher

!## 3D Particle pusher scheme by default!

  call pushc(1,ne,dte)

!##########################################

  call move(1,ne,dte)
  call pbcs(1,ne,vte,me,ipbc)

  call cputime(ts2,i1)
  tpush = tpush + ts2-ts1

  call eden

  call cputime(ts1,i1)
  tdens = tdens + ts1-ts2



  if ( mod(itime,isubi).eq.0 ) then

!#### Anupam & Bin 09/10: ion pusher scheme selection abandoned.
!#### choose 3D v push scheme "pushc" for all cases!

! Select particle pusher according to EM scheme
! Note: all ions pushed together

!   pusheri: select case(push_scheme)  
!      case(1)  ! Full EM, p-pol
!         call pushp(ion1,ni_tot,dti)  
!      case(2)
!         call pushs(ion1,ni_tot,dti)  !      s-pol
!      case(3)
!         call pushc(ion1,ni_tot,dti)   !      c-pol
!      case(4) !      quasi-ES pusher
!         call push_pond(ion1,ni_tot,dti)
!      case(5)  !      ES pusher
!         call push_pondes(ion1,ni_tot,dti)
!      case default
!         call pushp(ion1,ni_tot,dti)  
!      end select pusheri

     call pushc(ion1,ni_tot,dti)

!#### Anupam & Bin 09/10

     call move(ion1,ni_tot,dti)

!#### Anupam & Bin 09/10: separately deal with protons & heavy ions at the boundary

     call pbcs(ion1,np,vtp,mp,ipbc)
     call pbcs(ion1,ni_tot-np,vti,mi,ipbc)
!#### Anupam & Bin 2009/2010
     
     call cputime(ts2,i1)
     tpush = tpush + ts2-ts1
     call iden(ion1)
     call cputime(ts1,i1)
     tdens = tdens+ ts1-ts2
  endif

  endif

  !      currents
  call jeden
  call jiden(ion1)
  !      add currents -> jtot for field solver
  call addcurr
  call cputime(ts2,i1)

  tcur = tcur + ts2-ts1

  !      electrostatic field
  call esfield

  call rhoeav  ! Accumulate DC electron density (needs currents for lab frame)

  !      electromagnetic field

  em_mode: select case(em_scheme)  ! Select EM propag scheme
  case(2)
     call empond ! Helmholtz step solutions for Ez, By, Az (no epond)
  case(3)
     call fpond  ! Helmholtz step solution for Ez, By, Az, epond
  case(5)
     call em_helmholtz(itime,itav,nx,dx,dt,theta0,ilas,a0,w0,trise,tpulse,dcrhe,Az,Ez,By,Ey,Bz,Epond)
  case default
     call emfield
  end select em_mode  ! Select EM propag scheme

  call cputime(ts1,i1)
  tfield = tfield + ts1-ts2

  !      position shift new -> old
  call pshif
  call cputime(ts2,i1)
  tpush = tpush + ts2-ts1

  !      waveguide mode
  !      if (iembc.eq.3) call wguide

  call cputime(td1,i1)
  call diagnostics
  call pout
  call cputime(td2,i1)
  tdiag=tdiag+td2-td1
  call cputime(tim2,iflag)
  ttot=tim2-tim1

  !      estimated run time
  if (itime.eq.1) then
     call warn('Estimated run time            ')
     call ruform('Trun: ',ttot*nt/60.,'f12.3','min   ')
     call ruform('Trun: ',ttot*nt/3600.,'f12.3','hours ')
     call d0prn6('Trun m',ttot*nt/60.)
     call d0prn6('Trun h',ttot*nt/3600.)
     call blank
     call blank6
  endif

  end do


  !      Timing
  call cputime(tim2,iflag)
  call warn('loop time (s)                  ')
  ttot=tim2-tim1
  tdiag=tdiag
  tcyc=ttot/nt
  tnodia=(ttot-tdiag)
  call r0form('total ',ttot*1.d0,'f15.3')
  call r0form('push  ',tpush*1.d0,'f15.3')
  call r0form('dens  ',tdens*1.d0,'f15.3')
  call r0form('curr  ',tcur*1.d0,'f15.3')
  call r0form('field ',tfield*1.d0,'f15.3')
  call r0form('cycle ',tcyc*1.d0,'f15.3')
  call r0form('nodiag',tnodia*1.d0,'f15.3')
  call r0form6('Trun/s',ttot,'f15.2')
  call r0form6(' push ',tpush*1.d0,'f15.3')
  call r0form6(' dens ',tdens*1.d0,'f15.3')
  call r0form6(' curr ',tcur*1.d0,'f15.3')
  call r0form6(' field',tfield*1.d0,'f15.3')
  call r0form6(' diag',tdiag*1.d0,'f15.3')
  call r0form6('Tr-Td ',tnodia*1.d0,'f15.3')
  call r0form6('Tcycle',tnodia/nt*1.d0,'f16.4')
  write (15,'(a,f12.3,a)') &
       'Tcpu/particle/dt',1.e6*tnodia/nt/ne,' mu s'
  write (6,'(a,f12.3,a)') &
       'Tcpu/particle/dt',1.e6*tnodia/nt/ne,' mu s'

  if (itime.le.nt) then
     igr=itime
     call warn('Ran out of time!              ')
     call i0prnt('# dt  ',itime)
     call gsnap
  endif
  if (itropt.eq.2) then
     call plotrk
  else if (itropt.eq.3) then
     call scatrk
  endif
  if (nt.ge.itc) call ghist

  !      Time-integrated absorption data summary
  call finabs

  if(ldump) call savear

  !      Tidy up output files
  call fileclo
  call stamp(15,2)
  call stamp(6,2)

end program bops
