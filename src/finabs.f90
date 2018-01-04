
!  =======================
!
!   Final absorption data
!
!  =======================

subroutine finabs
  use bopsvars
  implicit none

  real*8 uabs, r, eabs, escabs, etav, erhlab, erinj
  real*8 uhote, fhotpe, uhoti, fhotpi, uemp, utmp, uesp, utotal, tlas
  real*8 uhotp, fhotpp

  integer netav, i

!  if (Uemin.eq.0 .or. ilas.eq.0) return

  if (ioboost.eq.0) then
! Rescale EM energy densities
   Uemin = gam0*Uemin
   Upoy_in=gam0*Upoy_in
   Upoy_out=gam0*Upoy_out
  endif

  !  Lab frame Poynting flux is Upoynt' x cos(thet0) 
  !  = power directed onto target
  
  ! Reflected flux
  Uemout = Uemin-Upoy_in

  ! Transmitted
  Uemtr = Upoy_out

  !  em energy remaining in sim box (lab frame)
!  Uemp = Uema(ntav) + Utma(ntav)
  Uemp = uem(ntc)


  !  Absorbed energy is difference between Poynting vectors - field energy on grid
  Uabs = Upoy_in - Upoy_out - Uemp

  if (uemin.eq.0) uemin=1.

  !  reflectivity
  R=Uemout/Uemin
  !  absorption fraction
  eabs=Uabs/Uemin
  !  fraction of energy into escapee electrons
  escabs=(erhb+elhb)/Uemin

  !  abs measured in sim frame
  etav=0.
  netav=0
  do i=1,ntav
     if (eta(i).gt.0) then
        netav=netav+1
        etav = etav + eta(i)
     endif
  end do
  netav = max(netav,1)
  etav=etav/netav

  !  lab frame escapee energy
  erhlab=0.
  do i=1,nesc
     erhlab=erhlab+uesc(i)*me/gam0**2
  end do

  !  lab frame reinjected electron energy
  erinj=0.
  do i=1,nesc
     erinj=erinj+uinj(i)*me/gam0**2
  end do

  do i=1,nesci
     erhlab=erhlab+uesci(i)*m(i)/gam0**2
  end do

  !  net absorbed particle energy
  erhlab = erhlab-erinj

  !  fraction hot electrons in plasma
!  Uhote=Uthea(ntav)-Uthea(1)
  Uhote = Uthe(ntc) - Uthe(1) 
  fhotpe=Uhote/Uemin

  !  fraction hot ions in plasma
!  Uhoti=Uthia(ntav)-Uthia(1)
  Uhoti = Uthi(ntc) - Uthi(1)
  fhotpi=Uhoti/Uemin

  !  fraction protons in plasma
  Uhotp=Uthpa(ntav)-Uthpa(1)
  fhotpp=Uhotp/Uemin

  !  es energy left in plasma
!  Uesp = Uesa(ntav) - Uesa(1)
  Uesp = Ues(ntc) - Ues(1)

  !  Total abs from components
  Utotal = erhb+elhb+uhote+uhoti+uhotp

  write (15,101) Uemin,Uemout,Uemtr,Uemp,Uabs &
       ,erhb,elhb,Uhote,Uhoti,Uhotp &
       ,Utotal &
       ,R,eabs*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

  !  write to graphics header and to stdout
  write (20,101) Uemin,Uemout,Uemtr,Uemp,Uabs &
       ,erhb,elhb,Uhote,Uhoti,Uhotp &
       ,Utotal &
       ,R,eabs*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

  write (6,101) Uemin,Uemout,Uemtr,Uemp,Uabs &
       ,erhb,elhb,Uhote,Uhoti,Uhotp &
       ,Utotal &
       ,R,eabs*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

101 format (//' Total EM energy in:               ',f10.3 &
       /' EM energy reflected:              ',f10.3 &
       /' EM energy transmitted:            ',f10.3 &
       /' Field energy still in box:        ',f10.3 &
       /'                                        -----' &
       /' EM energy absorbed:               ',f10.3 &
       /'                                        -----' &
       /' Energy out of RHB:                ',f10.3 &
       /' Energy out of LHB:                ',f10.3 &
       /' Plasma electron heating:          ',f10.3 &
       /' Ion heating:                      ',f10.3 &
       /' Proton heating:                   ',f10.3 &
       /'                                        -----' &
       /' Sum of particle energies:         ',f10.3 &
       /'                                        -----' &
       /' Time-integrated reflectivity:     ',f10.4 &
       /' Overall absorption (%):           ',f10.2 &
       /' Energy absorbed by escapee e- (%) ',f10.2 &
       /' Energy in hot plasma electrons (%)',f10.2 &
       /' Energy in hot plasma ions (%)     ',f10.2 &
       /' Energy in protons         (%)     ',f10.2 //)

  tlas=trise+tfall+2*xm1
  if (tlas.gt.nt*dt) then
     call warn('Laser not fully reflected     ')
  endif
end subroutine finabs
