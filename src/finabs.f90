
!  =======================
!
!   Final absorption data
!
!  =======================

subroutine finabs
  use bopsvars
  implicit none

  real*8 uabs, r, eabs, escabs, etav, erhlab, erinj
  real*8 uhote, fhotpe, uhoti, fhotpi, uemp, uesp, utotal, tlas
  real*8 uhotp, fhotpp
  integer netav, i

  if (Uemin.eq.0 .or. ilas.eq.0) return

  !  Lab frame i/p energy is Uem x cos(thet0) = power
  !  directed onto target

  Uemin=Uemin/gam0
  Uemout=Uemout/gam0
  Uemtr = Uemtr/gam0

  !  Absorbed energy
  Uabs = Uemin - (Uemout+Uemtr)
  if (uemin.eq.0) uemin=1.

  !  reflectivity
  R=Uemout/Uemin
  !  average absorption
  eabs=Uabs/Uemin
  !  fraction of energy into escapee electrons
  escabs=erhb/Uemin

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

  !  net absorped particle energy
  erhlab = erhlab-erinj

  !  fraction hot electrons in plasma
  Uhote=Uthea(ntav)-Uthea(1)
  fhotpe=Uhote/Uemin

  !  fraction hot ions in plasma
  Uhoti=Uthia(ntav)-Uthia(1)
  fhotpi=Uhoti/Uemin

  !  fraction protons in plasma
  Uhotp=Uthpa(ntav)-Uthpa(1)
  fhotpp=Uhotp/Uemin

  !  em energy remaining
  Uemp = Uema(ntav) - Uema(1)

  !  es energy left in plasma
  Uesp = Uesa(ntav) - Uesa(1)
  !  Total abs from components
  Utotal = erhb+uhote+uhoti+uhotp+uemp+uesp

  write (15,101) Uemin,Uemout,Uemtr,Uabs &
       ,erhb,erhlab,Uhote,Uhoti,Uhotp,Uemp,Uesp &
       ,Utotal &
       ,R,eabs*100.,etav*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

  !  write to graphics header and to stdout
  write (20,101) Uemin,Uemout,Uemtr,Uabs &
       ,erhb,erhlab,Uhote,Uhoti,Uhotp,Uemp,Uesp &
       ,Utotal &
       ,R,eabs*100.,etav*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

  write (6,101) Uemin,Uemout,Uemtr,Uabs &
       ,erhb,erhlab,Uhote,Uhoti,Uhotp,Uemp,Uesp &
       ,Utotal &
       ,R,eabs*100.,etav*100. &
       ,escabs*100.,fhotpe*100.,fhotpi*100.,fhotpp*100.

101 format (//' Total EM energy in:               ',f10.3 &
       /' EM energy reflected:              ',f10.3 &
       /' EM energy transmitted:            ',f10.3 &
       /' EM energy absorbed:               ',f10.3 &
       /' Energy out of rhb: (erhb)         ',f10.3 &
       /' Energy out of rhb: (Uesc)         ',f10.3 &
       /' Plasma electron heating:          ',f10.3 &
       /' Ion heating:                      ',f10.3 &
       /' Proton heating:                   ',f10.3 &
       /' Residual em energy:               ',f10.3 &
       /' Residual es energy:               ',f10.3 &
       /'                                         -----' &
       /' Sum of components:                ',f10.3 &
       /' Time-integrated reflectivity:     ',f10.2 &
       /' Overall absorption (%):           ',f10.2 &
       /' Sim. frame absorption <eta>(%):   ',f10.2 &
       /' Energy absorbed by escapee e- (%) ',f10.2 &
       /' Energy in hot plasma electrons (%)',f10.3 &
       /' Energy in hot plasma ions (%)     ',f10.3 &
       /' Energy in protons         (%)     ',f10.3 //)

  tlas=trise+tfall+2*xm1
  if (tlas.gt.nt*dt) then
     call warn('Laser not fully reflected     ')
  endif
end subroutine finabs
