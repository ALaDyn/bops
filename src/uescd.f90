!     ===========================================
!
!     Energy distribution of escapee particles
!     - electrons
!
!     $ Revision: $
!     ===========================================

subroutine uescd
  use bopsvars
  implicit none

  real*8 beta, ukm, du, dukev, ua, fescm, fluxcav, fmin, uinlab
  real*8 hotke, uht
  real*8 dline, pflux, fluxc, ekth, pflu1, sumt, u, sumf, hotflux
  integer i, l, i0, i3, is, iu1, iu2
  character*15 cu

  real*8, dimension(0:nvx+1) :: work0,work1,work2, work3, work4, work5,work6
  real*8 :: melab,qelab

  if (mod(itime,igr).ne.0 .or. itime.eq.0 .or. nesc.eq.0) return

  melab = me/gam0**2
  qelab = qe/gam0**2
  Te = 511*vte**2
  beta = sqrt(2/pi)

  !     energy distribution of escaped particles - in LAB FRAME
  do i=1,nvx
     work2(i)=0.
     work1(i)=0.
     work4(i)=0.
     work3(i)=0.
     work0(i)=0.
     work5(i)=0.
  end do

  ukm=0.
  do l=1,nesc
     ukm=max(ukm,uesc(l))
  end do


  du = umevmax/.511/nvx

  if (ukm.gt.nvx*du) then
     call warn(' hot electron E >  umevmax: increase umevmax ')
  endif
  dukeV=511*du

  do l=1,nesc
     ua=uesc(l)/du
     iu1 = ua+1
     i0=min(max(1,iu1),nvx)
     work1(i0)=work1(i0)+1
     !     reinjected distn
     ua=uinj(l)/du
     iu2 = ua+1
     i0=min(max(1,iu2),nvx)
     work4(i0)=work4(i0)+1
  end do

  !     distribution functions are:
  !
  !     f_tot(U) dU = work1
  !     f_in(U) dU  = work4
  !     f_hot(U) dU = work1 - work4

  !     get max
  fescm=0.
  do i=1,nvx
     fescm = max(fescm,work1(i))
  end do
  !     normalising flux of cold component, assuming it's Maxwellian
  fluxc = fescm*exp(1.)*Te/dukeV

  !     averaged estimate - integrate over cold part from 0 to 3 x max(Te,dukeV)
  i3 = 3*max(Te,dukeV)/511/du
  fluxcav=0.
  do i=1,i3
     fluxcav = fluxcav + 1.249*work1(i)
  end do

  write (15,*) 'fescm, Te, fluxc, fluxcav, du' &
       ,fescm,Te,fluxc,fluxcav,du

  fmin=.1
  !     lab frame line density
  dline = -rho0lab/qelab
  !     time integrated electron flux out of RHB
  pflux = .5*beta*dline*vte*itime*dt
  !     KE of thermals
  Ekth = 0.
  pflu1 = 0.
  sumt = 0.

  do i=1,nvx
     u=i*dukeV-dukeV/2
     !     flux distribution of thermals
     work5(i) = fluxcav/Te**2*u*exp(-u/Te)*dukeV
     Ekth = Ekth + u**2/Te**2*exp(-u/Te)*dukeV
     !     subtract off reinjected thermals
     sumt = sumt + work1(i)
     work0(i)=work1(i)-work5(i)
     if (work0(i).le.0) work0(i)=fmin
     if (work1(i).le.0) work1(i)=fmin
     if (work4(i).le.0) work4(i)=fmin
     pflu1 = pflu1 + work4(i)
     if (umevmax.ge.1) then
        xv(i)=(i*du-du/2)
        cu = '      U(MeV)   '
     else
        xv(i)=(i*du-du/2)*1.e3
        cu = '      U(keV)   '

     endif

  end do

  !     normalise to nesc
  sumf=0.

  do i=1,nvx
     sumf=sumf+work0(i)
  end do
  hotflux = nesc - fluxcav


  uht=0.
  do i=1,nvx
     u=du*i-du/2
     work2(i)=work0(i)*u*melab
     uht=uht+work2(i)
  end do

  work3(nvx+1)=0.
  do i=nvx,1,-1
     work3(i) = work3(i+1)+work2(i)
  end do

  !     - note scaling of Uemlab and me
  uinlab=uemin/gam0
  hotke = 511*uht/melab/sumf

  !     thermal energy distn
  !     of particles escaped from rh boundary
  !     plotted against k.e.  gamma-1  normalised to mc**2
  !

  call blank
  do is = 15,20,5
     write (is,*)
     write (is,*) ' Hot electron data:'
     write (is,*) '-------------------'
     write (is,1001) Te,Ekth,nesc,pflux,sumt,fluxc,fluxcav,uinlab &
          ,uht,erhb,hotke,sumf,hotflux,du
1001 format(//' Plasma Temp Te          ',f10.3 &
          /' Plasma K.E.             ',f10.3 &
          /' # escaped particles     ',i8 &
          /' Plasma flux  RHB        ',f10.1 &
          /' Total flux  Int Ues fes ',f10.1 &
          /' Flux of cold component  ',f10.1 &
          /'    ditto  averaged      ',f10.1 &
          /' EM energy in            ',f10.3 &
          /' Absorbed energy Uf(U)   ',f10.3 &
          /' Absorbed energy (erhb)  ',f10.3 &
          /' Hot electron K.E.       ',f10.3 &
          /' Hot electron flux       ',f10.1 &
          /' Phi_hot = nesc-Phi_cold ',f10.1 &
          /' du/keV                  ',f10.3//)
  end do

  !     TODO: Change MeV to keV if umevmax<1

  call grxy(xv,work1,nvx,6700+idc,1,2 &
       ,cu,'     fe(Uesc)  ','fues'//ctime(1:12) )
  call grxy(xv,work4,nvx,16800+idc,1,2 &
       ,cu,'     fe(Uinj)  ','fuin'//ctime(1:12) )
  call grxy(xv,work5,nvx,16900+idc,1,2 &
       ,cu,'     ge(Uin)   ','gein'//ctime(1:12) )
  call grxy(xv,work0,nvx,6800+idc,1,2 &
       ,cu,'     fh(U)     ','fhot'//ctime(1:12) )
  call grxy(xv,work2,nvx,7100+idc,1,1 &
       ,cu,'    U*fe(U)    ','ufue'//ctime(1:12) )
  call grxy(xv,work3,nvx,7200+idc,1,1 &
       ,cu,'   Q/Qtot      ','qoqt'//ctime(1:12) )
  !     reset absorbed particle counter
  !     nesc = 1
end subroutine uescd


