
!  =======================================================
!
!    EM/ES Power spectra:  TE mode
!
!   NB: Max freq  wmax = pi/(ift*dt)
!       Resolution  dw = 2*pi/(nf*dt*ift)
!
!  =======================================================

subroutine teft
  use bopsvars
  implicit none

  integer incall, ifbm, i, k, j, ifreq, nfg, igf, nhmax, iwband
  integer iw1, iw, ihar, nfg1, nfzoom
  real*8 uewm, wex, ze, uwmax, uwmin, wmax, wmin, dom, ureal, ufreq
  real*8 wmstep, wmzoom, uwstep
  complex cx(nf)
  real :: rw1(nf), rw2(6*nf+150)  ! FT work arrays - single precision
  integer :: iwk(6*nf+150)    !
  real*8 :: work1(0:nf), work2(0:nf), work3(0:nf)
  real*8 Pback(1000),Pfor(1000)
  save incall
  data incall/1/

  !  binning
  ifbm = max(1.d0,0.1/dw)
  if (ifbin.gt.ifbm) then
     ifbin=ifbm
     write (15,*) 'NB: Frequency binning reset to',ifbin
     write (6,*) 'NB: Frequency binning reset to',ifbin
  endif

  if (nt.lt.10) return

  !  forward em wave
  do i=1,nf
     rw1(i)=ffl(i)
     rw2(i)=0.
  end do


  !  IMSL FFT
  call fftrc(rw1,nf,cx,iwk,rw2)

  !  get remainder of coefficients
  do k=2,nfo2
     cx(nf+2-k) = conjg(cx(k))
  end do

  !  Power density
  do iw=2,nfo2
     work1(iw)=(abs(cx(iw))**2 + abs(cx(nf+2-iw))**2)/nf**2
  end do
  work1(1) = abs(cx(1))**2/nf**2

  !  rearrange power spectrum into bins
  uewm = 0.
  do j=1,nfo2/ifbin
     work3(j) = dw*ifbin*(j-0.5)
     ifreq = ifbin*(j-1)
     work2(j)=0.
     do i=1,ifbin
        work2(j)=work2(j)+work1(ifreq+i)/ifbin
     end do
     uewm=max(work2(j),uewm)
  end do

  wex=ze(uewm)+1
  uwmax=5*10**wex
  uwmin=10**(wex-6.)
  wmstep=.5
  uwstep = 10.

  work1(0)=uwmin
  wmax=omegm
  wmin=0.
  dom=5.
  nfg=nfo2/ifbin
  igf=1+nfg/1000

  call grxy(work3,work2,nfg,70000+incall,1,-2 &
       ,' w/w0          ',' F+(w)         ','tefs'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwstep)


  !  Power ratios of transmitted harmonics: average over bandwidth wband

  nhmax = omegm
  iwband = ifbin*2
  iw1 = 1./dw

  Pfor(1)=0.

  do iw = iw1-iwband,iw1+iwband
     Pfor(1) = Pfor(1) + work1(iw)
  end do

  if (Pfor(1).eq.0.) Pfor(1)=1.

  do ihar=2,nhmax
     iw1 = ihar/dw
     Pfor(ihar) = 0.
     do iw = iw1-iwband,iw1+iwband
        Pfor(ihar) = Pfor(ihar) + work1(iw)
     end do
  end do



  !   backward em wave
  !   ----------------


  do i=1,nf
     rw1(i)=fb0(i)
     cx(i)=0.
     rw2(i)=0.
  end do

  call fftrc(rw1,nf,cx,iwk,rw2)

  !  get remainder of coefficients
  do k=2,nfo2
     cx(nf+2-k) = conjg(cx(k))
  end do

  do iw=2,nfo2
     work1(iw)=(abs(cx(iw))**2 + abs(cx(nf+2-iw))**2)/nf**2
  end do

  work1(1) = abs(cx(1))**2/nf**2

  !  rearrange power spectrum into bins
  uewm = 0.

  do j=1,nfo2/ifbin
     work3(j) = dw*ifbin*(j-0.5)
     ifreq = ifbin*(j-1)
     work2(j)=0.
     do i=1,ifbin
        work2(j)=work2(j)+work1(ifreq+i)/ifbin
     end do
     uewm=max(work2(j),uewm)
  end do

  !  check Parseval Theorem
  ureal=0.
  ufreq=0.

  do j=1,nf
     ureal = ureal + fb0(j)**2
     ufreq = ufreq + cx(j)*conjg(cx(j))/nf
  end do

  call blank
  write (15,*) 'Parseval check:'
  write (15,*) 'time:',ureal
  write (15,*) 'freq:',ufreq
  call blank


  wex=ze(uewm)
  uwmax=5*10**wex
  uwmin=10**(wex-10.)
  work1(0)=uwmin
  wmax=omegm
  wmzoom=2.5
  wmin=0.
  dom=5.
  nfg1=nfo2/ifbin

  !  plot spectrum
  call grxy(work3,work2,nfg1,71000+incall,1,-2 &
       ,' w/w0          ',' F-(w)         ','tebs'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwstep)

  !  zoom of backscatter spectra (fund and 2nd harmonic)
  nfzoom=3.5/dw
  call grxy(work3,work2,nfzoom,71500+incall,1,-2 &
       ,' w             ',' F-(w)         ','tebz'//ctime(1:12) &
       ,wmin,wmzoom,wmstep,uwmin,uwmax,uwstep)


  !  Power ratios of reflected harmonics: average over bandwidth wband

  Pback(1)=0.
  iw1 = 1./dw

  do iw = iw1-iwband,iw1+iwband
     Pback(1) = Pback(1) + work1(iw)
  end do

  if (Pback(1).eq.0.) Pback(1)=1.

  do ihar=2,nhmax
     iw1 = ihar/dw
     Pback(ihar) = 0.
     do iw = iw1-iwband,iw1+iwband
        Pback(ihar) = Pback(ihar) + work1(iw)
     end do
  end do



  !  plasma wave at critical
  !  -----------------------


  do i=1,nf
     rw1(i)=fex(i)
     rw2(i)=0.
  end do

  call fftrc(rw1,nf,cx,iwk,rw2)

  !  get remainder of coefficients
  do k=2,nfo2
     cx(nf+2-k) = conjg(cx(k))
  end do

  do iw=2,nfo2
     work1(iw)=(abs(cx(iw))**2 + abs(cx(nf+2-iw))**2)/nf**2
  end do

  work1(1) = abs(cx(1))**2/nf**2

  !  rearrange power spectrum into bins
  uewm = 0.
  do j=1,nfo2/ifbin
     work3(j) = dw*ifbin*(j-0.5)
     ifreq = ifbin*(j-1)
     work2(j)=0.
     do i=1,ifbin
        work2(j)=work2(j)+work1(ifreq+i)/ifbin
     end do
     uewm=max(work2(j),uewm)
  end do
  wex=ze(uewm)
  uwmax=5*10**wex
  uwmin=10**(wex-8.)
  work1(0)=uwmin
  wmax=omegm
  wmin=0.
  dom=5.


  call grxy(work3,work2,nfg,72000+incall,1,-2 &
       ,'    w/w0       ','       Ex(w)   ','escr'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwstep)

  !  Electron current at critical


  do i=1,nf
     rw1(i)=fjy(i)
     rw2(i)=0.
  end do

  call fftrc(rw1,nf,cx,iwk,rw2)

  !  get remainder of coefficients
  do k=2,nfo2
     cx(nf+2-k) = conjg(cx(k))
  end do

  uewm=0.
  do iw=2,nfo2
     !  IMSL
     work1(iw)=(abs(cx(iw))**2 + abs(cx(nf+2-iw))**2)/nf**2
     uewm=max(work1(iw),uewm)
  end do

  work1(1) = abs(cx(1))**2/nf**2

  !  rearrange power spectrum into bins
  uewm = 0.
  do j=1,nfo2/ifbin
     work3(j) = dw*ifbin*(j-0.5)
     ifreq = ifbin*(j-1)
     work2(j)=0.
     do i=1,ifbin
        work2(j)=work2(j)+work1(ifreq+i)/ifbin
     end do
     uewm=max(work2(j),uewm)
  end do
  wex=ze(uewm)
  uwmax=5*10**wex
  uwmin=10**(wex-8.)
  work1(0)=uwmin
  wmax=omegm
  wmin=0.
  dom=5.

  call grxy(work3,work2,nfg,73000+incall,1,-2 &
       ,'    w/w0       ','       Jy(w)   ','jycr'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwstep)

  !  write harmonic powers to file

  if (cpolzn.eq.'P') then
     !       write (15,*) 'Power ratios'
     !       write (15,*) '------------'
     !       call blank
     open(61,file='powers_te.dat')
     write (61,*) '!  n, P_back(n), P_for(n)'
     do ihar=1,nhmax
        !          write (15,'(i6,2(1pe13.2))') ihar,Pback(ihar),Pfor(ihar)
        write (61,'(i6,2(1pe13.2))') ihar,Pback(ihar),Pfor(ihar)
     end do
     close(61)
  endif

  incall=incall+1
end subroutine teft
