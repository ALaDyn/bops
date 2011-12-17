
!  =======================================================
!
!    EM/ES Power spectra for TM mode
!
!   NB: Max freq  wmax = pi/(ift*dt)
!       Resolution  dw = 2*pi/(nf*dt*ift)
!
!  =======================================================

subroutine tmft
  use bopsvars
  implicit none

  integer incall, ifbm, i, k, iw, j, ifreq, nfg, igf, nhmax, iwband
  integer iw1, ihar, nfg1
  real*8 uewm, wex, ze, uwmax, uwmin, wmax, wmin, dom
  real*8 ureal, ufreq
  complex cx(nf)
  real :: rw1(nf), rw2(6*nf+150)  ! FT work arrays
  integer :: iwk(6*nf+150)    !
  real*8 :: work1(0:nf), work2(0:nf), work3(0:nf)

  real*8 P(1000)
  save incall
  data incall/1/

  !  binning
  ifbm = max(1.d0,0.1/dw)
  if (ifbin.gt.ifbm) then
     ifbin=ifbm
     write (15,*) 'Used binning freq',ifbin
  endif

  if (nt.lt.10) return

  !  forward em wave
  do i=1,nf
     rw1(i)=gfl(i)
     rw2(i)=0.
  end do


  !  FFT
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
  uwmin=10**(wex-5.)
  work1(0)=uwmin
  wmax=omegm
  wmin=0.
  dom=5.
  nfg=nfo2/ifbin
  igf=1+nfg/1000

  call grxy(work3,work2,nfg,75000+incall,1,-2 &
       ,' w/w0          ',' G+(w)         ','tmfs'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwmax/10.)


  !   backward em wave
  !   ----------------

  !      call filter1(gb0,nf)

  do i=1,nf
     rw1(i)=gb0(i)
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
     ureal = ureal + gb0(j)**2
     ufreq = ufreq + cx(j)*conjg(cx(j))/nf
  end do

  call blank
  write (15,*) 'Parseval check:'
  write (15,*) 'time:',ureal
  write (15,*) 'freq:',ufreq
  call blank

  !  Power ratios of harmonics: average over bandwidth wband

  nhmax = omegm
  iwband = ifbin*2
  iw1 = 1./dw
  P(1)=0.

  do iw = iw1-iwband,iw1+iwband
     P(1) = P(1) + work1(iw)
  end do

  if (P(1).eq.0.) P(1)=1.

  !      write (15,*) 'TM Power ratios'
  !      write (15,*) '--------------'
  !      call blank

  open(60,file='powers_tm.dat')
  write (60,'(i6,1pe13.2)') 1,P(1)

  do ihar=2,nhmax
     iw1 = ihar/dw
     P(ihar) = 0.
     do iw = iw1-iwband,iw1+iwband
        P(ihar) = P(ihar) + work1(iw)
     end do
     !        write (15,'(i6,1pe13.2)') ihar,P(ihar)
     write (60,'(i6,1pe13.2)') ihar,P(ihar)
  end do
  close(60)
  call blank

  wex=ze(uewm)
  uwmax=5*10**wex
  uwmin=10**(wex-7.)
  work1(0)=uwmin
  wmax=omegm
  wmin=0.
  dom=5.
  nfg1=nfo2/ifbin

  call grxy(work3,work2,nfg1,76000+incall,1,-2 &
       ,' w/w0          ',' G-(w)         ','tmbs'//ctime(1:12) &
       ,wmin,wmax,dom,uwmin,uwmax,uwmax/10.)

  incall=incall+1
end subroutine tmft
