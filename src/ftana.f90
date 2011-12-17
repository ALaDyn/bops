!  =======================================================
!
!    EM/ES Power spectra
!
!   NB: Max freq  wmax = pi/(ift*dt)
!       Resolution  dw = 2*pi/(nf*dt*ift)
!
!  =======================================================
      program ftana
      implicit none
      integer :: nffm
      real :: pi
      parameter(nffm=327680,pi=3.14159)
      real rw1(nffm),rw2(3*nffm+150),work1(3*nffm+250) &
      ,xwi(nffm),xlam(nffm),ampl(nffm),phase(nffm),efield(nffm)
      integer iwk(3*nffm+150)
      complex cx(nffm),cfilt(nffm)
      character cdir*80,cfile*80,cdf*80
      real :: P(10000), Ptotal, xt(nffm), ne(nffm), ni(nffm)
      real :: dt, dw, xlambda, tnorm
      real :: filter_min, filter_max, tspec, tstart, tend, tstart_n, tend_n
      integer :: ifilter_min, ifilter_max, nfr, nfilto2
      integer :: istart, iend, i, nf, nfo2, nfbin, ifbin, j, lf, n
      integer :: lench, iw, ifreq, iw1, iwband, ntc, nfg1, ihar, nhmax
      real :: pha, ureal, ufreq, ubinned, uewm, uwmin, uwmax, wmin, wmax, dom, wex, omegam
      real, external :: ze
      write (6,*) 'points, tstart, tend (in fs), filter min, max'
      read (5,*) nfr, tstart, tend, filter_min, filter_max
!      nfr=28720
!	nfr=33100
!      filter_min=60.
!      filter_max=200.
      omegam=1000.
      xlambda=0.8
      tnorm=xlambda/pi/0.6  ! conversion from t_fs to w0t
      tstart_n=tstart/tnorm
      tend_n=tend/tnorm
      ifbin=1
      write (6,*) 'filename:'
!      read (5,'(a)') cfile
      cfile='ey_back.t'
      lf=lench(cfile)
      cdf=cfile(1:lf)
      write(6,'(a)') cfile(1:lf)

!  open the read file
      open(50,file=cdf(1:lench(cdf) ) )
      read(50,*) (rw2(i),efield(i),i=1,nfr)
      rw2=rw2/tnorm
      tspec = rw2(nfr)-rw2(1)
      dt = rw2(2)-rw2(1) ! timestep

! Strip leading zeros
!      nf=0
!      do i=1,nfr
!	if (efield(i).ne.0) then
!	  nf=nf+1
!	endif
!     end do

      istart = tstart_n/dt  ! start and finish indices
      iend = tend_n/dt
      nf = iend - istart  ! # points in sample 
      tspec=tend_n-tstart_n ! duration of sample
      rw1(1:nf)=efield(1+istart:iend)  ! reduced sample
      nfo2=nf/2
      nfbin=nfo2/ifbin
      dw = 2*pi/tspec

      write (6,*) 'nf=',nf
      write (6,*) 'dt=',dt,' w0.tspec, tspec/fs',tspec,tspec*tnorm
      write (6,*) 'dw=',dw
      write (6,*) 'wmax=',dw*nfo2

      do i=1,nf
      cx(i)=0.
      end do

      call fftrc(rw1,nf,cx,iwk,rw2)

!  get remainder of coefficients
      do j=2,nfo2
        cx(nf+2-j) = conjg(cx(j))
      end do



      do iw=2,nfo2
      work1(iw)=(abs(cx(iw))**2 + abs(cx(nf+2-iw))**2)/nf**2
      end do

      work1(1) = abs(cx(1))**2/nf**2

!  rearrange power spectrum into bins

      do j=1,nfbin
         xwi(j) = dw*ifbin*(j-0.5)
         xlam(j) = xlambda/xwi(j)
         ifreq = ifbin*(j-1)
         ampl(j)=0.
         phase(j)=0.
         do i=1,ifbin
           iw=ifreq+i
           ampl(j)=ampl(j)+ work1(iw)/ifbin
           phase(j)=phase(j)+mod(pha(cx(iw))/ifbin,2*pi)
         end do
      end do


!  Power ratios of harmonics: average over bandwidth wband

      nhmax = omegam
      iwband = ifbin*2
      iw1 = 1./dw
      P(1)=0.
 	Ptotal = SUM(work1(1:nfo2))

      do iw = iw1-iwband,iw1+iwband
        P(1) = P(1) + work1(iw)
      end do

      if (P(1).eq.0.) P(1)=1.

      open(15,file='ratios.data')
      write (15,*) '! Power ratios'
      write (15,*) '! ------------'
      call blank

      do ihar=2,nhmax
        iw1 = ihar/dw
        P(ihar) = 0.
        do iw = iw1-iwband,iw1+iwband
          P(ihar) = P(ihar) + work1(iw)
        end do
      end do

! write out powers, normalising to integrated spectral power
      write (15,'(i6,1pe13.2)') (ihar,P(ihar)/Ptotal,ihar=1,nhmax)

      close(15)
      call blank

      wex=ze(uewm)
      uwmax=5*10**wex
      uwmin=10**(wex-7.)
      work1(1)=uwmin
      wmax=omegam
      wmin=0.
      dom=5.
      nfg1=nfo2/ifbin

!  amplitudes
!  need to rescale amplitudes for wavelength plot
      write(*,*) 'Writing ampl_w.spec, ampl_l.spec, phase.spec'
      open(20,file='ampl_w.spec')
      open(21,file='ampl_l.spec')
      write (20,'(2(1pe15.4))') (xwi(i),ampl(i),i=1,nfbin)
      write (21,'(2(1pe15.4))') (xlam(i),ampl(i),i=nfbin,1,-1)
      close(20)
      close(21)

!  phase
      open(20,file='phase.spec')
      write (20,'(2(1pe15.4))') (xwi(i),phase(i),i=1,nfbin)
      close(20)

!  check Parseval Theorem
      ureal=0.
      ufreq=0.

      do j=1,nf
        ureal = ureal + rw1(j)**2
        ufreq = ufreq + cx(j)*conjg(cx(j))/nf
      end do

      ubinned = 0.
      do i=1,nfbin
        ubinned=ubinned+ampl(i)*nf*ifbin
      end do

      call blank
      write (6,*) 'Parseval check:'
      write (6,*) 'time:  ',ureal
      write (6,*) 'freq:  ',ufreq
      write (6,*) 'binned:',ubinned
      write (6,*) 'Ptotal:',ptotal
      call blank

! Filter for as pulse
      ifilter_min = filter_min/dw
      ifilter_max = filter_max/dw
      write(6,*) 'filter min, max, nfreq',filter_min,filter_max,wmax
!      cfilt(ifilter_min:ifilter_max) = cx(ifilter_min:ifilter_max)
	cfilt(1:nf)=cx(1:nf)


!  Copy filtered coefficients
      do j=1,ifilter_min-1
       cfilt(j) = (0.,0.)
       cfilt(nf+2-j) = (0.,0.)
      end do
      do j=ifilter_max+1,nfo2
       cfilt(j) = (0.,0.)
       cfilt(nf+2-j) = (0.,0.)
      end do

!  Conjugate for inverse FT
     cfilt(1:nf) = conjg(cfilt(1:nf))

!  Inverse transform to get time series
      call fftcc(cfilt,nf,iwk,rw2)
      do i=1,nf
        work1(i)=real(conjg(cfilt(i))/n)
      end do
      open(20,file='aspulse.t')
      write (20,'(4(1pe15.4))') (i*dt*tnorm,work1(i),efield(i+istart),cfilt(i),i=1,nf)
      close(20)
      ntc=676
      open(20,file='nehi00.xy')
      open(21,file='nihi00.xy')
      open(30,file='ne.t')
      read (20,*) (xt(i),ne(i),i=1,ntc)
      read (21,*) (xt(i),ni(i),i=1,ntc)
      write(30,'(3(1pe15.5))') (xt(i)-13.3,ne(i),ni(i),i=1,ntc)  ! Offset to match efield
      close(20)
      close(21)
      close(30)
      end








