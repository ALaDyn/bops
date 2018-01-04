
! ================
!
!    Distribution of escaped ions
!
! ================

subroutine uescdi
  use bopsvars
  implicit none

  integer i, l, i0, iu
  real*8 ukm, du, ua, sumf, uht, u, fmin, fmax
  real*8, dimension(0:nvx+1) :: work1,work2, work3

  if (mod(itime,igr).ne.0 .or. itime.eq.0 .or. nesci.eq.0 ) return

  !  energy distribution of escaped particles - in LAB FRAME
  do i=1,nvx
     work2(i)=0.
     work1(i)=0.
  end do

  ukm=0.
  do l=1,nesci
     ukm=max(ukm,uesci(l))
  end do

  du=umevmax/.511/miome/nvx

  do l=1,nesci
     ua=uesci(l)/du
     iu = ua+1
     i0=min(max(1,iu),nvx)
     work1(i0)=work1(i0)+1
  end do

  !  normalise to nesci
  sumf=0.

  do i=1,nvx
     sumf=sumf+work1(i)*du
  end do

  if (sumf.eq.0) sumf=1

  do i=1,nvx
     work2(i)=work1(i)/sumf*nesci
  end do

  !  1st moment - gives absorption
  !  k.e. distn in work2, cumulative in work3
  work3(0)=0.

  !  renorm to total em energy in lab frame
  uht=0.
  do i=1,nvx
     u=du*i-du/2
     work2(i)=work2(i)*u*du*mi**2
     uht=uht+work2(i)
  end do

  !  hot ion temp
  uht=uht/nesci

  !   thermal energy distn
  !  of particles escaped from rh boundary
  !  plotted against k.e.  gamma-1  normalised to mc**2
  !
  fmin=work1(1)
  fmax=work1(1)
  do i=1,nvx
     if (work1(i).ne.0) fmin=min(fmin,work1(i))
     if (work1(i).ne.0) fmax=max(fmax,work1(i))
  end do

  do i=1,nvx
     if (work1(i).le.0) work1(i)=fmin
     u=i*du-du/2
     xv(i)=u*.511*miome
  end do

  call r0form('Thot  ',uht*511*miome,'f12.2')
  call grxy(xv,work1,nvx,6900+idc,1,2 &
       ,'   Uesc(MeV)   ','     fi(Uesc)  ','fuis'//ctime(1:12) &
	     ,-1.0,1.0,.2,1.,1e3,10.)
end subroutine uescdi


















