

!  =============================
!
!  Random scramble of velocities
!
!  =============================

subroutine scramb(i1,n)
  use bopsvars
  implicit none

  integer p, time, l, idum1, idum2, idum3, n, n1, i, j, k, i1
  real*8 rano
  real*8 temp1,temp2,temp3
  data idum1,idum2,idum3/-11,-7,-23/
  save idum1,idum2,idum3

  if (iran0.eq.0) then
     !  initialise seeds according to clock time: really random
     !        call warn('Rand. particles with clock time')
     !      write(15,*) time()
     !        idum1 = -time()
     !      idum2 = -time()-23
     !      idum3 = -time()-47
  else
     !  leave as defined above - enables repeat simulation with same particle config
     call warn('Rand. particles with preset seeds')
  endif

  n1=n
  !  exclude odd one out
  if (mod(n,2).ne.0) then
     n1=n1-1
!#### Anupam & Bin 2009/2010: for the exclusion one should also pass xo to xn
!CHECK(???)
     xn(i1+n-1)=xo(i1+n-1)
  endif

  !  scramble indices to remove correlation between (ux,uy,uz),x
  do i=1,n1
     p=i+i1-1
     j=n1*rano(idum1)+i1
     k=n1*rano(idum1)+i1
     l=n1*rano(idum1)+i1
     temp1=ux(p)
     temp2=uy(p)
     temp3=uz(p)
     ux(p)=ux(j)
     uy(p)=uy(k)
     uz(p)=uz(l)
     ux(j)=temp1
     uy(k)=temp2
     uz(l)=temp3
     xn(p)=xo(p)
  end do


end subroutine scramb

