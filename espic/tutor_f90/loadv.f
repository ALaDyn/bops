!  ======================================
!
!   Set up particle velocity distribution
!
!  ======================================

subroutine loadv

include 'es.h'

 iseed = 76523; idum1 = 137; idum2 = 45126

 do l=1,ne

!  inverted 2v-distribution - amplitude

    vm=vte*(-2.*alog((l-0.5)/ne))**0.5

!  random angle
    rs = ran(iseed)
    theta=2*pi*rs

!  x-component of v
    vx(l)=vm*sin(theta)

 end do


!  scramble particle indicies to remove correlations
!  between x and vx


!        call warn('Rand. particles with clock time')
!	write(15,*) time()
!        idum1 = -time()
!	idum2 = -time()-23
!	idum3 = -time()-47

!  exclude odd one out

      n1=ne
      if (mod(ne,2).ne.0) then
	n1=n1-1
      endif

      do i=1,n1
	j=n1*ran(idum1)+1

	temp1 = vx(i)        ! switch i,j
	vx(i) = vx(j)
	vx(j) = temp1

      end do

end
