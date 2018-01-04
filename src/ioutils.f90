!  ====================
!  **  i/o routines  **
! ====================

!i0prnt
      subroutine i0prnt(a,i0)
      parameter(ichan=15)
      character a*6
      write(15,10) a,i0
!      write(6,10) a,i0
   10 format(1x,a6,i8)
      return
      end


!     ==========
!i0prn6
      subroutine i0prn6(a,i0)
      character a*6
      write(6,10) a,i0
   10 format(1x,a6,i8)
      return
      end


!     ==========

      subroutine blank6
      write (6,10)
  10  format(/)
      end

!     ==========


      subroutine blank
      parameter(ichan=15)
      write (ichan,10)
  10  format(/)
      end


!r0form
!
!  formatted real number

      subroutine r0form(a,r0,cfm)
      real*8 r0
      character a*6,cform*20
      character (len=*) cfm
      if (cfm(1:1).eq.'f') then
      cform='(1x,a6,'//cfm(1:5)//')'
      lfm=13
      else
      cform='(1x,a6,1pe18.8)'
      lfm=15
      endif
      write(15,cform(1:lfm)) a,r0
!      write(6,cform(1:lfm)) a,r0
      return
      end

      subroutine r0form6(a,r0,cfm)
      real*8 r0
      character a*6,cform*20
      character (len=*) cfm
      if (cfm(1:1).eq.'f') then
      cform='(1x,a6,'//cfm(1:5)//')'
      lfm=13
      else
      cform='(1x,a6,1pe18.8)'
      lfm=15
      endif
      write(6,cform(1:lfm)) a,r0
      return
      end

      subroutine ruform(a,r0,cfm,u)
      real*8 r0
      character a*6,cform*25,u
      character (len=*) cfm
      if (cfm(1:1).eq.'f') then
      cform='(1x,a6,'//cfm(1:5)//',1x,a6)'
      lfm=19
      else
      cform='(1x,a6,1pe18.8,1x,a6)'
      lfm=21
      endif
      write(15,cform(1:lfm)) a,r0,u
      return
      end

!r0prn6
      subroutine r0prn6(a,r0)
      character a*6
      real r0
      write(6,10) a,r0
   10 format(1x,a6,10(1pe12.4))
      return
      end

!d0prn6
      subroutine d0prn6(a,r0)
      character a*6
      real*8 r0
      write(6,10) a,r0
   10 format(1x,a6,10(1pe12.4))
      return
      end

!r0prnt
      subroutine r0prnt(a,r0)
      character a*6
      real r0
      write(15,10) a,r0
   10 format(1x,a6,10(1pe12.4))
      return
      end

!d0prnt
      subroutine d0prnt(a,r0)
      character a*6
      real*8 r0
      write(15,10) a,r0
   10 format(1x,a6,10(1pe12.3))
      return
      end


!r1prnt
      subroutine r1prnt(a,r1,m)
      character a*6
      dimension r1(m)
      write(15,10) a,(r1(i),i=1,m)
   10 format(1x/1x,a6/(1x,6(1pe12.3)))
      return
      end


!i1prnt
      subroutine i1prnt(a,i1,m)
      parameter(ichan=15)
      character a*6
      dimension i1(m)
      write(ichan,10) a,(i1(i),i=1,m)
   10 format(1x,a6/(1x,8i8))
      return
      end


!i3prnt
      subroutine i1prnt3(a,i1,m)
      parameter(ichan=15)
      character a*6
      dimension i1(m)
      write(ichan,10) a,(i1(i),i=1,m)
   10 format(1x/1x,a6/(1x,3i8))
      return
      end


      subroutine i2list(a,i1,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      dimension i1(m1:m2)
      write(ichan,10) a,(i,i1(i),i=is,ie)
   10 format(1x/1x,a6/(1x,2i8))
      return
      end


      subroutine i3list(a1,i1,a2,i2,m1,m2,is,ie)
      parameter(ichan=15)
      character a1*6,a2*6
      dimension i1(m1:m2),i2(m1:m2)
      write(ichan,10) a1,a2,(i,i1(i),i2(i),i=is,ie)
   10 format(1x/12x,2(2x,a6)/(1x,3i8))
      return
      end

      subroutine r2list(a,r1,r2,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      real r1(m1:m2),r2(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,2(1pe12.3))))
      return
      end


      subroutine r3list(a,r1,r2,r3,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      real r1(m1:m2),r2(m1:m2),r3(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,3(1pe12.3))))
      return
      end


      subroutine d2list(a,r1,r2,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      double precision r1(m1:m2),r2(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,2(f12.3))))
      return
      end

      subroutine d3list(a,r1,r2,r3,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      double precision r1(m1:m2),r2(m1:m2),r3(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,3(f12.3))))
      return
      end


      subroutine r4list(a,r1,r2,r3,r4,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      integer r1
      dimension r1(m1:m2),r2(m1:m2),r3(m1:m2),r4(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),r4(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(2i8,3(f12.3))))
      return
      end


!c0prnt
      subroutine c0prnt(a,c0)
      complex c0
      real a
      c0r=real(c0)
      c0i=aimag(c0)
      write(6,10) a,c0r,c0i
   10 format(1x,a6,'   (',1pe12.3,',',1pe12.3,')')
      return
      end


!     c1prnt
      subroutine c1prnt(a,c1,m)
      complex c1
      real a
      dimension c1(m)
      dimension c1r(300),c1i(300)
      do  i=1,m
         c1r(i)=real(c1(i))
         c1i(i)=aimag(c1(i))
      end do
      write(15,10) a,(c1r(i),i=1,m)
      write(15,11) a,(c1i(i),i=1,m)
 10   format(1x/1x,a6,' (real part)'/(1x,6(1pe12.3)))
 11   format(1x/1x,a6,' (imag part)'/(1x,6(1pe12.3)))
      return
      end



!d1prnt
      subroutine d1prnt(a,r1,m)
      real a,r1
      dimension r1(m)
      write(15,10) a,(r1(i),i=1,m)
   10 format(1x/1x,a6/(1x,6(1pe12.3)))
      return
      end

!  ========================================
!
!      function   PHA
!
!  Returns phase of complex number -pi -> pi
!
!  ========================================

      real function pha(w)
      parameter (pi=3.1415926)
      complex w
      x=real(w)
      y=aimag(w)
      sx=sign(1.,x)
      sy=sign(1.,y)
      itx=1
      if (sx.eq.1.) itx=0
      if ((x.ne.0.).and.(y.ne.0)) then
       pha=atan(y/x)+itx*sy*pi
      else
       pha=0.
      endif
      return
      end

!  ========================================
!
!         function ZE
!
!   returns exponent of v
!
!  ========================================


      real function ze(v)
      implicit none
      integer ia
      real*8 v,q
      call manex(v,q,ia)
      ze=ia
      end

!  ========================================
!
!       MANEX
!
!   returns   mantissa  q
!             exponent  a
!
!  ========================================

      subroutine manex(z,q,ia)
      implicit none
      integer ia
      real*8 z,q,xlog,x
      real*8 b,s,a
      s=sign(1.d0,z)
      x=abs(z)
      if (x.eq.0) then
       q=0.
       a=0.
       return
      endif
      xlog=log10(x*1.000001d0)
      ia=xlog
      if (xlog.lt.0) ia=ia-1
      if (ia.ne.xlog) then
      b=xlog-ia
      else
      b=xlog
      endif
      q=s*10**b
      end


!  ========================================
!
!         ICSIZ
!
!  Returns length of character string
!
!  ========================================

      integer function icsiz(cm)
      character cm*40,ch*2
      l=0
20    l=l+1
      ch=cm(l:l)
      if ((ch.ne.'$').and.(l.lt.40)) goto 20
      if (ch.eq.' ') l=0
      icsiz=max(l-1,1)
      end


!     ==========

      subroutine warn(a)
      parameter(ichan=15)
      character a*30
      write (ichan,10) a
      write (6,10) a
  10  format(/1x,a30/1x,'------------------------------'/)
      end


!  ======================================
!
!            CHR
!
!    converts real number to character
!
!  ======================================

      subroutine chr(z,ndp,ch,l)
      implicit none
      real z
      integer ndp,l
      real za,e
      integer ia, ie, iz, idp, iunit, idigit, ic, is
      character ch*40,chnum*10
      data chnum/'0123456789'/
!
!   z     = real number
!   ndp   = # decimal points (0-7)
!   ch*40 = returned character string
!   l     = length of string
!
!  multiply up to remove decimal points
      za = max(1.,abs(z))
      ia=log10(za)
      e=float(ia)
      ie=10**ndp
      iz=int(abs(z)*ie+0.5)
!  add - sign for negatives
      is=sign(1.,z)
      l=0
      if (is.lt.0) then
      l=l+1
      ch(l:l)='-'
      endif
!  left of decimal point
      do idp=ia,0,-1
      l=l+1
      iunit=10**(idp+ndp)
      idigit=iz/iunit
      ic=max(min(idigit+1,10),1)
      ch(l:l)=chnum(ic:ic)
        iz=iz-idigit*iunit
      end do
!  decimal point
      if (ndp.gt.0) then
      l=l+1
      ch(l:l)='.'
      endif
!  right of decimal point
      do idp=-1,-ndp,-1
      l=l+1
      iunit=10**(idp+ndp)
      idigit=iz/iunit
      ic=max(min(idigit+1,10),1)
      ch(l:l)=chnum(ic:ic)
        iz=iz-idigit*iunit
      end do
      end


      function lench(c1)
      character (len=*) c1
      character c0*1
      l=80
  10  l=l-1
      c0=c1(l:l)
      if ((c0.eq.' ').and.(l.gt.1)) goto 10
      if (c0.eq.' ') l=0
      lench=l
      end
!     Random number scrambler
!     =======================
!
!     called with:
!     x=rano(iseed)
!
!     returns real number in the interval (0.0 - 1.0)
!     set iseed = -ve integer before using call
!
!     Example of calling routine:
!
!     subroutine coords
!     include 'common.h'
!
!     iseed1 = -11
!     iseed2 = -7
!
!
!     do i=1,n
!     x(i)=xlen*rano(iseed1)
!     y(i)=ylen*rano(iseed2)
!     end do
!
!
!     end




!
      real*8 function rano(idum)

      implicit none
      integer idum
      real*8  dseed, dum
      real*8 v(97), y
      real*8 genr32
      integer  iff, icall, j,i
      data iff,icall/0,0/
      save dseed, v, y, iff, icall

      if (idum.lt.0.or.iff.eq.0) then
       iff = 1
       dseed=abs(idum)*1.d0
       idum=1
       do j=1,97
            dum=genr32(dseed)
         end do

       do j=1,97
         v(j)=genr32(dseed)
         end do
       y=genr32(dseed)
      endif

!  next index - make sure we don't overstep array bounds if
!  generator returns a 0.0 or 1.0

      j=max(mod(1+int(97.*y),98),1)
      if(j.gt.97.or.j.lt.1) then
      write (6,*) 'Call: ',icall
        write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
      write (6,*) 'Random No. generator not initialised properly'
      write (6,*) 'dummy =',dum,' dseed=',dseed
      write (6,100) (i,v(i),i=1,97)
 100    format (i4,f10.6)
      stop
      endif
!  get next variate and generate new one to fill gap

      y=v(j)
      rano=y
      v(j)=genr32(dseed)
      icall = icall + 1

      return
      end

!  genr32
!
!-----------------------------------------------------------------------
!

!   purpose             - basic uniform (0,1) random number generator -
!                          32 bit

      real*8 function genr32(dseed)
!                                  specifications for arguments
      real*8   dseed
!                                  specifications for local variables
      real*8   d2p31m,d2p31
!                                  d2p31m=(2**31) - 1
!                                  d2p31 =(2**31)(or an adjusted value)
      data               d2p31m/2147483647.d0/
      data               d2p31 /2147483648.d0/
!                                  first executable statement
      dseed = dmod(16807.d0*dseed,d2p31m)
      genr32 = dseed / d2p31
      return
      end




