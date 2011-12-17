
!  ======================
!
!   'Movie' snapshots,
!      frequency igmovie
!
!  ======================

      subroutine movie
      use bopsvars
      implicit none

      integer i, ngx, lcn, lcf, i1, i2
      character cform*40,cn*40
      real*8, dimension(0:nx+1) :: work1

      i1 = xcur1/dx
      i2 = xcur2/dx
      ngx = i2-i1+1

!  lab frame ion density
!  ---------------------

      do i=i1,i2
        work1(i)=(rhoi(i) - vy0*jyi(i))/gam0
      end do

      write (75,'((2i6,1pe15.4))') (i,itime,work1(i),i=i1,i2)


!  lab frame electron density: high res zoom
!  ---------------------

      do i=i1,i2
        work1(i)=abs(rhoe(i) - vy0*jye(i))/gam0
      end do

      write (76,'((2i6,1pe15.4))') (i,itime,work1(i),i=i1,i2)


!  lab frame B-field
!  ---------------------

      do i=i1,i2
        if (cpolzn.eq.'P') then
          work1(i)=(bz(i) + vy0*ex(i))
        else
          work1(i)=ez(i)
        endif
      end do

      write (77,'((2i6,1pe15.4))') (i,itime,work1(i),i=i1,i2)

!  lab frame transverse current
!  ----------------------------

      do i=i1,i2
        if (cpolzn.eq.'P') then
          work1(i)=((jye(i)-vy0*rhoe(i))/gam0)
        else
          work1(i)=jze(i)
        endif
      end do

      write (78,'((2i6,1pe15.4))') (i,itime,work1(i),i=i1,i2)



      end






