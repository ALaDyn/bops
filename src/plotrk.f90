!  ========================================================
!
!       Plot tracked particle orbits
!
!  ========================================================

subroutine plotrk
  use bopsvars
  implicit none

  real*8 tstart, tend, xmax, xmin, uboost, uymin, uymax
  integer i,j, lsym, itr, lct
  character cfmt*80,cp*80,csym*40,cxmin*80,cxmax*80,cymin*80 &
       ,cymax*80,cmima*60,cfile*80,csnap*40
  real*8 :: work1(ittrk), work2(ittrk)

  if (ntrack.eq.0.or.itropt.eq.3) return

 write (*,*) 'Writing out particle tracks'
 open(60,file='track_xvst.xy')
 open(61,file='track_uxvst.xy')
 open(62,file='track_uyvst.xy')
 open(63,file='track_uzvst.xy')

  !  real space: x vs t
  !  ------------------

  tstart = dt*itstart
  tend = itend*dt
  if (ioboost.eq.1) then
    xmax=xsol+2
    xmin=xsol-2
  else
    xmax = gam0*(xsol+2.)
    xmin = gam0*(xsol-2.)
  endif

  call chr(1.*ntrack,0,csnap,lct)
  cfmt = '(f12.3,'//csnap(1:lct)//'(f12.3))'
  write(*,'(a20)') cfmt

  ! Tracking data written as 
  !  t1, x(1,1), x(2,1), x(3,1), ... x(ntrack,1)
  !  t2, x(2,1), x(2,2) ....         x(ntrack,2)
  !  ..
  !  t_ittrk,  x(ittrk,1) ....       x(ntrack,ittrk)

  !  real space: x(t)
  !  ------------------
  do i=1,ittrk
     twork(i)=i*dt*itsk*ttrans
     write (60,cfmt) twork(i),(xtrk(itr,i),itr=1,ntrack)
     write (61,cfmt) twork(i),(uxtrk(itr,i),itr=1,ntrack)
     write (62,cfmt) twork(i),(uytrk(itr,i),itr=1,ntrack)
     write (63,cfmt) twork(i),(uztrk(itr,i),itr=1,ntrack)
  end do

  close(60)
  close(61)
  close(62)
  close(63)
 return


  !  phase space: ux vs x
  !  --------------------

  do i=1,ntrack
     do j=1,ittrk
        work2(j)=uxtrk(i,j)
        work1(j)=xtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      ux       ','uvsx            ')
  end do

 
  !  real space: y vs x
  !  ------------------


  tstart = dt*itstart
  tend = itend*dt
  xmax = gam0*(xsol+2.)
  xmin = gam0*(xsol-5.)


  do i=1,ntrack
     do j=1,ittrk
	  work1(j)=xtrk(i,j)
          work2(j)=ytrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','       y       ','yvsx            ')

  end do

 !   uy vs x
  !  --------------

  uboost=vy0*gam0
  uymin = uboost-2*a0
  uymax = uboost+2*a0


  do i=1,ntrack
     do j=1,ittrk
        work2(j)=uytrk(i,j)
        work1(j)=xtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      uy       ','uyvx            ')

   end do


  do i=1,ntrack
     do j=1,ittrk
        work1(j)=xtrk(i,j)
        work2(j)=axtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      ax       ','axvx            ')

  end do

  do i=1,ntrack
     do j=1,ittrk
        work2(j)=axtrk(i,j)
     end do
     call gtrak(twork,work2,ittrk,i,1,0 &
          ,'       t       ','      ax       ','axvt            ')

  end do

end subroutine plotrk
