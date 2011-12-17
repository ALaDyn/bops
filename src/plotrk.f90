!  ========================================================
!
!       Plot tracked particle orbits
!
!  ========================================================

subroutine plotrk
  use bopsvars
  implicit none

  real*8 tstart, tend, xmax, xmin, uboost, uymin, uymax
  integer i,j, lsym
  character chead*80,cp*80,csym*40,cxmin*80,cxmax*80,cymin*80 &
       ,cymax*80,cmima*60,cfile*80,csnap*80
  real*8 :: work1(ittrk), work2(ittrk)

  if (ntrack.eq.0.or.itropt.eq.3) return

  !  cplot file header
  cfile='track.cp'
  open(60,file=cfile)

  !  initialise graphics filter

  write (60,'(a/a/a)') 'zi psfilter track.ps ','cd','cd pic'

  !  real space: x vs t
  !  ------------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'
  tstart = dt*itstart
  tend = itend*dt
  xmax = gam0*(xsol+2.)
  xmin = gam0*(xsol-2.)
  write (60,'(a,4f13.4)') 'ra ',tstart,min(tend,nt*dt),xmin,xmax

  !  Title

  write (60,'(a,9(/a))') &
       'tx' &
       , 'real space orbits' &
       , ' t' &
       , ' ' &
       , ' x' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .' &
       , 'zalt',' '


  do i=1,ittrk
     twork(i)=i*dt*itsk
  end do

  do i=1,ntrack
     do j=1,ittrk
        work1(j)=gam0*xtrk(i,j)
     end do
     call gtrak(twork,work1,ittrk,i,1,0 &
          ,'       t       ','       x       ','xvst            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'xvst/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     write (60,'(a/a)') 'zp',' '

  end do
  !  new page
  write (60,'(a/a)') 'zx','re'


  !  real space: y vs x
  !  ------------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'
  tstart = dt*itstart
  tend = itend*dt
  xmax = gam0*(xsol+2.)
  xmin = gam0*(xsol-5.)
  write (60,'(a,2f13.4)') 'rax ',xmin,xmax

  !  Title

  write (60,'(a,8(/a))') &
       'tx' &
       , 'real space orbits' &
       , ' x' &
       , ' ' &
       , ' y' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .'
  !     :, 'z1ta',' '


  do i=1,ntrack
     do j=1,ittrk
        work1(j)=gam0*xtrk(i,j)
        work2(j)=gam0**2*(xytrk(i,j) - vy0*twork(j))
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','       y       ','yvsx            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'yvsx/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     if (i.eq.1) then
        write (60,'(a/a)') 'zaltp',' '
     else
        write (60,'(a/a)') 'zp',' '
     endif

  end do
  !  new page
  write (60,'(a/a)') 'zx','re'

  !  phase space: ux vs x
  !  --------------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'

  write (60,'(a,4f13.4)') 'ra ',xmin,xmax,-3*a0,3*a0

  !  Title

  write (60,'(a,9(/a))') &
       'tx' &
       , 'phase space (px)' &
       , ' x' &
       , ' ' &
       , ' px' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .' &
       , 'zalt',' '


  do i=1,ntrack
     do j=1,ittrk
        work2(j)=uxtrk(i,j)
        work1(j)=gam0*xtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      ux       ','uvsx            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'uvsx/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     write (60,'(a/a)') 'zp',' '

  end do
  !  new page
  write (60,'(a/a)') 'zx','re'


  !   uy vs x
  !  --------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'

  uboost=vy0*gam0
  uymin = uboost-2*a0
  uymax = uboost+2*a0
  write (60,'(a,4f13.4)') 'ra ',xmin,xmax,uymin,uymax

  !  Title

  write (60,'(a,9(/a))') &
       'tx' &
       , 'phase space (py)' &
       , ' x' &
       , ' ' &
       , ' py' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .' &
       , 'zalt',' '

  do i=1,ntrack
     do j=1,ittrk
        work2(j)=uytrk(i,j)
        work1(j)=gam0*xtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      uy       ','uyvx            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'uyvx/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     write (60,'(a/a)') 'zp',' '
  end do

  !  new page
  write (60,'(a/a)') 'zx','re'


  !  Acceleration vs x
  !  --------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'

  !      write (60,'(a,2f13.4)') 'rax ',tstart,min(tend,nt*dt)
  write (60,'(a,4f13.4)') 'ra ',xmin,xmax,-a0,2*a0

  !  Title

  write (60,'(a,8(/a))') &
       'tx' &
       , 'acceleration' &
       , ' x' &
       , ' ' &
       , ' ax' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .'
  !     :, 'zalt',' '

  do i=1,ntrack
     do j=1,ittrk
        work1(j)=gam0*xtrk(i,j)
        work2(j)=axtrk(i,j)
     end do
     call gtrak(work1,work2,ittrk,i,1,0 &
          ,'       x       ','      ax       ','axvx            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'axvx/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     if (i.eq.1) then
        write (60,'(a/a)') 'zaltp',' '
     else
        write (60,'(a/a)') 'zp',' '
     endif
  end do
  !  new page
  write (60,'(a/a)') 'zx','re'


  !  Acceleration vs t
  !  --------------

  write (60,'(5(/a))') 'zw','re','tu','ft 3','wi 5. 5. 11. 7.'

  write (60,'(a,4f13.4)') 'ra ',tstart,min(tend,nt*dt),-a0,2*a0

  !  Title

  write (60,'(a,8(/a))') &
       'tx' &
       , 'acceleration' &
       , ' t' &
       , ' ' &
       , ' ax' &
       , ' ' &
       , 'sy L' & !  lines
       , 'ty +6 +6 .'
  !     :, 'zalt',' '

  do i=1,ntrack
     do j=1,ittrk
        work2(j)=axtrk(i,j)
     end do
     call gtrak(twork,work2,ittrk,i,1,0 &
          ,'       t       ','      ax       ','axvt            ')

     !  write filename to cplot header
     call chr(i*1.,0,csym,lsym)
     csnap = 'axvt/'//csym(1:lsym)
     write (60,'(a)') 'gd 2 '//csnap(1:5+lsym)
     !  display
     if (i.eq.1) then
        write (60,'(a/a)') 'zaltp',' '
     else
        write (60,'(a/a)') 'zp',' '
     endif
  end do

  !  close filter
  write (60,'(a/a)') 'ex','y'

  close(60)

end subroutine plotrk
