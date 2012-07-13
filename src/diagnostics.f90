
!  ==========================
!
!   Main I/O diagnostic routine
!
!  ==========================
!
      subroutine diagnostics
      use bopsvars
      implicit none

      integer :: imovc, lct, lctime, l, i
      real :: tcurrent

      character :: ct*40,cbl*16
      data imovc/1/
      save imovc


!  convert cpu time to char variable
      cbl='                '
      if (iunits.eq.2) then
         tcurrent = itime*dt/tcfs*ttrans
         call chr(tcurrent,1,ct,lct)
         lctime=12-lct
         ctime='t_{fs}='//ct(1:lct)

      else if (iunits.eq.1) then
!  timestamp in wp^-1
         tcurrent = itime*dt/tconv*ttrans
         call chr(tcurrent,1,ct,lct)
         lctime=12-lct
         ctime='w_pt='//ct(1:lct)

      else
!  timestamp in w0^-1
         tcurrent = itime*dt*ttrans
         call chr(tcurrent,0,ct,lct)
         lctime=12-lct
         ctime='w_0t='//ct(1:lct)

      endif



!  graphical snapshots every igr timesteps
      if (mod(itime,igr).eq.0) then
      write (*,*) 'Writing out snapshots at t=',tcurrent
      call gsnap
      endif


!  2D cumulative plots every igmovie timesteps, starting
!  ncyc laser  cycles before graphics dump
!      if (igmovie.eq.-999) then
        if (itime.eq.0) then
          call chr(1.*imovc,0,ct,l)
           open(75,file='idens'//ct(1:l)//'.2D')
          open(76,file='edens'//ct(1:l)//'.2D')
          open(77,file='Bz'//ct(1:l)//'.2D')
          open(78,file='Jy'//ct(1:l)//'.2D')
          open(79,file='Jyzoom'//ct(1:l)//'.2D')
          call blank
          write(15,*) 'Dimensions of 2D plots (x,t): ',(xcur2-xcur1)/dx/igx2d &
      ,'  x  ',(iplot2d-xm1/dt)/igmovie
          write(15,*) 'Plot area PA=',xcur1,xcur2,' grid points ',xcur1/dx/igx2d,xcur2/dx/igx2d

          write(15,*) 'Start recording 2D plots at timestep ',xm1/dx
          write(*,*) 'Dimensions of 2D plots (x,t): ',(xcur2-xcur1)/dx/igx2d &
      ,'  x  ',(iplot2d-xm1/dt)/igmovie
               write(15,*) 'Zoom: ',(xcur2-xcur1)/dx+1 &
      ,'  x  ',itav*ncyc/igmovie
          call blank
        else if (itime.gt.0.and.mod(itime,iplot2d).eq.0) then
!  close/reopen 2D data files at graphics dump time
          if (itime.gt.0) then
            do i=0,3
               close(75+i)
            end do
            imovc = imovc+1
          endif
!          if (itime+iplot2d.le.nt) then
!            call chr(1.*imovc,0,ct,l)
!              open(75,file='idens'//ct(1:l)//'.2D')
!          open(76,file='edens'//ct(1:l)//'.2D')
!          open(77,file='Bz'//ct(1:l)//'.2D')
!          open(78,file='Jy'//ct(1:l)//'.2D')
!          open(79,file='Jyzoom'//ct(1:l)//'.2D')
!          endif
!        else  if (mod(itime,iplot2d).ge.iplot2d-ncyc*itav) then
!  calculate 1D profiles at time t and dump to file
        else if (mod(itime,igmovie).eq.0 .and. itime*dt.gt.xm1) then
	    call movie
        endif
!      endif

!  plots at special times in isp() eg: phase-space movie
!      call splot

!  FT sampling and dump reflected em waves
      call emstore

!      if (nt.lt.10) return
      call avlas
      call fieldave
      call rhoiav
      call vdist
      call udist
      call uidist
      call uescd
!      call f2v
!     call uescdi
!  increment snapshot counter
      if (mod(itime,igr).eq.0 .and. itime.ne.0) idc=idc+1
!      call gcon
 !     call currft

!  particle tracking (labels read from file)
      if (itropt.eq.2) then
        call stotrk

!  track particles just exiting plasma
      else if (itropt.eq.3) then
        call phatrk
      endif

!  non cycle-averaged history O/P
      if (mod(itime,itc).eq.0) then
       call energy
       call surfie
       ntc=ntc+1
      endif



      end
