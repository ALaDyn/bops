
!     ==========

subroutine pout
  use bopsvars
  implicit none


  if (mod(itime,iout).ne.0.or.itime.eq.0) return
  !      call blank

  if (iunits.eq.2) then
     write (15,101) 'timestep',itime,' w0t =',itime*dt*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'
     write (6,101) 'timestep',itime,' w0t =',itime*dt*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'

  else if (iunits.eq.1) then
     write (15,101) 'timestep',itime,' wpt =',itime*dt/tconv*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'
     write (6,101) 'timestep',itime,' wpt =',itime*dt/tconv*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'

  else
     write (15,101) 'timestep',itime,' w0t =',itime*dt*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'
     write (6,101) 'timestep',itime,' w0t =',itime*dt*ttrans,' t/fs =' &
          ,itime*dt/tcfs*ttrans,' (',100*itime*dt/trun,'% runtime)'
  endif
101 format(a,i7,3x,a,f8.1,3x,a,f8.1,3x,a,f5.1,a)

end subroutine pout
