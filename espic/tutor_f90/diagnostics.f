
!  ============================================
!
!    Diagnostics - graphics and printed output
!
!  ============================================

subroutine diagnostics

include 'es.h'


!  write run information to terminal

 if (mod(itime,iout).eq.0) then
   write (6,*) 'timestep:', itime
 endif

!  do graphics snapshots

 if (mod(itime,igraph).eq.0) then
    call plots
 endif

!  write out time-dep. quantities to file
 if (mod(itime,ihist).eq.0) then
    call histories
 endif

end
