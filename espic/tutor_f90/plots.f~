
!  =========================================
!
!   Graphical snapshots
!
!  =========================================

subroutine plots
include 'es.h'

 real :: xgrid(0:nxm)             ! grid work array

 character*40 cfile
 character*1 cid

 integer :: isnap                ! counts number of calls to routine 'plots'

 data isnap/0/ 
 save isnap                      ! retains value of isnap for next call

 xgrid(0:nx) = (/(i*dx,i=0,nx)/)   ! set up x-values for grid plots


 cid = achar(mod(isnap,10)+48)   !  convert counter to ASCII character '0-9'

!  electron px-x phase space
 
 cfile = 'phase'//cid//'.data'    ! build filename from components
 open (50,file=cfile(1:11))      !  open data file

 write (50,'(2f10.5)') (x(i),vx(i),i=1,ne)
 close(50)



!  electron density

 cfile = 'edens'//cid//'.data'    ! build filename from components
 open (50,file=cfile(1:11))      !  open data file

 write (50,'(2f10.5)') (xgrid(j),-rhoe(j),j=0,nx)
 close(50)


!  electrostatic field

 cfile = 'field'//cid//'.data'    ! build filename from components
 open (50,file=cfile(1:11))      !  open data file

 write (50,'(2f10.5)') (xgrid(j),Ex(j),j=0,nx)
 close(50)


 isnap = isnap + 1
     
end
