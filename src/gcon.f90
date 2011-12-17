
!     ==========

subroutine gcon
  use bopsvars
  implicit none

  integer ncx, ncy, jc, jcon, i, icon

  !  contours in target frame

  if (itime.eq.0) return
  ncx=nx/igx2d
  if (mod(nx,igx2d).ne.0) ncx=ncx+1
  ncy=itav/isgy
  if (mod(itav,isgy).ne.0) ncy=ncy+1
  if (mod(itime,igr).ge.igr-itav) then
     jc=mod(itime+itav,igr)
     jcon=(jc+isgy-1)/isgy

     do i=1,nx,igx2d
        icon=(i+igx2d-1)/igx2d
        !        cbz(icon,jcon)=gam0*(bz(i)+vy0*ex(i))
        !        cph(icon,jcon)=phi(i)
     end do

  endif
  if (mod(itime,igr).eq.0) then
     !        call grcon(lun,cbz,ncx,ncy,5000+idc)
     !        call grcon(lun,cph,ncx,ncy,5100+idc)
  endif
end subroutine gcon













