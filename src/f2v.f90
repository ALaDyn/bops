
!  ==================

subroutine f2v
  use bopsvars
  implicit none

  integer icycav, i, j, l, i1, i2, j1, j2
  real*8 umax, vxl, vyl, vzl, dv, fma, vax, vay
  real*8 fx2, fx1, fy2, fy1
  real*8 uxi,uyi,uzi,gami
  real*8 fv(nvx,nvy)

  if (nt.lt.itav) then
     icycav=nt
  else
     icycav=ncyc*itav
  endif


  if (itime.eq.0) then
     !  zero distn functions
     do i=1,nvx
        do j=1,nvx
           fv(i,j)=0.
        end do
     end do

     !  find max velocity
     umax=0.
     do l=1,ne

        !  boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)
        uzi=uz(l)

        vxl=uxi/gami
        vyl=uyi/gami
        vzl=uzi/gami
        umax=max(umax,vxl,vyl,vzl)
        dv=2*umax/nvx
     end do
  endif

  if (mod(itime,igr).eq.0) then

     if (itime.gt.0) then
        !  plot distn fn and re-zero;

        !  find max
        fma=0.
        do i=1,nvx
           do j=1,nvx
              fv(i,j)=fv(i,j)/icycav
              fma=max(fma,fv(i,j))
           end do
        end do
        open(95,file="fv.z")
        write(95,*) '! nx ',nvx/2,' ny ',nvx/2 &
             ,' xmax=',umax,' xmin= ',0.,' ymax= ',umax,' ymin= ',0.
        write(6,*) 'f(v) contour max',fma
        write(95,101) ((fv(i,j),i=nvx/2+1,nvx),j=nvx/2+1,nvx)
101     format(100f12.5)
        close(95)
     endif

     !  zero distn functions
     do i=1,nvx
        do j=1,nvx
           fv(i,j)=0.
        end do
     end do

     !  find max velocity
     umax=0.
     do l=1,ne

        !  boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)
        uzi=uz(l)

        vxl=uxi/gami
        vyl=uyi/gami
        vzl=uzi/gami
        umax=max(umax,vxl,vyl,vzl)
        dv=2*umax/nvx
     end do

  else if (mod(itime,igr).ge.igr-icycav) then
     !  running average of distn
     do l=1,ne

        !  boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)

        vxl=uxi/gami
        vyl=uyi/gami
        vax=abs(vxl+umax)/dv
        vay=abs(vyl+umax)/dv
        i1=vax+1
        i2=i1+1
        fx2=vax-(i1-1)
        fx1=1.-fx2
        j1=vay+1
        j2=j1+1
        fy2=vay-(j1-1)
        fy1=1.-fy2

        if (i1.le.nvx.and.j1.le.nvx) fv(i1,j1)=fv(i1,j1)+fx1*fy1
        if (i1.le.nvx.and.j2.le.nvx) fv(i1,j2)=fv(i1,j2)+fx1*fy2
        if (i2.le.nvx.and.j1.le.nvx) fv(i2,j1)=fv(i2,j1)+fx2*fy1
        if (i2.le.nvx.and.j2.le.nvx) fv(i2,j2)=fv(i2,j2)+fx2*fy2

     end do
  endif


end subroutine f2v
