
!  ==================

subroutine vdist
  use bopsvars
  implicit none

  integer icycav, i,l, i0, j0, k0, isp1
  real*8 umax, vxl, vyl, vzl, dv, vax, vay, vaz, ampp, xlol
  real*8 uyi,uxi,uzi,gami
  real*8 :: work1(0:nvx), work2(0:nvx), work3(0:nvx)

  if (nt.lt.itav) then
     icycav=nt
  else
     icycav=ncyc*itav
  endif

  if (itime.eq.0) then
     !  zero distn functions
     do i=1,nvx
        fvx(i)=0.
        fvy(i)=0.
        fvz(i)=0.
     end do

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

  if (mod(itime,igr).ge.igr-icycav) then
     !  running average of distn
     do l=1,ne

        !  boost back to lab
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyi,gami)
        uxi=ux(l)
        uzi=uz(l)

        vxl=uxi/gami
        vyl=uyi/gami
        vzl=uzi/gami
        vax=abs(vxl+umax)/dv
        vay=abs(vyl+umax)/dv
        vaz=abs(vzl+umax)/dv
        i0=vax+1
        j0=vay+1
        k0=vaz+1
        i0=max(1,min(i0,nvx))
        j0=max(1,min(j0,nvx))
        k0=max(1,min(k0,nvx))
        fvx(i0)=fvx(i0)+1
        fvy(j0)=fvy(j0)+1
        fvz(k0)=fvz(k0)+1
     end do
  endif


  if (mod(itime,igr).eq.0) then

     if (itime.gt.0) then
        !  plot distn fn and re-zero;
        do i=1,nvx
           xv(i)=i*dv-dv/2.-umax
           work1(i)=fvx(i)/icycav
           work2(i)=fvy(i)/icycav
           work3(i)=fvz(i)/icycav
        end do

        call grxy(xv,work1,nvx,6000+idc,1,1 &
             ,'      vx       ','     f(vx)     ','fvxe'//ctime(1:12) )
        call grxy(xv,work2,nvy,6200+idc,1,1 &
             ,'      vy       ','     f(vy)     ','fvye'//ctime(1:12) )
        call grxy(xv,work3,nvy,6300+idc,1,1 &
             ,'      vz       ','     f(vz)     ','fvze'//ctime(1:12) )
     endif

     !  zero distn functions
     do i=1,nvx
        fvx(i)=0.
        fvy(i)=0.
        fvz(i)=0.
     end do

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
end subroutine vdist























