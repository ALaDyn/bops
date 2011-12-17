

!  ===========================
!
!   Graphical output of time-ave quantities
!
!  ===========================


subroutine grav
  use bopsvars
  implicit none

  integer i, icycdel
  real*8 xinc
  real*8, dimension(0:nx+1) :: work1,work2, work3
  real*8 :: emdensity(0:nx+1), gradw(0:nx+1),grad2w(0:nx+1)
  character chx*15

  write (*,*) 'Writing out time-averages'
  if (iunits.eq.2) then
     chx = '     x/mu      '
  else if (iunits.eq.1) then
     chx = '     k_px      '
  else
     chx = '     k_0x      '
  endif


  !  Lab frame oscillatory Ex
  do i=1,nx+1
     work1(i)=sqrt(abs(avex(i) - dcex(i)**2))
  end do

  call grxy(xx,work1,nx+1,20000+idc,igxs,1 &
       ,chx,'     Erms      ','exrm'//ctime(1:12) )

  !  Lab frame dc Ex
  call grxy(xx,dcex,nx+1,20500+idc,igxs,1 &
       ,chx,'     <Ex>      ','exdc'//ctime(1:12) )


  !  Lab frame oscillatory potential
  do i=1,nx+1
     work1(i)=sqrt(abs(avphi(i) - dcphi(i)**2))
  end do

  call grxy(xx,work1,nx+1,23000+idc,igxs,1 &
       ,chx,'   <Phi>rms    ','phrm'//ctime(1:12) )
  call grxy(xx,dcphi,nx+1,23500+idc,igxs,1 &
       ,chx,'     <Phi>     ','phdc'//ctime(1:12) )


  !  TE fields
!#### Anupam & Bin 2009/2010 : output all for either TE or TM mode
  !      if (ppol.ne.0) then

  !  Lab frame Ey
  do i=1,nx+1
     work1(i)=avey(i)**0.5
  end do

  call grxy(xx,work1,nx+1,21500+idc,igxs,1 &
       ,chx,'     <Ey>rms   ','eyrm'//ctime(1:12) )

  !  Lab frame oscillatory Bz
  do i=1,nx+1
     work1(i)=sqrt(abs(avbz(i) - dcbz(i)**2))
  end do

  call grxy(xx,work1,nx+1,22000+idc,igxs,1 &
       ,chx,'     <Bz>rms   ','bzrm'//ctime(1:12) )
  call grxy(xx,dcbz,nx+1,22500+idc,igxs,1 &
       ,chx,'     <Bz>      ','bzdc'//ctime(1:12) )

  !  AC current
  do i=1,nx+1
     work1(i)=sqrt(abs(avjy(i) - dcjy(i)**2))
  end do

  call grxy(xx,work1,nx+1,29000+idc,igxs,1 &
       ,chx,'     <Jy>rms   ','jyrm'//ctime(1:12) )

  !  DC current
  call grxy(xx,dcjy,nx+1,28500+idc,igxs,1 &
       ,chx,'     <Jy>      ','jydc'//ctime(1:12) )
  !  DC current jx
  call grxy(xx,dcjx,nx+1,33500+idc,igxs,1 &
       ,chx,'     <Jy>      ','jydc'//ctime(1:12) )

  !  DC charge density
  call grxy(xx,dcrhe,nx+1,29500+idc,igxs,1 &
       ,chx,'     <rhoe>    ','redc'//ctime(1:12) )


  !  pond force  vy x Bz; vy from py=ay+p0
  call grxy(xx,vxb,nx+1,26000+idc,igxs,1 &
       ,chx,'   <vy x Bz>   ','vxbp'//ctime(1:12) )

  !  pond force  vz x By; vz from pz=az
  call grxy(xx,vxb2,nx+1,27000+idc,igxs,1 &
       ,chx,'   <vz x By>   ','vxbs'//ctime(1:12) )

  !  pond force  <vy> x <Bz>
  do i=1,nx+1
     work1(i) = dcvy(i)*dcbz(i)
  end do
  call grxy(xx,work1,nx+1,26500+idc,igxs,1 &
       ,chx,'   <vy>x<Bz>   ','vxbd'//ctime(1:12) )

  !  pond force  Epond from pond. laser model

  call grxy(xx,avfp,nx+1,28000+idc,igxs,1 &
       ,chx,'     <fp>      ','fpav'//ctime(1:12) )

  !  E.J
  icycdel=(2*xm1/tav+0.5)
  xinc=uinc(ntav-icycdel)
  !     if(xinc.ne.0) then

  do i=1,nx+1
     work1(i)=edotj(i)
  end do

  call grxy(xx,work1,nx+1,25000+idc,igxs,1 &
       ,chx,'     <E.j>rms  ','edoj'//ctime(1:12))
  !     endif

  !      endif		!#### Anupam & Bin 2009/2010


  !  TM fields
!#### Anupam & Bin 2009/2010: output all for either TM or TE mode
 
!  if (spol.ne.0) then		!#### Anupam & Bin 2009/2010

     !  Lab frame Ez
     do i=1,nx+1
        work1(i)=avez(i)**0.5
     end do

     call grxy(xx,work1,nx+1,21000+idc,igxs,1 &
          ,chx,'     <Ez>rms   ','ezrm'//ctime(1:12) )

     !  Lab frame oscillatory By
     do i=1,nx+1
        work1(i)=sqrt(abs(avby(i) - dcby(i)**2))
     end do

     call grxy(xx,work1,nx+1,30000+idc,igxs,1 &
          ,chx,'     <By>rms   ','byrm'//ctime(1:12) )
     call grxy(xx,dcby,nx+1,30500+idc,igxs,1 &
          ,chx,'     <By>      ','bydc'//ctime(1:12) )

     !  Lab frame oscillatory Bx
     do i=1,nx+1
        work1(i)=sqrt(abs(avbx(i) - dcbx(i)**2))
     end do

     call grxy(xx,work1,nx+1,31000+idc,igxs,1 &
          ,chx,'     <Bx>rms   ','bxrm'//ctime(1:12) )
     call grxy(xx,dcbx,nx+1,31500+idc,igxs,1 &
          ,chx,'     <Bx>      ','bxdc'//ctime(1:12) )

     !  AC current jz
     do i=1,nx+1
        work1(i)=sqrt(abs(avjz(i) - dcjz(i)**2))
     end do

     call grxy(xx,work1,nx+1,32000+idc,igxs,1 &
          ,chx,'     <Jz>rms   ','jzrm'//ctime(1:12) )

     !  DC current
     call grxy(xx,dcjy,nx+1,32500+idc,igxs,1 &
          ,chx,'     <Jy>      ','jzdc'//ctime(1:12) )


!  endif		!#### Anupam & Bin 2009/2010

! EM energy density - sum all field components (normal incidence)
  emdensity(1:nx+1) = avez(1:nx+1) + avey(1:nx+1) + avby(1:nx+1) + avbz(1:nx+1)

  do i=1,nx
    gradw(i) = (emdensity(i+1)-emdensity(i))/dx
  end do
  gradw(nx+1) = gradw(nx)
  call filter1(gradw,nx+1)

  do i=1,nx
    grad2w(i) = (gradw(i+1)-gradw(i))/dx
  end do
  call filter1(grad2w,nx+1)

  call grxy(xx,emdensity,nx+1,36000+idc,igxs,1 &
          ,chx,'       W       ','emde'//ctime(1:12) )
  call grxy(xx,gradW,nx+1,36500+idc,igxs,1 &
          ,chx,'   gradW       ','emgr'//ctime(1:12) )
  call grxy(xx,grad2W,nx+1,36600+idc,igxs,1 &
          ,chx,'   grad2W       ','emg2'//ctime(1:12) )

end subroutine grav
