
! ===========================================
!
!  Energy distribution of plasma electrons
!    - Thot plot
!
!  $ Revision: $
! ===========================================

subroutine udist
  use bopsvars
  implicit none

  real*8 uxl,uyl,uzl,gl
  real*8 du, u2, u, fmin
  real*8 :: work1(nvx), work2(nvx)
  integer i, l, i0, iu
  character*15 cu

  !  energy interval in mc^2
  du=umevmax/.511/nvx

  if (itime.eq.0) then
     do i=1,nvx
        fu(i)=0.
     end do


  endif



  !  running average of distn
  if (mod(itime,igr).ge.igr-nuav) then
     do i=1,nvx
        work1(i)=0.
     end do

     !  relativistic thermal energy distn (lab frame)
     do l=1,ne
        !  lab frame k.e.
        call sim2lab(vy0,gam0,uy(l),gamma(l),uyl,gl)
        u2=ux(l)**2+uz(l)**2+uyl**2
        u=u+u2/(gl+1.)
        u=gam0*( gamma(l)-vy0*uy(l) ) - 1.0
        iu = u/du+1
        i0=min(max(1,iu),nvx)
        work1(i0)=work1(i0)+1
     end do

     fmin=work1(1)

     do i=1,nvx
        if (work1(i).ne.0) fmin=min(fmin,work1(i))
     end do

     if (fmin.le.0) fmin=1

     do i=1,nvx
        if (work1(i).ne.0) then
           fu(i)=fu(i)+work1(i)/nuav
        else
           fu(i)=fmin
        endif
     end do
  endif


  if (mod(itime,igr).eq.0) then
     if (itime.gt.0) then
        !  plot distn fn and re-zero; Te in keV
        do i=1,nvx
           if (umevmax.ge.1) then
              xv(i)=(i*du-du/2)*.511
              cu = '      U(MeV)   '
           else
              xv(i)=(i*du-du/2)*.511e+3
              cu = '      U(keV)   '

           endif
           work2(i) = fu(i)
        end do

        call grxy(xv,fu,nvx,6500+idc,1,2 &
             ,cu,'     fe(U)     ','fuep'//ctime(1:12) )
     endif
     do i=1,nvx
        fu(i)=0.
     end do


  endif
end subroutine udist





