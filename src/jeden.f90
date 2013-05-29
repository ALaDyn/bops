
!  ==================================
!
!    Electron current density
!
!    23/8/95  -  upgrade for p/s-pol
!
!  ==================================

subroutine jeden
  use bopsvars
  implicit none

  real*8 :: re, xan, xao, fm1, fm2, fp1, fp2, rdx, jy, jz
  real*8 :: jye_pr(0:nx+1), jyem_pr(0:nx+1), jze_pr(0:nx+1), jzem_pr(0:nx+1)
  integer :: i, ip, io1, io2, in1, in2, iwl, iwr

  rdx = 1./dx

  !$omp parallel
  !$omp& private(jye_pr,jze_pr,jyem_pr,jzem_pr)
  !$omp& private(io1,io2,in1,in2,fm2,fm1,fp1,fp2,xan,jy,jz)
  !$omp& default(shared)
  do i=1,nx+1
     jye(i)=0.
     jyem(i)=0.
     jze(i)=0.
     jzem(i)=0.
     jye_pr(i)=0.
     jyem_pr(i)=0.
     jze_pr(i)=0.
     jzem_pr(i)=0.
  end do

  !$omp do
  do ip=1,ne
     xan = xn(ip)*rdx
     xao = xo(ip)*rdx
     re=q(ip)/dx
     jy = re*uy(ip)/gamma(ip)
     jz = re*uz(ip)/gamma(ip)
     io1 = xao
     io2 = io1+1  ! io1 < xao < io2
     in1 = xan
     in2 = in1+1
     fm2=xao-io1
     fm1=1.-fm2
     fp2=xan-in1
     fp1=1.-fp2

     !  TE currents

     !   j-
     jyem_pr(io1+1) = jyem_pr(io1+1)+jy*fm1
     jyem_pr(io2+1) = jyem_pr(io2+1)+jy*fm2
     !   j+
     jye_pr(in1+1) = jye_pr(in1+1)+jy*fp1
     jye_pr(in2+1) = jye_pr(in2+1)+jy*fp2

     !  TM currents

     !   j-
     jzem_pr(io1+1) = jzem_pr(io1+1)+jz*fm1
     jzem_pr(io2+1) = jzem_pr(io2+1)+jz*fm2
     !   j+
     jze_pr(in1+1) = jze_pr(in1+1)+jz*fp1
     jze_pr(in2+1) = jze_pr(in2+1)+jz*fp2

  end do
  !$omp end do

  !  Now do OMP reduce on grid: thread-local copies of jye_pr, etc
  !  contain particle current density sums

  do i=0,nx+1
     !$omp atomic
     jye(i) = jye(i) + jye_pr(i)
     !$omp atomic
     jze(i) = jze(i) + jze_pr(i)
     !$omp atomic
     jyem(i) = jyem(i) + jyem_pr(i)
     !$omp atomic
     jzem(i) = jzem(i) + jzem_pr(i)
  end do

  !$omp end parallel


  !  boundaries: fold electron currents
  iwl=wl/dx+1
  iwr = wr/dx+1

  !  TE mode

  !   j-
  jyem(iwl+1) = jyem(iwl+1) + jyem(iwl)
  jyem(iwl)=0.
  !   j+
  jye(iwl+1) = jye(iwl+1) + jye(iwl)
  jye(iwl)=0.

  !   j-
  jyem(iwr) = jyem(iwr) + jyem(iwr+1)
  jyem(iwr+1)=0.
  !   j+
  jye(iwr) = jye(iwr) + jye(iwr+1)
  jye(iwr+1)=0.

  !  TM mode

  !   j-
  jzem(iwl+1) = jzem(iwl+1) + jzem(iwl)
  jzem(iwl)=0.
  !   j+
  jze(iwl+1) = jze(iwl+1) + jze(iwl)
  jze(iwl)=0.

  !   j-
  jzem(iwr) = jzem(iwr) + jzem(iwr+1)
  jzem(iwr+1)=0.
  !   j+
  jze(iwr) = jze(iwr) + jze(iwr+1)
  jze(iwr+1)=0.

end subroutine jeden
