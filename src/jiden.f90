
!  ==================================
!
!    Ion current density
!
!  ==================================

subroutine jiden(i1)
  use bopsvars
  implicit none

  integer i,l, ip, i1, io1, io2, in1, in2, iwl, iwr
  real*8 ri, xan, xao, fm1, fm2, fp1, fp2, rdx, jy, jz

  rdx = 1./dx

  if (ni_tot.eq.0 ) then
     !   add neutralising current
     do i=1,nx+1
        jyim(i) = z*rhoi(i)*vy0
        jyi(i) = z*rhoi(i)*vy0
        jzim(i) = 0.
        jzi(i) = 0.
     end do
     return
  else
     do i=1,nx+1
        jyim(i)=0.
        jyi(i)=0.
        jzim(i)=0.
        jzi(i)=0.
     end do
  endif


  do l=1,ni_tot
     ip = i1+l-1
     xan = xn(ip)*rdx
     xao = xo(ip)*rdx
     ri=q(ip)/dx
     jy = ri*uy(ip)/gamma(ip)
     jz = ri*uz(ip)/gamma(ip)
     io1 = xao
     io2 = io1+1
     in1 = xan
     in2 = in1+1
     fm2 = xao-io1
     fm1 = 1.-fm2
     fp2 = xan-in1
     fp1 = 1.-fp2
     !  TE mode
     !   j-
     jyim(io1+1) = jyim(io1+1) + jy*fm1
     jyim(io2+1) = jyim(io2+1) + jy*fm2
     !   j+
     jyi(in1+1) = jyi(in1+1) + jy*fp1
     jyi(in2+1) = jyi(in2+1) + jy*fp2

     !  TM mode
     !   j-
     jzim(io1+1) = jzim(io1+1) + jz*fm1
     jzim(io2+1) = jzim(io2+1) + jz*fm2
     !   j+
     jzi(in1+1) = jzi(in1+1) + jz*fp1
     jzi(in2+1) = jzi(in2+1) + jz*fp2
  end do

  !  boundaries: fold ion currents
  iwl=wl/dx+1
  iwr = wr/dx+1

  !  TE mode

  !   j-
  jyim(iwl+1) = jyim(iwl+1) + jyim(iwl)
  jyim(iwl)=0.
  !   j+
  jyi(iwl+1) = jyi(iwl+1) + jyi(iwl)
  jyi(iwl) =0.

  !   j-
  jyim(iwr) = jyim(iwr) + jyim(iwr+1)
  jyim(iwr+1) =0.
  !   j+
  jyi(iwr) = jyi(iwr) + jyi(iwr+1)
  jyi(iwr+1)=0.

  !  TM mode

  !   j-
  jzim(iwl+1) = jzim(iwl+1) + jzim(iwl)
  jzim(iwl)=0.
  !   j+
  jzi(iwl+1) = jzi(iwl+1) + jzi(iwl)
  jzi(iwl) =0.

  !   j-
  jzim(iwr) = jzim(iwr) + jzim(iwr+1)
  jzim(iwr+1) =0.
  !   j+
  jzi(iwr) = jzi(iwr) + jzi(iwr+1)
  jzi(iwr+1)=0.

end subroutine jiden
