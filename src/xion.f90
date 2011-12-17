
!  ========================================================
!
!       Initialise ion positions
!        -  placed on top of electrons
!
!  ========================================================

subroutine xion
  use bopsvars
  implicit none

  integer iz, l, ipi
  if (ni.eq.0) return
  iz=Z
  l=1
  !  place ions at every Zth electron
  do ipi=ne+1,npart
     xo(ipi)=xo(l)
     l=l+iz
  end do

end subroutine xion
