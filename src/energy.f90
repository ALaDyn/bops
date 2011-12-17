
!  ======================
!
!    Energy diagnostics
!
!   $Revision: 26 $
!  ======================

subroutine energy
  use bopsvars
  implicit none

  real*8 ute,uti,vvx,vvy,temp1,g,vsex,vsey,vsix,vsiy,utes,utis
  real*8 dxlab, ese, xmelab, xmilab, emte, emtm
  real*8 bzl, bxl, byl, exl, eyl, ezl, r1
  integer i, ip
  logical :: debug_out=.false.
  complex yes

  uth(ntc)=0.

  !   potential energy (includes dc) - lab frame
  dxlab=dx*gam0
  ese=0.
  do i=1,nx+1
     ese = ese + dxlab*(ex(i)+vy0*bz(i))**2
  end do

  ues(ntc)=ese

  !  electron and ion thermal energies

  xmelab=me/gam0**2
  xmilab=mi/gam0**2
  ute=0.d0
  do ip=1,ne
     ute=ute+xmelab*(gam0*(gamma(ip)-vy0*uy(ip)) - 1.0)
  end do

  uthe(ntc)=ute
  uti=0.d0
  utis=0.d0

  if (ni.gt.0) then
     do ip=ne+1,ne+ni
        uti=uti+xmilab*(gam0*(gamma(ip)-vy0*uy(ip)) - 1.0)
     end do
  endif

  uthi(ntc)=uti
  uth(ntc)=ute + uti

  !  field energy - lab frame: includes es and em components
  emte=0.
  emtm = 0.
  do i=1,nx
     bzl = bz(i)+vy0*ex(i)
     bxl = -vy0*ez(i)
     byl = by(i)/gam0
     eyl = ey(i)/gam0
     exl = ex(i)+vy0*bz(i)
     ezl = ez(i)
     emte = emte + 0.5*( eyl**2 + exl**2 + bzl**2)*dxlab
     emtm = emtm + 0.5*( ezl**2 + bxl**2 + byl**2)*dxlab

  end do
  uem(ntc)=emte + emtm
  !
  !   energy lost to RH boundary
  erh(ntc)=erhb
  !   energy lost to LH boundary
  elh(ntc)=elhb
  !
  !   total
  usys(ntc)=uth(ntc)+uem(ntc)
  utot(ntc)=usys(ntc)+erhb+elhb
  !
  !   ratio field to thermal
  r1=0.
  if (uth(ntc).ne.0) r1=ues(ntc)/uth(ntc)
  esoth(ntc)=r1
  !
  !  em absorption
  !  poynting flux at left hand boundary
  eta1(ntc)=(ff(1)*ff(1)-fb(1)*fb(1))/a0/a0

! Debug energies
     if (debug_out) write (*,'(i6,3(a,1pe12.5),a,2(1pe12.3))') &
       itime,' thermal:',uth(ntc), &
       ' field:',uem(ntc),' rhb ', erh(ntc),' laser:',ff(1),gf(1)

end subroutine energy
