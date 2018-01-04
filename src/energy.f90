
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
  real*8 :: gl, u2, xmplab, utp, uyl
  integer i, ip, l
  logical :: debug_out=.false.
  complex yes

  uth(ntc)=0.

  !   potential energy (includes dc) - lab frame
  dxlab=dx*gam0
  ese=0.
  do i=1,nx+1
     ese = ese + 0.5*dxlab*(ex(i)+vy0*bz(i))**2
  end do

  ues(ntc)=ese

  !  electron and ion thermal energies

  if (ioboost==0) then
   xmelab=me/gam0**2
   xmilab=mi/gam0**2
  else
   xmelab=me
   xmilab=mi
  endif

  ute=0.d0
  do ip=1,ne
	if (ioboost==0) then
         call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyl,gl)
	else
	 uyl=uy(ip)
        endif 
     u2=ux(ip)**2+uz(ip)**2+uyl**2
     ute=ute+xmelab*u2/(gl+1.d0)
  end do

  uthe(ntc)=ute
  uti=0.d0
  utis=0.d0
  utp=0.d0

     !  ion thermal energies - lab frame
     uti=0.d0
     utp=0.d0

  if (ni_tot.gt.0) then
     do l=1,ni_tot
        ip=ne+l
	if (ioboost==0) then
         call sim2lab(vy0,gam0,uy(ip),gamma(ip),uyl,gl)
	else
	 uyl=uy(ip)
        endif        
	u2=ux(ip)**2+uz(ip)**2+uyl**2
 	if (species(ip).eq.2) then
          uti=uti+xmilab*u2/(gl+1.d0)  ! heavy ions
  	else if (species(ip).eq.3) then
          utp=utp+xmplab*u2/(gl+1.d0)  ! protons
	endif
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
    if (ioboost.eq.0) then
     emte = emte + 0.5*( eyl**2 + exl**2 + bzl**2)*dxlab
     emtm = emtm + 0.5*( ezl**2 + bxl**2 + byl**2)*dxlab
    else
     emte = emte + 0.5*( ex(i)**2 + ey(i)**2 + bz(i)**2)*dx
     emtm = emtm + 0.5*( bx(i)**2 + by(i)**2 + ez(i)**2)*dx
    endif
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
