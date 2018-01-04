subroutine setup_arrays

  use bopsvars	
  implicit none

  integer :: ntm,nxm,npm,nvm,ncxm,ncym,nffm,ntrkm  ! array li
  integer :: ncycm, nspm, io
  real :: tdelay
  ntm = 2*nt/itc+10
  nxm = nx+10
  npm = ne+ni_tot+100       !#### Anupam & Bin 2009: change ni with ni_tot including protons
  nvm = nvx+10
  nffm = nf+10
  nptrk = 50
  ncxm = nx+10
  ncym = max(itav+1,2)
  ntrkm = nt/itsk+10
  ncycm = 2*nt/(min(itc,itav))+1
  nspm = 4*itav+1  ! movie snapshots

  if (debug.gt.1) then
   do io=6,15,9
    write(io,*) 'Array sizes:'
    write(io,'(8(a20,i8/))') &
	'max nx:',nxm, &
	'max npart:',npm, &
	'max vel:',nvm, &
	'max frequency:',nffm, &
	'# cycles:',ncycm, &
	'# dt/cycle:',itav, &
	'sample frequency:',itc, &
        'time history:',ntm, &
	'movie snapshots:',nspm
    end do
  endif


!    parameter(pi=3.141592654,c=3.e8 &
!      ,ntm=1000,nxm=25002,npm=10000002,nvm=300,nptrk=npm/100 &
!      ,nfm=10024,ncxm=50,ncym=256,nffm=15000,ntrkm=400 &
!      )

     allocate (gamma(npm),ux(npm),uy(npm),uz(npm),xo(npm),xn(npm),ax(npm), q(npm), m(npm), species(npm) )

      allocate ( ff(0:nxm),fb(0:nxm),gf(0:nxm),gb(0:nxm) &
      ,rhoe(0:nxm),rhoi(0:nxm), rhop(0:nxm) &
      ,ex(0:nxm),ey(0:nxm),ez(0:nxm) &
      ,bx(0:nxm),by(0:nxm),bz(0:nxm) &
      ,jm(0:nxm),jp(0:nxm) &
      ,jyi(0:nxm),jye(0:nxm),jyim(0:nxm),jyem(0:nxm) &
      ,jzm(0:nxm),jzp(0:nxm) &
      ,jzi(0:nxm),jze(0:nxm),jzim(0:nxm),jzem(0:nxm) &
      ,rhot(0:nxm),rk2(nxm),sm(nxm) &
      ,phi(0:nxm),ay(0:nxm),az(0:nxm) &
      , epond(0:nxm) )

! Diagnostic arrays
      allocate (uux(0:npm),uuy(0:npm),uuz(0:npm),xnlab(0:npm), &
                 vx(npm),vy(npm),vz(npm),gamlab(npm) )

      allocate ( yrhok(nxm),yphik(nxm) )

      allocate ( xx(0:nxm),xk(0:nxm), xw(0:nffm), & 
	vxt(0:ntm),vyt(0:ntm),xv(0:nvm),xpt(0:ntm) &
	,exsurf(0:ntm), eysurf(0:ntm),bzsurf(0:ntm) &
	,idenmax(0:ntm), edenmax(0:ntm) )
 
! cycle-averaged histories
      allocate ( uesa(0:ncycm),utha(0:ncycm),uema(0:ncycm),urha(0:ncycm) &
      ,utota(0:ncycm),u1a(0:ncycm),usysa(0:ncycm) &
      ,uinc(-ncycm/2:ncycm),uback(0:ncycm),xi0(-ncycm/2:ncycm) &
      ,uthea(0:ncycm),uthia(0:ncycm),uthpa(0:ncycm),utma(0:ncycm),rhoim(0:ncycm),xcrne(0:ncycm) &
      ,xcrni(0:ncycm),utrans(0:ncycm),eta(0:ncycm),ulha(0:ncycm))

       allocate (fvx(0:nvm),fvy(0:nvm),fvz(0:nvm),fu(0:nvm),fion(0:nvm),fproton(0:nvm) &
      ,ffl(0:nffm),fb0(0:nffm),fex(0:nffm),fjy(0:nffm) &
      ,gfl(0:nffm),gb0(0:nffm) &
      ,utot(0:ntm),uth(0:ntm),ues(0:ntm),udr(0:ntm),esoth(0:ntm) &
      ,uthe(0:ntm),uthi(0:ntm) &
      ,ufb(0:ntm),uff(0:ntm),uem(0:ntm),usys(0:ntm) &
      ,erh(0:ntm),elh(0:ntm),eta1(0:ntm), rhoe_last(0:nxm) &
      ,avex(0:nxm),avey(0:nxm),avbz(0:nxm),edotj(0:nxm),edotjx(0:nxm) &
      ,avez(0:nxm),avby(0:nxm),avbx(0:nxm),avfp(0:nxm) &
      ,avni(0:nxm),dcbz(0:nxm),avext(0:nxm),avjm(0:nxm),vxb(0:nxm) &
      ,avphi(0:nxm),dcex(0:nxm),exlab(0:nxm),bzlab(0:nxm),dcphi(0:nxm) &
      ,dcvy(0:nxm),dcjy(0:nxm),dcrhe(0:nxm),vxb2(0:nxm),avjy(0:nxm) &
      ,avjz(0:nxm),dcez(0:nxm),dcbx(0:nxm),dcby(0:nxm),dcjz(0:nxm) &
      ,dcjx(0:nxm),dcey(0:nxm) ) !#### Anupam & Bin 2009/2010: added dcey


      allocate (uesc(0:10*npm),uesci(0:10*npm) &
      ,uinj(0:10*npm),uest(0:nptrk),phat(0:nptrk) )

      allocate ( iesci(0:10*npm), iesc(0:10*npm) )

      allocate ( xtrk(40,ntrkm),twork(0:ntrkm),uxtrk(40,ntrkm) &
      , uytrk(40,ntrkm),uztrk(40,ntrkm),axtrk(40,ntrkm),ytrk(40,0:ntrkm) &
      ,itrack(nptrk) )
!      allocate ( isp(5000) )
      eta(1)=0.
      uinc(-ncycm/2:0) = 0.

end subroutine setup_arrays
