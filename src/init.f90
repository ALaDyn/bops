!
!    INIT
!
!     Initialization routine

!     $Revision: 31 $
!
!     ========================================================

subroutine init
  use bopsvars
  implicit none

  real*8 tp, umax
  integer nwkstn, jcon, ntrsto

  include 'names.h'
  !
  !   default data
  lrstrt=.false.
  ldump=.false.
  nt=1
  nx=10
  ne=10
  ni=0
!#### Anupam & Bin 2009/2010:specie density intitialization to zero  
  np=0           ! necessary to initial np ==0
  n_pp1=0
  n_pp2=0        
  n_pp=0
!#### Anupam & Bin 2009/2010
  dt=0.1
  ipbc=1
  ifbc=2
  iembc=2
  ilas=1
  xl=1.
  z=1
  qome=-1.
  miome=1836.
  mpome=1836.
  amass=1.
  a0max=0.
  theta0=0.
  vte=0.001            !#### probably can't be set to zero!
  vti=0.
  Te=0.
  Ti=0.
  Tproton=0.            !#### Anupam 2010
  nonc=1.5
  rho_layer=0.		!NOTE: default rho_layer ==0 and x_layer == 0
  x_layer=0.		!NOTE: ===> n_pp=0 ===> np=0 
  xlolam=0.
  umax=100.
  umevmax = .5
  uimax = 0.1
  upmax = 0.1
  fcrit=1.
  vxm=0.2
  vym=0.1
  uxma=1.0
  uyma=1.0
  duin = 0.
  nvx=400
  nvy=400
  ifreeze=0

  xlambda=1.06
  sigma = 5.0
  wp=0.1
  w0=1.0
  w02=1.0
  tpulse=1.0
  trise=1.0
  tfall=0.
  a02 = 0.
  cpolzn='P'
  bchirp = 0.

  xm1=0.
  xm2=1.0
  xsol=xm2
  erhb=0.
  elhb=0.
  idc=0
  itc=1
  ntc=0
  idia=1
  igr=1
  iout=1
  ncyc_total = 1
  igmovie=0
  inprof=1
  ncyc=3
  ioboost = 0
  dfoil=1.
  ! filter ion and electron density
  ifilt=1
  ipskip=1
  nesc=0
  nesci=0
  nuav=1
  ndenav=1
  igx2d=1
  isgy=1
  jcon=0
  ncycle=1
  icft = 1
  tcstart=100000.
  ittrk=0
  itstart=-1
  itend=0
  nftem=500000
  ift = 1
  omegm=20.
  ifbin=2
  rhotrak=1.5
  ntrack=0
  ntrsto=0
  uhot=0.
  itropt=1
  itsk=1
  abo=0.5
  !      abo=0.
  asm=0.5
  !      asm=0.0
  ampl=0.
  yi=(0.,1.)
  imaxw=1
  iran0 = 0
  !  units
  iunits=0
  tconv=1.
  xconv=1.
  !  Time-integrated diagnostics
  Uemin=0.
  Uemout=0.
  Uemtr=0.
  open(10,file='bops.indata')
  read (10,NML=picohd)
end subroutine init
