module bopsvars
  real*8, parameter :: pi=3.141592654,c=3.0e8             ! constants

  real*8, allocatable :: gamma(:), ux(:), uy(:), uz(:)    ! velocities, gamma
  real*8, allocatable :: xo(:), xn(:)  ! particle positions
  real*8, allocatable :: ax(:)      ! particle acceleration
  real*8, allocatable :: q(:), m(:)      ! particle charge, mass
  integer, allocatable :: species(:)  ! particle species label 1=electron, 2=ion, 3=proton

  real*8, allocatable :: ff(:), fb(:), gf(:), gb(:)  ! forward & backward em fields
  real*8, allocatable :: rhoe(:)    ! electron density
  real*8, allocatable :: rhoi(:)    ! ion density
  real*8, allocatable :: rhop(:)    ! proton density
  real*8, allocatable :: rhot(:)    ! rhoe+rhoi+rhop
  real*8, allocatable :: ex(:), ey(:), ez(:), bx(:), by(:), bz(:), ay(:), az(:), epond(:), phi(:)   ! fields
  real*8, allocatable :: jm(:), jp(:), jyi(:), jye(:), jyim(:), jyem(:), jzi(:), jze(:), jzim(:), jzem(:), jzp(:), jzm(:)
  real*8, allocatable :: ffl(:), fb0(:), fex(:), fjy(:), gfl(:), gb0(:)   ! Field FTs
  real*8, allocatable :: xw(:), rk2(:)
  complex, allocatable :: yrhok(:), yphik(:)  ! FT arrays

  real*8, allocatable :: xx(:), xk(:), vxt(:), vyt(:), xpt(:), xv(:), sm(:)   ! grid arrays
  real*8, allocatable :: fvx(:), fvy(:), fvz(:), fu(:), fion(:), fproton(:)  ! distribution functions
  real*8, allocatable :: utot(:), uth(:), ues(:), udr(:), esoth(:), uthe(:), uthi(:), ulha(:), idenmax(:), edenmax(:)
  real*8, allocatable :: ufb(:), uff(:), eta(:), uem(:), usys(:), erh(:), elh(:), eta1(:)
  real*8, allocatable :: avex(:), avey(:), avbz(:), edotj(:), edotjx(:), avez(:), avby(:), avbx(:), avfp(:)
  real*8, allocatable :: avni(:), dcbz(:), avext(:), avjm(:), vxb(:), avphi(:), dcex(:), exlab(:), bylab(:)
  real*8, allocatable :: dcvy(:), dcjy(:), dcrhe(:), vxb2(:), avjy(:), avjz(:), dcey(:), dcbx(:), bzlab(:), rhoe_last(:)
  real*8, allocatable :: dcez(:), dcby(:), dcjz(:), dcjx(:), dcphi(:)
  real*8, allocatable :: xtrk(:,:), twork(:), itrack(:), uxtrk(:,:), uytrk(:,:), uztrk(:,:), axtrk(:,:), ytrk(:,:) ! tracking arrays

  real*8, allocatable :: uinj(:), uest(:), phat(:), utrans(:)

  real*8, allocatable :: uesa(:), utha(:), uema(:), urha(:), uesc(:), uesci(:)        ! time-averaged
  real*8, allocatable :: utota(:), u1a(:), usysa(:), uinc(:), uback(:), xi0(:), uthea(:), uthia(:), uthpa(:), utma(:)
  real*8, allocatable :: exsurf(:), eysurf(:), ezsurf(:), bzsurf(:), rhoim(:), xcrne(:), xcrni(:)

  real*8, allocatable :: uux(:),uuy(:),uuz(:),xnlab(:)
  real*8, allocatable :: vy(:),vz(:),vx(:),gamlab(:) !#### Anupam & Bin 2009/2010:added vx

  integer, allocatable :: iesci(:), iesc(:)     ! bookkeeping for particles
  integer :: isp(10000)
  logical :: lun,lshif,lrnorm,lmovie,lrstrt,ldump,lpow2

! grid stuff
  integer :: nx,   nt, nxo2, itime,   nvx,   nvy,   nf,  isubi, isube &
     , nfo2, itim0

  real*8 :: xl,    dx,   xn0,    dkx,  dt,   dvx,   dvy,   vxm &
      ,   vym,  dto2,   dw,   dti,  uxma,  uyma &
      ,   dte,  ttot, trun,  tconv,  xconv, tcfs, xmu

  complex :: yi      ! sqrt(-1)


!   plasma parameters
  integer :: target_config=0  ! Target type switch 
				! 0=fixed ions
				! 1=single species ions
				! >2= multi-layer - see manual
  integer :: inprof  ! Density profile type
  real*8 :: qe    ! electron macro-charge
  real*8 :: qi    ! ion macro-charge
  real*8 :: qp    ! proton macro-charge
  real*8 :: me    ! electron macro-mass
  real*8 :: mi    ! ion macro-mass
  real*8 :: mp    ! proton macro-mass
  real*8 :: qome  ! electron charge/mass ratio (= -1)
  real*8 :: qomi  ! ion charge/mass ratio 
  real*8 :: qomp  ! proton charge/mass ratio 
  real*8 :: miome=3600. ! ion/electron mass ratio #### Anupam & Bin 2009/2010:given by input
  real*8 :: mpome=1836. ! proton/electron mass ratio #### Anupam & Bin 2009/2010:given by input
  real*8 :: z     ! ionization state of ions/ charge state of ions (species==2)
  real*8 :: vte   ! electron thermal velocity
  real*8 :: vti   ! ion thermal velocity
  real*8 :: vtp   ! proton thermal velocity #### Anupam & Bin 2009/2010
  real*8 :: xpert,   upe, vpert, xload, xllab &
      ,     ncell, rhoex, Qgrid &
      ,     Te,    Ti, Tproton, amass, dfoil,  gap_foil, xcrit, rhotrak !#### Anupam & Bin 2010: added Tproton

  real*8 :: umevmax ! Max electron energy for spectral plots
  real*8 :: uimax   ! Max ion energy
  real*8 :: upmin   ! Max proton energy towards laser
  real*8 :: upmax   ! Max proton energy

  real*8 :: theta0, the0, gamma0, uy0, vy0,  gam0        ! boost veloci

  integer ::  ne 	 ! # electrons
  integer ::  ni_tot  	 ! total # ions including protons
  integer ::  np     !=0   ! # protons #### Anupam & Bin 2009/2010: should be calculated in the code
  integer ::  ni     ! # heavy ions
  integer ::  npart  ! total # particles
  integer ::  nesc   ! # electrons escaped from RHB
  integer ::  nesci  ! # ions escaped from RHB
  integer ::  ion1   ! start # of 1st ion

  integer ::  mode
  integer :: nfoil=1  ! # foils for config inprof=17 

!   laser parameters
  real*8 ::  a0, a0max,   w0, trise,  tdel, tpulse, delfac &
      , tfall,  ppol,  spol,xlambda, sigma, theta_opt &
      , a02, trise2, tfall2, tpulse2, tsep,  w02, bchirp,w0r

  integer :: ilas  	     ! Laser pulse shape
  integer :: em_scheme=1 ! Scheme for EM fields 
  integer :: push_scheme=1 ! Scheme for particle pusher 

  real*8 :: rho0  ! Max normalized initial plasma density in code frame 
  real*8 :: rho0lab ! Max density in lab frame (n0/nc)
  real*8 :: rho_layer   ! Density of additional proton layer (/nc)
  real*8 :: x_layer  ! width of proton layer

  real*8 ::    wp,   v01,   v02,   xm1,   xm2,  xsol &
      ,   abo,   asm,  ampl, xsol2,   tav,  vyav,  gnon,  xgin &
      ,  erhb,  elhb,    wl,    wr,  rsv2,  thp, uthe0, uthi0, uthp0 &
      ,xdebye, xdodx, xlolam,  nonc, uemin, uemout, uemtr, omegm &
      , wrf, denmin, Upoy_in, Upoy_out

  integer :: ifbin


!    switches
  integer ::  ipbc,  ifbc, imaxw, ifilt &
      , iembc, ioboost, iunits, iran0, ifreeze
  integer :: debug=2 ! Debug level switch
  logical :: lcycave=.false. ! Write out time-averaged energies
!   i/o
  integer ::  iout,  igr, idia,  ipskip,   itc,  itav,   nuav &
      ,  igx2d,  isgy,   igxs,  nsp, igmovie, icstart,  icft,ncycle &
      , nftem, ndenav, ncyc, ift, idc,   ntc,  ntav, &
	ncyc_total, & ! total # laser cycles in simulation 
        ntrack, ittrk, itstart, itend, itsk, itropt, nptrk
  integer :: iplot2d=2000  ! frequency of 2D plots
  real :: tgr=-1
  character :: cw*40, ctime*16, cpolzn*1


  real*8 ::  xcur1,  xcur2, tcstart, fcrit, duin, xpint, xpstart, uhot, ttrans

! #### Anupam & Bin 2009/2010
! n_pp1: the proton particles for the first layer ; n_pp2 : the proton particles for the second layer

  integer :: n_pp, n_pp1,n_pp2  !for double/three layer target the number of proton particles (n_pp) and ion numbers (ni-n_pp)

end module bopsvars





