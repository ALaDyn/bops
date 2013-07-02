#!/bin/sh 

RUN=multirobinson

#put here the same name as the .id file that you want to use 
#from the tools folder
SIM_TYPE=???


BOPS=$HOME/bops/src/bops.exe
ODPP=$HOME/bops/tools/gle/odpp.sh
TOOLS_ID="../../tools/id"


if [ -d $RUN ] 
then
  echo "Run directory $RUN already exists"
else
  echo "Creating run directory $RUN"
  mkdir $RUN 
fi
cd $RUN 
echo "Cleaning up.."
rm -f *.xy *.2D plots.tar plots.tar.gz *.t

# Inputs get copied into bops.indata

cat <<'EOF'>bops.indata


 &picohd
!  trun=500
  trun=50
!  nt=1
  nx=10000
  ne=96000
  ni=16000
  iunits= 2    ! input times in fs; lengths in microns
  em_scheme=1

  a0=2e21     ! laser intensity (W/cm²)
  xlambda=1.0  ! laser wavelength (microns)
  tpulse=33. ! pulse length 
  theta0=0.    ! incidence angle
  cpolzn='C'   ! polarization
  ilas=4       ! sin² 
  trise=1.
  tfall=3.3333
  tdel =3.3333

  miome=22032.  ! ion/electron mass ratio
  mpome=1836.  ! proton/electron mass ratio
  Z=6.
  nonc=126.  ! Electron density to match heavy ion layer
  Te=5.0
  Ti=0.0
  Tproton=0.0
  fcrit=0.0

!  target_config=3  ! foil+proton layers on both rear and front side
  target_config=4  ! foil+proton layer at rear side
!  target_config=5  ! foil+proton layer at front side
  inprof=10         ! 9: three layer target profile 
                    ! 10: multispecies
  dfoil=20.      ! Foil width
  xm1=10.0       ! Foil position
  x_layer=0.0    ! width of proton layer
  rho_layer=42.  ! Proton density np/nc
  xsol=3.008
  xl=50.0 ! Grid length
  xlolam=0.
  rhotrak=10.

  lrstrt=.false.
  ldump=.false.

  isubi=1
  ioboost=0
  iout=5
  igr=10
  igmovie=1
  itc=1
  ncyc=1

  vxm=1.0
  vym=1.0

  ipbc=3
  ifbc=1
  nsp=0
  isp=5600

  ntrack=0
  itropt=1
  uhot=0.1
  xpint=0.01
  xpstart=10.0
  itstart=0
  itend=6000

  umevmax=800.  ! max energy in particle spectra
  uimax=800.
  upmax=800.
  itsk=30
  igxs=1
  ipskip=1
  nuav=10
  nftem= 1024
  ift= 2
  omegm= 20.
  ifbin=3 /
 &end
Glossary
========

Density profile
	       nonc   = n/nc
	       xlolam = L/lambda
	       inprof:  1 = uniform
			2 = linear from xm1 to xl
			3 = linear (xm1 to xsol) + flat top (xsol to
			xl)
			4 = linear + flat top + trailing ramp
			(xsol-xm2)
			5 = exp(-y/L), y=xsol-x + solid from xsol->xl
			6 = tanh(ay), y=x-xc, xc=xm1+(xsol-xm1)/2 + solid
			7 = foil target in centre of grid
			8 = layered target

Laser
	       a0     = vosc/c (lab frame)
	       ilas:    1 = uniform sinusoid (trise=rise-time)
			2 = gaussian
			3 = beat-wave
			4 = triangular (trise, tfall)
			5 = sin²
	       theta0 = obliquity (from normal)
               cpolzn = polarization:  'P'  (default), 'S' or 'C'

Plasma
	       ne     = no. electrons
	       ni     = no. ions (0 for fixed ions)
	       miome  = mass ratio
	       Z      = ion charge
	       amass  = ion atomic weight (multiplies miome)
	       Te     = electron temp (keV)
	       Ti     = ion temp (keV)

Diagnostics
	       igr    = graphical snapshots
	       itc    = store for history plots
	       iout   = printed output
	       igmovie = time-sequence snapshot freq.
               nav    = # cycles for time-average plots

EOF

#
#
######################################################

echo 'Running bops ...'
$BOPS
#$ODPP ${TOOLS_ID} $RUN ${SIM_TYPE} ? ?
#nb: vedere altri .sh per esempi riguardo gli ultimi due parametri
echo 'Finished run!'

