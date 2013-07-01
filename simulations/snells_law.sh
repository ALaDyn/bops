#! /bin/sh
#
#  Snell's Law demo 
#

# Run directory

RUN=snell

#put here the same name as the .id file that you want to use 
#from the tools folder
SIM_TYPE=?????


# Top directory
BOPS=$HOME/bops/src/bops.exe
ODPP=$HOME/bops/tools/gle/odpp.sh
TOOLS_ID="../../tools/id"

if [ -d $RUN ] 
then
  echo "Run directory" $RUN "already exists"
else
  echo "Creating run directory" $RUN
  mkdir $RUN
fi
cd $RUN

# Inputs get copied into bops.indata

cat <<'EOF'>bops.indata


 &picohd
  trun = 500.
  nx = 500
  ne = 10000
  ni = 0
  iunits = 0

  a0=1.e14
  xlambda = 1.06
  tpulse = 500.
  ilas = 1
  trise = 10.
  theta0 = 45.
  cpolzn = 'S'

  miome = 3000.
  nonc = 0.3
  fcrit = 0.25
  Te = 0.5
  Ti = 0.1
  xl = 40.
  inprof = 3
  xsol = 4.5
  xm1 = 20.
  xlolam = 0.001
  xsol2=17.
  xm2=8.

  lrstrt=.false.
  ldump=.false.

  isubi = 5
  ioboost = 0
  iout = 50
  igr = 500
  igmovie = 1
  itc=20
  ncyc=2
  igx2d = 2
  igxs = 1
  
  Z=1.
  amass=1.
  w0=1.0
  wp=1.4142
  tfall=30.
  vxm=0.2
  vym=0.5

  ipbc=4
  ifbc=2
  nsp=0
  isp=5600

  ntrack=0
  itropt=1
  uhot=0.1
  xpint=0.01
  xpstart=10.0
  itstart=0
  itend=6000

  itsk=30
  ipskip=2
  nuav=10
  nftem = 1024
  ift = 2
  omegm = 20.
  ifbin=3/

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
echo 'Running bops ..'
$BOPS
#$ODPP ${TOOLS_ID} $RUN ${SIM_TYPE} ? ?
#nb: vedere altri .sh per esempi riguardo gli ultimi due parametri
echo 'Finished run'
