#! /bin/sh 
#
#  Script to produce 2D plot of nonlinear surface oscillations
# 
# Run directory

RUN=nl_surf

#put here the same name as the .id file that you want to use 
#from the tools folder
SIM_TYPE=?????


# Top directory
BOPS=$HOME/bops/src/bops.exe
ODPP=$HOME/bops/tools/gle/odpp.sh
TOOLS_ID="../../tools/id"

# Check & cleanup 
if [ -d $RUN ] 
then
  echo "Run directory" $RUN "already exists"
  rm $RUN/movie/*
else
  echo "Creating run directory" $RUN
  mkdir $RUN
  mkdir $RUN/movie
fi
cd $RUN

# Inputs get copied into bops.indata
cat <<'EOF'>bops.indata
 &picohd
  nx=2000
  ne=40000
  ni=0
  trun = 80.

  a0=2.e19
  tpulse=100.
  theta0=30.
  cpolzn='P'
  miome=1836.
  nonc=15
  Te=1.0
  Ti=0.1

  xl=40.
  inprof=3
  xsol=11.
  xm1=30.
  xlolam=0.02
  xsol2=17.
  xm2=8.
  xcur1 = 25.
  xcur2 = 35.

  lrstrt=.false.
  ldump=.false.

  isubi=5
  ioboost=0
  iout=10
  igr=80
  igmovie=2
  igx2d = 1
  itc=10
  ncyc=4

  Z=1.
  amass=1.
  ilas=1
  w0=1.0
  wp=1.4142
  trise=6.
  tfall=30.
  vxm=0.2
  vym=0.5

  ipbc=4
  ifbc=2
  nsp=-1
  isp=5600

  ntrack=0
  itropt=1
  uhot=0.1
  xpint=0.01
  xpstart=10.0
  itstart=0
  itend=6000

  itsk=30
  igxs=2
  ipskip=9
  nuav=50
  nftem = 16384
  ift = 10
  omegm = 20.
  ifbin=5/

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

echo 'Running BOPS ..'
$BOPS
#$ODPP ${TOOLS_ID} $RUN ${SIM_TYPE} ? ?
#nb: vedere altri .sh per esempi sugli ultimi due parametri
echo 'Finished run'
