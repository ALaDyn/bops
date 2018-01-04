#! /bin/sh
#
#  Foil demo 
#

# Run directory
RUN=foil1

# Top directory
BOPS=`pwd`/src/bops

if [ -d $RUN ] 
then
  echo "Run directory" $RUN "already exists"
else
  echo "Creating run directory" $RUN
  mkdir $RUN
fi
cd $RUN
echo "Cleaning up.."
rm -f *.xy *.2D plots.tar plots.tar.gz *.t

# Inputs get copied into bops.indata

cat <<'EOF'>bops.indata


 &picohd
  trun=600
  nx=2000
  ne=50000
  ni=50000
  iunits= 1

  a0=2.e18
  xlambda=0.8 
  tpulse=90.
  theta0=0.
  cpolzn='P'
  miome=1836.
  nonc=10.
  fcrit=0.25
  Te=1.0
  Ti=0.1
  xl=150.
  inprof=57
  dfoil=30.
  xm1=20.
  xsol=90
  xlolam=1.
  rhotrak=2.

  lrstrt=.false.
  ldump=.false.

  isubi=2
  ioboost=0
  iout=10
  igr=200
  igmovie=300
  itc=40
  ncyc=2

  Z=1.
  amass=1.
  ilas=1
  trise=30.
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

  umevmax=10.
  itsk=30
  igxs=4
  ipskip=2
  nuav=10
  nftem= 300
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
# llrun -p1 $BOPS
$BOPS

echo 'Finished run'
#echo 'Packing plots..'
#rm -f plots.tar plots.tar.gz
#tar cvf plots.tar *.xy bops.*
#gzip plots.tar
