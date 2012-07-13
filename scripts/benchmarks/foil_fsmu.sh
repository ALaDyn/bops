#! /bin/sh
#
#  Isolated foil demo
# 
#  fixed ions
#  'Experimental mode' iunits=2:
#   times, pulse length in fs
#   lengths in microns 
#

# Run directory
RUN=foil_fs

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
  trun=100
  nx=4000
  ne=50000
  ni=0
  iunits= 2

  a0 = 5.e18
  xlambda=0.8 
  tpulse=50.
  theta0=45.0
  cpolzn='P'
  bchirp=0.
  miome=6000.
  nonc= 20.
  fcrit=0.25
  Te=2.0
  Ti=0.1
  xl=10.
  inprof=7
  dfoil=.5
  xm1=5.5
  xsol=3.5
  xlolam=0.2

  isubi=1

  iout=10
  igr=50
  igmovie=50
  itc=40
  ncyc=1

  Z=1.
  amass=2.
  ilas=1
  trise=20.
  tfall=260.
  vxm=0.2
  vym=0.5

  ipbc=4
  ifbc=2
  nsp=0
  isp=5600

  umevmax=2.0
  itsk=30
  igxs=5
  ipskip=2
  nuav=10
  nftem= 120
  ift= 1
  omegm= 20.
  ifbin=3/ 
EOF
#
#
echo 'Running bops ..'
# llrun -p1 $BOPS
$BOPS

echo 'Finished run'
cd ..
./odpp_summary foil_fs
