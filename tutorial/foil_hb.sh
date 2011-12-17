#! /bin/sh
#
#  Isolated foil 
# 
#  iunits=2:
#   times, pulse length in fs
#   lengths in microns 
#

# Run directory
RUN=foil_hb

# Top directory
BOPS=`pwd`/../src/bops

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
  trun=150 ! Run time in fs
  nx=10000 ! # grid points
  ne=50000 ! # electrons
  ni=50000 ! # ions
  iunits= 2

  a0 = 5.e19  ! intensity         .......................
  xlambda=0.8 ! wavelength
  tpulse=150. ! pulse length in fs.......................
  trise=5.    ! rise time for flat-top profile
  theta0=0.   ! angle of incidence
  cpolzn='C'  ! polarization
  miome=1836. ! mass ratio
  nonc= 36. ! frozen H            .......................
  Z=1.      !                      .......................
  amass=1. ! atomic mass          .......................
  Te=1.0
  Ti=0.01
  xl=20.
  target_config=1 ! Single species ions
  inprof=7  ! foil config
  dfoil=2.  ! thickness (microns)  ......................
  xm1=10.   ! foil position
  xsol=3.5  ! position of max. density
  xlolam=0.2 ! density scalelength

  iout=10 ! printed o/p frequency
  igr=50  ! snapshot frequency
  igmovie=50 ! 'movie' snapshots
  itc=40
  ncyc=1

  ilas=1 ! squared                 .....................

  ipbc=4 ! BC: particles reflective
  ifbc=2 ! BC: em wave transmissive

  umevmax=10.0 ! Max electron energy in spectra plots
  uimax=10.0 ! Max ion energy
  upmax=10.0 ! Max proton energy
  igxs=5 ! plot resolution fields
  ipskip=1 ! plot resolution particles
  nuav=10 ! time ave. for energy spectra
  nftem= 120 ! HHG control
  ift= 1
  omegm= 20.
  ifbin=3/ 
EOF
#
#
echo 'Running bops ..'
$BOPS

echo 'Finished run'
cd ..
