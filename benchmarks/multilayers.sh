#!/bin/bash -x
#MSUB -N pic1d
#MSUB -l walltime=1:00:00
#MSUB -l nodes=1:ppn=1

### Start  of jobscript

cd ..
if [ -d layers] 
then
  echo "Run directory layers already exists"
else
  echo "Creating run directory layers"
  mkdir layers
fi
cd layers
echo "Cleaning up.."
rm -f *.xy *.2D plots.tar plots.tar.gz *.t

# Inputs get copied into bops.indata

cat <<'EOF'>bops.indata


 &picohd
  trun=500
!  nt=1
  nx=150000
  ne=9600
  ni=1600
  iunits= 2    ! input times in fs; lengths in microns
  em_scheme=1

  a0=3.45e19     ! laser intensity (W/cm�)
  xlambda=1.0  ! laser wavelength (microns)
  tpulse=126.6667 ! pulse length 
  theta0=0.    ! incidence angle
  cpolzn='C'   ! polarization
  ilas=4       ! sin� 
  trise=3.3333
  tfall=3.3333
  tdel =3.3333

  miome=22032.  ! ion/electron mass ratio
  mpome=1836.  ! proton/electron mass ratio
  Z=6.
  amass=1.
  Te=0.0
  Ti=0.0
  Tproton=0.0
  nonc=195.918
  fcrit=0.0

  target_config=3  ! foil+proton layers on both rear and front side
!  target_config=4  ! foil+proton layer at rear side
!  target_config=5  ! foil+proton layer at front side
  inprof=9         ! 9: three layer target profile 
                    ! 10: multispecies
  dfoil=0.008      ! Foil width
  xm1=3.0       ! Foil position
  x_layer=0.008    ! width of proton layer
  rho_layer=8.164  ! np/nc
  xsol=3.008
  xl=15.0 ! Grid length
  xlolam=0.
  rhotrak=10.

  lrstrt=.false.
  ldump=.false.

  isubi=1
  ioboost=0
  iout=1
  igr=1
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
			5 = sin�
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
../src/bops

echo 'Finished run!'
###end
