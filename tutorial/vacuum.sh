#! /bin/sh
#
#  Vacuum test of boost fields 
#  Set cpolzn='S', 'P' or 'C' for s-/p-/c-polarized light respectively
#  Change theta0=0 .. 89 to alter incidence angle
#  Set iboost=0 for lab-frame, 1 for boost frame fields 

# Run directory
RUN=vacuum

# Location of executable from tutorial/ directory
BOPS=`pwd`/../src/bops
# Location of executable from top directory
#BOPS=`pwd`/src/bops

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
! Run time
  trun = 100.
! Grid points  
  nx = 1000
! # electrons  
  ne = 50
! # ions  
  ni = 0
!  iunits = 1 ! wp, kp
  iunits = 0 ! Default

! Laser amplitude (-ve) or intensity in W/cm**2 (+ve)
  a0=-1.0
! Laser wavelength in microns  
  xlambda = 0.8
! Pulse length in 1/w0  
  tpulse = 40.
! Pulse shape (see below) 
  ilas = 5
! Pulse delay  
  trise = 0.

! Polarization
  cpolzn = 'P'
! Angle of incidence
  theta0=45.
! Output fields in lab frame (1 for boost frame)
  ioboost=0
    
! Mass ratio ion/electron
  miome = 3000.
! Initial electron temperature  
  Te = 0.1
! Initial ion temperature  
  Ti = 0.0
! Plasma density ne/nc  
  nonc = 1.e-5
! Density profile  
  inprof = 2
! Grid length in c/w0  
  xl = 200.
! Start of LH plasma
  xm1 = 199
! Ramp end  
  xsol = 200
! RH plasma ramp start  
  xsol2=440.
! RH plasma edge  
  xm2=450.

! Frequency of printed output
  iout = 10
! Frequency of 1D plots  
  igr = 100
  igmovie = 100
! Frequency of time sequence stores  
  itc=5

! Particle boundary conditions
  ipbc=4
! field boundary conditions (transmit)  
  ifbc=2
! Resolution of spatialplots  
  igxs=5
! Resolution of phase space plots  
  ipskip=1
! Max energy for spectra  
  umevmax=20. 
  nuav=1/

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
			5 = sin-squared pulse (tpulse)
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
cd ..
./odpp_vacuum vacuum
echo 'Finished run'
echo
