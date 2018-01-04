#
#  RPA of nm foil
#
# Run directory
RUN=$HOME/bops/tutorial/rpa1

# Top directory
BOPS=$HOME/bops/src/bops

if [ -d $RUN] 
then
  echo "Run directory " $RUN " already exists!"
else
  echo "Creating run directory" $RUN
  mkdir $RUN
fi
cd $RUN
#
echo "Cleaning up.."
rm -f *.xy *.2D plots.tar plots.tar.gz *.t

# Inputs get copied into bops.indata

cat <<'EOF'>bops.indata


 &picohd
  trun=300
!  nt=1
  nx=96000
  ne=96000
  ni=16000
  iunits=2    ! input times in fs; lengths in microns
  em_scheme=1

  a0=8e19      ! laser intensity (W/cm**-2)
  xlambda=0.8  ! laser wavelength (microns)
  tpulse=27.33 ! pulse length (fs)
  theta0=0.    ! incidence angle
  cpolzn='C'   ! polarization
  ilas=2       ! gaussian 
  trise=3.3333
  tfall=3.3333
  tdel =150 ! delay time (fs)

  miome=22032.  ! ion/electron mass ratio for C12
  mpome=1836.   ! proton/electron mass ratio
  Z=6.  ! C^6+ ions
  amass=1.
  Te=0.0  ! Start cold
  Ti=0.0
  Tproton=0.0
  nonc=400.
  fcrit=0.0


  target_config= 1  ! single ion species
!  target_config=3  ! foil+proton layers on both rear and front side
!  target_config=5  ! foil+proton layer at front side

  inprof=7          ! 9: three layer target profile 
                    ! 10: multispecies
  dfoil=0.005       ! Foil width (microns)
  xm1=6.0           ! Foil position (microns)
  x_layer=0.0       ! width of proton layer
!  rho_layer=4.082  ! np/nc
  rho_layer=0.
!  xsol=3.008
  xl=30.0 ! Grid length (microns)
  xlolam=0.
  rhotrak=10.

  lrstrt=.false.
  ldump=.false.

  isubi=1
  ioboost=0
  iout=10
  igr=10
  igmovie=10
  itc=1

  ipbc=3
  ifbc=1

  umevmax=4.  ! max electron energy in particle spectra (MeV)
  uimax=40.
  upmax=40.
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

echo 'Finished run!'
###end
