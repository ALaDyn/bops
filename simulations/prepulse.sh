#!/bin/bash
#
#  Studio del prepulse: lunghissima onda piana sinusoidale di bassa intensita' 
#  a0 = sqrt( I*(lambda)^2 / (1.38*10^18) )
#  Pertanto se noi abbiamo un prepulso di 1 ns di intensita' 10^10 W/cm^2, 
#  con un laser avente lambda=0.8, questo corrisponde ad
#  un a0 = 0.0000681 = 6.81 * 10^(-5)


#  Nell'input file e' definito iunits=2; significa che tempi e lunghezza
#  dell'impulso sono in femtosecondi, le altre lunghezze in micrometri

#  Nome simulazione
RUN=prepulse

#put here the same name as the .id file that you want to use 
#from the tools folder
SIM_TYPE=foil


# Top directory
BOPS=$HOME/bops/src/bops.exe
ODPP=$HOME/bops/tools/gle/odpp.sh
TOOLS_ID="../../tools/id"
INPUT_FILE="bops.indata"


###############################
#    simulation parameters    #
###############################

TMAX_FEMTOSECONDS=3000
GRID_POINTS=50000
TOTAL_NUMBER_OF_ELECTRONS=1000000
TOTAL_NUMBER_OF_IONS=1000000
iunits=2

a0=6.81e-5
LASER_WAVELENGTH_MICROMETERS=0.8
LASER_PULSE_DURATION_FEMTOSECONDS=1000
LASER_ANGLE_OF_INCIDENCE=0.0
LASER_POLARIZATION="P"

#NB: LA MASSA DI UN ELETTRONE E' MEZZO MeV, QUELLA DI UN PROTONE 938, PERTANTO IL LORO RAPPORTO E' POCO MENO DI 2000
MASS_RATIO_ION_OVER_ELECTRON=1836.153

#NB: LA DENSITA' DELL'IDROGENO "SOLIDO" (CONGELATO) E' 36 nc. USIAMO IDROGENO CONGELATO VISTO CHE GLI IONI SONO PROTONI PER ORA
DENSITY_IN_nc_UNITS=36

#CARATTERISTICHE IONI
Z=1.0
A=1.0

#CONDIZIONI INIZIALI PLASMA
INITIAL_ELECTRON_TEMPERATURE_KEV=1.0
INITIAL_ION_TEMPERATURE_KEV=0.01


#BOUNDARY CONDITIONS FOR PARTICLES (ipbc)
# 1 periodic particles
# 2 reflective particles
# 3 absorb/reemit at both sides
# 4 absorb ions at LHB, electrons only if charged
# 5 absorb electrons, reflect ions at RHB
PARTICLE_BOUNDARYCONDITIONS=4

#BOUNDARY CONDITIONS FOR FIELDS (ifbc)
# 1 periodic fields
# 2 bounded fields { reflective at solid, RHB)
FIELD_BOUNDARYCONDITIONS=2



#target_config 0=fixed ions, 1=single ion species, 2=additional ion layers (protons)
TARGET_CONFIG=1

NUMBER_OF_OUTPUTS=100
SNAPSHOT_FREQUENCY=50
MOVIE_SNAPSHOTS=50
NUMBER_OF_HISTORY_PLOTS=50
NUMBER_OF_CYCLES_FOR_TIMEAVERAGED_PLOTS=1

#PLOT PARAMETERS
ELECTRON_EMAX_MEV=50.0
ION_EMAX_MEV=10.0
PROTON_EMAX_MEV=10.0
SUBSAMPLE_GRID_POINTS=100
SUBSAMPLE_PARTICLES=50
NUMBER_OF_TIMESTEPS_FOR_AVERAGED_DISTRIBUTION_FUNCTIONS=100
FREQUENCY_OF_FOURIER_TRANSFORM_PLOTS=120
SKIP_FACTOR_FT=1
FREQUENCY_MAX_FT_PLOTS=50.0
BINNING_FT_PLOTS="3/"


#Laser type:
#1 uniform sinusoid
#2 gaussian I(t) = I0 exp{-(t-td)^2/tp^2)}  (tp = t_fwhm/2*sqrt(log 2))
#3 beat-wave (2-frequency pump)
#4 triangular (t_rise, t_fall)
#5 I(t) = I0 sin^2(pi*t/tp)      (tp = 2*t_fwhm)
LASER_TYPE=1


#inprof
# 1   uniform profile
# 2   linear ramp from xm1 to xl
# 3   linear ramp (xm1 - xsol) + flat top (xsol - xl) scalelength xlolam
# 4   linear + flat top + trailing ramp (xsol - xm2)
# 5   exponential ramp xm1 to xsol, scalelength xlolam
# 6   tanh ramp (xm1 - xsol), scalelength xlolam
# 7   foil, thickness dfoil, starting at xm1
# 8   2 uniform layers with densities nlayer (xm1-xm2) and n0 (xsol-xsol2)
# 57  foil thickness dfoil, with exponential ramp starting at xm1, scalelength xlolam
INPROF=7

############################

#### INPROF = 4

#          ___________
#         /           \
#        /             \
#       /               \
# ____ /_________________\________
# ^    ^  ^          ^   ^       ^ 
# x0  x1  x2         x3  x4     x5

#x0 = 0 per definizione
#x1
#### x1 si chiama xm1 nell'input file
#x2
#### x2 si chiama xsol nell'input file
#x3
#### x3 si chiama xsol2 nell'input file
#x4
#### x4 si chiama xm2 nell'input file
#x5
#
#
#Va inoltre definito un parametro di scaling della densita'
#che nell'input file si chiama xlolam (ad esempio xlolam=0.2)

##############################

#### INPROF = 7

#       ___________
#      |           |
#      |           |
#      |           |
# ____ |___________|_____
# ^    ^           ^    ^ 
# x0  x1          x2    x3

#x0 = 0 per definizione

#x1
TARGET_POSITION_MICROMETERS=350.0

#x2 = x1+target thickness
THICKNESS_MICROMETERS=2.0

#x3 = x0+grid length
TOTAL_GRID_LENGTH_MICROMETERS=400


##########################################
##########################################
#        don't touch from here on        #
##########################################
##########################################

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

rm -f ${INPUT_FILE}
touch ${INPUT_FILE}

printf " &picohd\n" >> ${INPUT_FILE}

printf "  trun=${TMAX_FEMTOSECONDS}\n"       >> ${INPUT_FILE}
printf "  nx=${GRID_POINTS}\n"               >> ${INPUT_FILE}
printf "  ne=${TOTAL_NUMBER_OF_ELECTRONS}\n" >> ${INPUT_FILE}
printf "  ni=${TOTAL_NUMBER_OF_IONS}\n"      >> ${INPUT_FILE}
printf "  iunits=${iunits}\n"                >> ${INPUT_FILE}

printf "  a0=${a0}\n"     >> ${INPUT_FILE}
printf "  xlambda=${LASER_WAVELENGTH_MICROMETERS}\n"     >> ${INPUT_FILE}
printf "  tpulse=${LASER_PULSE_DURATION_FEMTOSECONDS}\n" >> ${INPUT_FILE}
printf "  theta0=${LASER_ANGLE_OF_INCIDENCE}\n"          >> ${INPUT_FILE}
printf "  cpolzn='${LASER_POLARIZATION}'\n"              >> ${INPUT_FILE}

printf "  miome=${MASS_RATIO_ION_OVER_ELECTRON}\n"   >> ${INPUT_FILE}
printf "  nonc=${DENSITY_IN_nc_UNITS}\n"             >> ${INPUT_FILE}
printf "  Z=${Z}\n"                                  >> ${INPUT_FILE}
printf "  amass=${A}\n"                              >> ${INPUT_FILE}
printf "  Te=${INITIAL_ELECTRON_TEMPERATURE_KEV}\n"  >> ${INPUT_FILE}
printf "  Ti=${INITIAL_ION_TEMPERATURE_KEV}\n"       >> ${INPUT_FILE}

printf "  xl=${TOTAL_GRID_LENGTH_MICROMETERS}\n" >> ${INPUT_FILE}
printf "  target_config=${TARGET_CONFIG}\n"      >> ${INPUT_FILE}
printf "  inprof=${INPROF}\n"                    >> ${INPUT_FILE}
printf "  dfoil=${THICKNESS_MICROMETERS}\n"      >> ${INPUT_FILE}
printf "  xm1=${TARGET_POSITION_MICROMETERS}\n"  >> ${INPUT_FILE}

printf "  iout=${NUMBER_OF_OUTPUTS}\n"                        >> ${INPUT_FILE}
printf "  igr=${SNAPSHOT_FREQUENCY}\n"                        >> ${INPUT_FILE}
printf "  igmovie=${MOVIE_SNAPSHOTS}\n"                       >> ${INPUT_FILE}
printf "  itc=${NUMBER_OF_HISTORY_PLOTS}\n"                   >> ${INPUT_FILE}
printf "  ncyc=${NUMBER_OF_CYCLES_FOR_TIMEAVERAGED_PLOTS}\n"  >> ${INPUT_FILE}
printf "  ilas=${LASER_TYPE}\n"                               >> ${INPUT_FILE}

printf "  ipbc=${PARTICLE_BOUNDARYCONDITIONS}\n"  >> ${INPUT_FILE}
printf "  ifbc=${FIELD_BOUNDARYCONDITIONS}\n"     >> ${INPUT_FILE}

printf "  umevmax=${ELECTRON_EMAX_MEV}\n"   >> ${INPUT_FILE}
printf "  uimax=${ION_EMAX_MEV}\n"          >> ${INPUT_FILE}
printf "  upmax=${PROTON_EMAX_MEV}\n"       >> ${INPUT_FILE}
printf "  igxs=${SUBSAMPLE_GRID_POINTS}\n"  >> ${INPUT_FILE}
printf "  ipskip=${SUBSAMPLE_PARTICLES}\n"  >> ${INPUT_FILE}

printf "  nuav=${NUMBER_OF_TIMESTEPS_FOR_AVERAGED_DISTRIBUTION_FUNCTIONS}\n"  >> ${INPUT_FILE}
printf "  nftem=${FREQUENCY_OF_FOURIER_TRANSFORM_PLOTS}\n"                    >> ${INPUT_FILE}
printf "  ift=${SKIP_FACTOR_FT}\n"                                            >> ${INPUT_FILE}
printf "  omegm=${FREQUENCY_MAX_FT_PLOTS}\n"                                  >> ${INPUT_FILE}
printf "  ifbin=${BINNING_FT_PLOTS}\n"                                        >> ${INPUT_FILE} 



echo 'Running bops ..'
$BOPS
#$ODPP ${TOOLS_ID} $RUN ${SIM_TYPE} 9 y
echo 'Finished run'
