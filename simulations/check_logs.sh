#! /bin/bash
declare -a logs=( LOG_foil_fsmu.txt LOG_foil_hb.txt LOG_foil_ls.txt LOG_foil_ramp.txt LOG_foil_tnsa.txt LOG_foil.txt LOG_multiplayer.txt LOG_multi_robinson.txt LOG_nl_surface.txt LOG_resabs.txt LOG_rpa_nanofoil.txt LOG_snells_law.txt LOG_vacuum.txt LOG_wake.txt )
for ((i = 0; i < 14 ; i++)) ; do
#tail -2 ${logs[i]} | head -1
tail -1 ${logs[i]}
done
