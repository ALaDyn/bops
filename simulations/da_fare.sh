#!/bin/bash
./foil.sh > LOG_foil.txt 2> LOG_foil.txt &
./foil_fsmu.sh > LOG_foil_fsmu.txt 2> LOG_foil_fsmu.txt &
./foil_hb.sh > LOG_foil_hb.txt 2> LOG_foil_hb.txt &
./foil_ls.sh > LOG_foil_ls.txt 2> LOG_foil_ls.txt &
./foil_ramp.sh > LOG_foil_ramp.txt 2> LOG_foil_ramp.txt &
./foil_tnsa.sh > LOG_foil_tnsa.txt 2> LOG_foil_tnsa.txt &
./multi_robinson.sh > LOG_multi_robinson.txt 2> LOG_multi_robinson.txt &
./multilayers.sh > LOG_multiplayer.txt 2> LOG_multiplayer.txt &
./nl_surface.sh > LOG_nl_surface.txt 2> LOG_nl_surface.txt &
./resabs.sh > LOG_resabs.txt 2> LOG_resabs.txt &
./rpa_nanofoil.sh > LOG_rpa_nanofoil.txt2 > LOG_rpa_nanofoil.txt &
./snells_law.sh > LOG_snells_law.txt 2> LOG_snells_law.txt &
./vacuum.sh > LOG_vacuum.txt 2> LOG_vacuum.txt &
./wake.sh > LOG_wake.txt 2> LOG_wake.txt &
