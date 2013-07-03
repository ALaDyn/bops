#! /bin/bash
if [ $# != 1 ]
 then
 echo "scrivere ./launcher.sh seguito dal nome della simulazione (senza estensione .sh)"
 exit
fi

./$1.sh > LOG_$1.txt 2> LOG_$1.txt 
