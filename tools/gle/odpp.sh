BOPS=$HOME/bops/src/bops.exe
OD2GLE=$HOME/bops/tools/gle/od2gle.exe
GLE="/c/MinGW/gle/bin/gle.exe"
#GLE=gle

ID_FOLDER=$1
RUN_NAME=$2
ID_NAME=$3
PARAM=$4
TEST=$5

echo 'Postprocessing ...'
mkdir plots
echo "Running odplot"
echo 'odplot -> gle format'
rm -f *.gle
printf "${ID_FOLDER}\n" > od.in
printf "${RUN_NAME}\n" >> od.in
printf "${ID_NAME}\n" >> od.in
printf "$PARAM\n" >> od.in
printf "$TEST\n" >> od.in

$OD2GLE < ./od.in
echo '... done'
echo
echo 'Making plots ...'
rm -f plots/*
ls *.gle  > glelist
    cat glelist | while read PLOT
    do
      $GLE -d jpg -o plots/$PLOT $PLOT
 #     $GLE -d pdf -dpi 200 -o plots/$PLOT $PLOT
    done


