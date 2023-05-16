#!/usr/bin/bash
#SBATCH -p part1
#SBATCH -J make_bc
#SBATCH -N 1
#SBATCH --output=/ceaza/lucas/CROCO-CEAZAMAR/make_bc.log
#
# This shell script purpose is to run all the python preprocessing scripts for creating the bulk and bry files of the hindcast and forecast run.
#
maindir='/ceaza/lucas/CROCO-CEAZAMAR/'
cd $maindir
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CREATING CROCO BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

echo 'Checking if global forecasts have arrived...'
ohindcast=$maindir/DATA/MERCATOR/HINDCAST/$(date -d '6 day ago' +%F -u).nc
ahindcast=$maindir/DATA/ERA5/$(date -d '6 day ago' +%F -u).nc
oforecast=$maindir/DATA/MERCATOR/$(date +%F -u).nc
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
SECONDS=0
now=$(date +%F\ %H:%M:%S -u)
if [ ! -f $ahindcast ]; then
    echo "Atmospheric hindcast file $ahindcast doesnt exists!"
else
    echo $now Creating hindcast atmospheric boundary conditions...
    ./make_blk_hindcast.py
fi
printf "\n\n"
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
now=$(date +%F\ %H:%M:%S -u)
if [ ! -f $ohindcast ]; then
    echo "Ocean hindcast file $ohindcast doesnt exists!"
else
    echo $now Creating hindcast oceanic boundary conditions...
    ./make_bry_hindcast.py
fi
printf "\n\n"
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo 'Cleaning hindcast scratch directory...'
rm -rf $maindir/HINDCAST/SCRATCH/*
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
now=$(date +%F\ %H:%M:%S -u)
echo $now  Creating forecast atmospheric boundary conditions...
./make_blk_forecast.py
printf "\n\n"
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
now=$(date +%F\ %H:%M:%S -u)
if [ ! -f $oforecast ]; then
    echo "Ocean forecast file $oforecast doesnt exists!"
else
    echo $now Creating forecast oceanic boundary conditions...
    ./make_bry_forecast.py
fi
printf "\n\n"
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo 'Cleaning forecast scratch directory...'
rm -rf $maindir/FORECAST/SCRATCH/*

echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo "Done"
echo "Elapsed time: $(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
exit