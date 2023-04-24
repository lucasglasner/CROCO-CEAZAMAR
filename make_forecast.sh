#!/usr/bin/bash
# This shell script purpose is to run all the python preprocessing scripts for creating the bulk and bry files of the hindcast and forecast run.

maindir='/ceaza/lucas/CROCO-CEAZAMAR/'
cd $maindir

now=$(date +%F\ %H:%M:%S)
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo $now  Creating hindcast atmospheric boundary conditions...
./make_blk_hindcast.py
printf '\n\n'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
now=$(date +%F\ %H:%M:%S)
echo $now  Creating forecast atmospheric boundary conditions...
./make_blk_forecast.py
printf '\n\n'
echo Done
exit