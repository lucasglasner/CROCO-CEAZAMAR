#!/usr/bin/bash
#SBATCH -p part1
#SBATCH -J spinup
#SBATCH -n 60
#SBATCH --output=run_crocodaily.log


source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
export NETCDF=/opt/wrf/LIBS/netcdf
export PATH=$NETCDF/bin::${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${NETCDF}/lib

#################################################################################################################################
#                                                    USER PARAMETERS                                                            #
#################################################################################################################################
SIMNAME='crococeazah'                                                           # Simulation prefix name (e.g croco)
MAINDIR=/ceaza/lucas/CROCO-CEAZAMAR                                             # Directory of this script
TIMESTEP=120                                                                    # Desired time step in seconds
NUMTIMES=$(expr 86400 / $TIMESTEP)                                              # Number of timesteps to run (1 day)
RUNCMD='mpirun '                                                                # Command for running the model executable
cd $MAINDIR

INIDATE='2020-01-02'
ENDDATE='2020-12-31'
SCRATCHDIR=${MAINDIR}/HINDCAST/OUTPUT/spinup                                    # Directory where model is run
OUTPUTDIR=${MAINDIR}/HINDCAST/OUTPUT                                            # Directory where model outputs are saved
CROCOFILESDIR=${MAINDIR}/HINDCAST/CROCO_FILES                                   # Directory where croco forcing are stored
CROCOEXEC=${MAINDIR}/HINDCAST/croco                                             # croco executable name

CROCOIN=${MAINDIR}/HINDCAST/${SIMNAME}.in                                       # croco.in textfile path
CROCOGRID=${CROCOFILESDIR}/${SIMNAME}_grd.nc                                    # path to model grid                                
NCCOPY="/ceaza/lucas/miniconda3/envs/main/bin/nccopy -k classic"                # path to nccopy binary (netcdf package)
RESTART=1                                                                       # 0 - 1 for use restart or ini file
#################################################################################################################################
#                                                    END OF USER SELECTIONS                                                     #
#################################################################################################################################
printf -- '#%.0s' {1..160}
printf "\n"
printf -- ' %.0s' {1..60}
echo RUNNING CROCO MODEL FROM DAILY FORCINGS
printf -- '#%.0s' {1..160}
printf "\n\n"
# Copy model executable to run directory
if test -f "$CROCOEXEC"; then
    cp $CROCOEXEC $SCRATCHDIR
else
    echo croco executable not found !! Dont forget to compile croco first !!
    exit
fi


# Getting number of days between initial and final date
NDAYS=$(echo \($(date -d "$ENDDATE" +%s) - $(date -d "$INIDATE" +%s)\) / 86400 | bc) 
INIYEAR=$(echo $INIDATE | cut -c -4)
# MODEL PREFFIX
MODEL=${SCRATCHDIR}/${SIMNAME}
for i in $(seq 0 $NDAYS); do
    echo Running $(date -d "${INIDATE} +$i days" +%F)...
    TIME=$(date -d "${INIDATE} +$i days" +%Y%m%d)
    # Getting forcing file names
    CROCOBLK=${CROCOFILESDIR}/${SIMNAME}_blk_${TIME}.nc
    CROCOBRY=${CROCOFILESDIR}/${SIMNAME}_bry_${TIME}.nc
    CROCOFRC=${CROCOFILESDIR}/${SIMNAME}_frc_${TIME}.nc
    if [ $RESTART -eq "0" ]; then
        CROCOINI=${CROCOFILESDIR}/${SIMNAME}_ini_${TIME}.nc
    else
        CROCOINI=${MODEL}_rst_$(date -d "${INIDATE} +$(expr $i - 1) days" +%Y%m%d).nc
    fi

    # Copy model input files to run directory
    echo Copy grd file: cp $CROCOGRID ${MODEL}_grd.nc
    $NCCOPY $CROCOGRID ${MODEL}_grd.nc
    echo Copy blk file: cp $CROCOBLK ${MODEL}_blk.nc
    $NCCOPY $CROCOBLK ${MODEL}_blk.nc
    echo Copy bry file: cp $CROCOBRY ${MODEL}_bry.nc
    $NCCOPY $CROCOBRY ${MODEL}_bry.nc
    echo Copy frc file: cp $CROCOFRC ${MODEL}_frc.nc
    $NCCOPY $CROCOFRC ${MODEL}_frc.nc
    echo Copy ini file: cp $CROCOINI ${MODEL}_ini.nc
    $NCCOPY $CROCOINI ${MODEL}_ini.nc
    echo Copy .in file: cp $CROCOIN ${MODEL}_${TIME}.in
    cp $CROCOIN ${MODEL}_${TIME}.in
    
    # Edit croco.in file with desired parameters (timestep and number of runtimes)
    echo "Edit ${MODEL}_${TIME}.in: TIMESTEP = ${TIMESTEP}s ; NUMTIMES = ${NUMTIMES}"
    sed -i "s/%NUMTIMES%/${NUMTIMES}/g" ${MODEL}_${TIME}.in
    sed -i "s/%TIMESTEP%/${TIMESTEP}/g" ${MODEL}_${TIME}.in
    
    # Run the model
    echo Running the model...
    cd $SCRATCHDIR
    echo "${RUNCMD}./croco ${SIMNAME}_${TIME}.in > ${SIMNAME}_${TIME}.out"
    ${RUNCMD}./croco ${SIMNAME}_${TIME}.in > ${SIMNAME}_${TIME}.out
    cd $MAINDIR
    printf "\n"

    # Test if the run has finised properly
    if test -f ${MODEL}_${TIME}.out; then
        echo "Test ${MODEL}_${TIME}.out"
        status=`tail -2 ${MODEL}_${TIME}.out | grep DONE | wc -l`
        if [[ $status == 1 ]]; then
            echo "¡ All good !"
        else
            echo "Warning: run not finished properly"
            tail -20 ${MODEL}_${TIME}.out
            echo "Run for ${TIME} didnt work"
            exit 1
        fi
    fi
    echo Saving outputs...
    mv -f ${MODEL}_rst.nc ${MODEL}_rst_${TIME}.nc
    mv -f ${MODEL}_avg.nc ${OUTPUTDIR}/${SIMNAME}_avg_${TIME}.nc
    mv -f ${MODEL}_his.nc ${OUTPUTDIR}/${SIMNAME}_his_${TIME}.nc

    # Clean netcdf from scratch directory
    rm -f ${MODEL}_rst.nc
    rm -f ${MODEL}_blk.nc
    rm -f ${MODEL}_bry.nc
    rm -f ${MODEL}_frc.nc
    rm -f ${MODEL}_ini.nc

    echo Done

    printf -- '-%.0s' {1..160}
    printf "\n"
done
exit