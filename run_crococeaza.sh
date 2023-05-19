#!/usr/bin/bash
#SBATCH -p part2
#SBATCH -J crocoF
#SBATCH -n 60
#SBATCH --output=run_crococeaza.log

# set -e
source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
export NETCDF=/opt/wrf/LIBS/netcdf
export PATH=$NETCDF/bin::${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${NETCDF}/lib

#################################################################################################################################
#                                                    USER PARAMETERS                                                            #
#################################################################################################################################
RUNHINDCAST=0                                                                   # Initialize forecast from OGCM or CROCOHINDCAST
HSIMNAME='crococeazah'                                                          # Hindcast simulation prefix name (e.g crocoh)
FSIMNAME='crococeazaf'                                                          # Forecast simulation prefix name (e.g crocoh)
MAINDIR=/ceaza/lucas/CROCO-CEAZAMAR                                             # Directory of this script
TIMESTEP=120                                                                    # Desired time step in seconds
HDAYS=6
FDAYS=10
RUNCMD='mpirun '                                                                # Command for running the model executable
cd $MAINDIR

RUNDATE=$(date -u +%Y%m%d)
INIYEAR=$(echo $RUNDATE | cut -c -4)

FNUMTIMES=$((( 86400 / $TIMESTEP )*( $HDAYS + $FDAYS )))                        # Number of timesteps to run
HNUMTIMES=$(( 86400 / $TIMESTEP ))

SCRATCHDIR=${MAINDIR}/FORECAST/OUTPUT/SCRATCH                                   # Directory where model is run
FOUTPUTDIR=${MAINDIR}/FORECAST/OUTPUT                                           # Directory where forecast model outputs are saved
HOUTPUTDIR=${MAINDIR}/HINDCAST/OUTPUT                                           # Directory where hindcast model outputs are saved

FCROCOFILESDIR=${MAINDIR}/FORECAST/CROCO_FILES                                  # Directory where croco forecast forcing are stored
HCROCOFILESDIR=${MAINDIR}/HINDCAST/CROCO_FILES                                  # Directory where croco hindcast forcing are stored

FCROCOEXEC=${MAINDIR}/FORECAST/croco                                            # croco forecast executable path
HCROCOEXEC=${MAINDIR}/HINDCAST/croco                                            # croco hindcast executable path

FCROCOIN=${MAINDIR}/FORECAST/${FSIMNAME}.in                                     # croco.in forecast textfile path
HCROCOIN=${MAINDIR}/HINDCAST/${HSIMNAME}.in                                     # croco.in hindcast textfile path
                             
NCCOPY="/ceaza/lucas/miniconda3/envs/main/bin/nccopy -k classic"                # path to nccopy binary (netcdf package)
#################################################################################################################################
#                                                    END OF USER SELECTIONS                                                     #
#################################################################################################################################

if [ $RUNHINDCAST -eq "1" ]; then
    printf -- '#%.0s' {1..160}
    printf "\n"
    printf -- ' %.0s' {1..60}
    echo RUNNING CROCO HINDCAST MODEL
    printf -- '#%.0s' {1..160}
    printf "\n\n"
    # Copy model executable to run directory
    if test -f "$HCROCOEXEC"; then
        cp $HCROCOEXEC ${HOUTPUTDIR}/
    else
        echo croco executable not found !! Dont forget to compile croco first !!
        exit
    fi


    # MODEL PREFFIX
    MODEL=${HOUTPUTDIR}/${HSIMNAME}
    echo Running $(date -d "${RUNDATE} -$HDAYS days" +%F)...
    TIME=$(date -d "${RUNDATE} -$HDAYS days" +%Y%m%d)
    # Getting forcing file names
    CROCOBLK=${HCROCOFILESDIR}/${HSIMNAME}_blk_${TIME}.nc
    CROCOBRY=${HCROCOFILESDIR}/${HSIMNAME}_bry_${TIME}.nc
    CROCOFRC=${HCROCOFILESDIR}/${HSIMNAME}_frc_${TIME}.nc
    CROCOGRD=${HCROCOFILESDIR}/${HSIMNAME}_grd.nc
    CROCOINI=${MODEL}_rst_$(date -d "${TIME} -1 days" +%Y%m%d).nc

    # Copy model input files to run directory
    echo Copy grd file: cp $CROCOGRD ${MODEL}_grd.nc
    $NCCOPY $CROCOGRD ${MODEL}_grd.nc
    echo Copy blk file: cp $CROCOBLK ${MODEL}_blk.nc
    $NCCOPY $CROCOBLK ${MODEL}_blk.nc
    echo Copy bry file: cp $CROCOBRY ${MODEL}_bry.nc
    $NCCOPY $CROCOBRY ${MODEL}_bry.nc
    echo Copy frc file: cp $CROCOFRC ${MODEL}_frc.nc
    $NCCOPY $CROCOFRC ${MODEL}_frc.nc
    echo Copy ini file: cp $CROCOINI ${MODEL}_ini.nc
    $NCCOPY $CROCOINI ${MODEL}_ini.nc
    echo Copy .in file: cp $HCROCOIN ${MODEL}_${TIME}.in
    cp $HCROCOIN ${MODEL}_${TIME}.in

    # Edit croco.in file with desired parameters (timestep and number of runtimes)
    echo "Edit ${MODEL}_${TIME}.in: TIMESTEP = ${TIMESTEP}s ; NUMTIMES = ${HNUMTIMES}"
    sed -i "s/%NUMTIMES%/${HNUMTIMES}/g" ${MODEL}_${TIME}.in
    sed -i "s/%TIMESTEP%/${TIMESTEP}/g" ${MODEL}_${TIME}.in

    # Run the model
    echo Running the model...
    cd $HOUTPUTDIR
    echo "cd $HOUTPUTDIR"
    echo "${RUNCMD}./croco ${HSIMNAME}_${TIME}.in > ${HSIMNAME}_${TIME}.out"
    # ${RUNCMD}./croco ${HSIMNAME}_${TIME}.in > ${HSIMNAME}_${TIME}.out
    echo "cd $MAINDIR"
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
    # mv -f ${MODEL}_rst.nc ${MODEL}_rst_${TIME}.nc
    # mv -f ${MODEL}_avg.nc ${MODEL}_avg_${TIME}.nc
    # mv -f ${MODEL}_his.nc ${MODEL}_his_${TIME}.nc

    # Clean netcdf from scratch directory
    rm -f ${MODEL}_rst.nc
    rm -f ${MODEL}_blk.nc
    rm -f ${MODEL}_bry.nc
    rm -f ${MODEL}_frc.nc
    rm -f ${MODEL}_ini.nc
    rm -f ${MODEL}_grd.nc

    echo Done

    printf -- '-%.0s' {1..160}
    printf "\n"
fi

#################################################################################################################################
#                                                    DONE WITH HINDCAST SIMULATION                                              #
#################################################################################################################################

printf -- '#%.0s' {1..160}
printf "\n"
printf -- ' %.0s' {1..60}
echo RUNNING CROCO FORECAST MODEL
printf -- '#%.0s' {1..160}
printf "\n\n"
# Copy model executable to run directory
if test -f "$FCROCOEXEC"; then
    cp $FCROCOEXEC ${FOUTPUTDIR}/
else
    echo croco executable not found !! Dont forget to compile croco first !!
    exit
fi


# MODEL PREFFIX
MODEL=${FOUTPUTDIR}/${FSIMNAME}
echo Running forecast for: $(date -d "${RUNDATE}" +%F)...
TIME=$(date -d "${RUNDATE}" +%Y%m%d)
# Getting forcing file names
CROCOBLK=${FCROCOFILESDIR}/${FSIMNAME}_blk_${TIME}.nc
CROCOBRY=${FCROCOFILESDIR}/${FSIMNAME}_bry_${TIME}.nc
CROCOFRC=${FCROCOFILESDIR}/${FSIMNAME}_frc_${TIME}.nc
CROCOGRD=${FCROCOFILESDIR}/${FSIMNAME}_grd.nc
if [ $RUNHINDCAST -eq "1" ]; then
    CROCOINI=${HOUTPUTDIR}/${HSIMNAME}_rst_$(date -d "${RUNDATE} -$HDAYS days" +%Y%m%d).nc
else
    CROCOINI=${FCROCOFILESDIR}/${FSIMNAME}_ini_${TIME}.nc
fi

# Copy model input files to run directory
echo Copy grd file: cp $CROCOGRD ${MODEL}_grd.nc
$NCCOPY $CROCOGRD ${MODEL}_grd.nc
echo Copy blk file: cp $CROCOBLK ${MODEL}_blk.nc
$NCCOPY $CROCOBLK ${MODEL}_blk.nc
echo Copy bry file: cp $CROCOBRY ${MODEL}_bry.nc
$NCCOPY $CROCOBRY ${MODEL}_bry.nc
echo Copy frc file: cp $CROCOFRC ${MODEL}_frc.nc
$NCCOPY $CROCOFRC ${MODEL}_frc.nc
echo Copy ini file: cp $CROCOINI ${MODEL}_ini.nc
$NCCOPY $CROCOINI ${MODEL}_ini.nc
echo Copy .in file: cp $FCROCOIN ${MODEL}_${TIME}.in
cp $FCROCOIN ${MODEL}_${TIME}.in

# Edit croco.in file with desired parameters (timestep and number of runtimes)
echo "Edit ${MODEL}_${TIME}.in: TIMESTEP = ${TIMESTEP}s ; NUMTIMES = ${FNUMTIMES}"
sed -i "s/%NUMTIMES%/${FNUMTIMES}/g" ${MODEL}_${TIME}.in
sed -i "s/%TIMESTEP%/${TIMESTEP}/g" ${MODEL}_${TIME}.in

 # Run the model
echo Running the model...
cd $FOUTPUTDIR
echo "cd $FOUTPUTDIR"
echo "${RUNCMD}./croco ${FSIMNAME}_${TIME}.in > ${FSIMNAME}_${TIME}.out"
${RUNCMD}./croco ${FSIMNAME}_${TIME}.in > ${FSIMNAME}_${TIME}.out
echo "cd $MAINDIR"
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
mv -f ${MODEL}_avg.nc ${MODEL}_avg_${TIME}.nc
# mv -f ${MODEL}_his.nc ${MODEL}_his_${TIME}.nc

# Clean netcdf from scratch directory
rm -f ${MODEL}_rst.nc
rm -f ${MODEL}_blk.nc
rm -f ${MODEL}_bry.nc
rm -f ${MODEL}_frc.nc
rm -f ${MODEL}_ini.nc
rm -f ${MODEL}_grd.nc

echo Done

printf -- '-%.0s' {1..160}
printf "\n"

exit
