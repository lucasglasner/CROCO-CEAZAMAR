#!/ceaza/lucas/miniconda3/envs/main/bin/python3
'''
 # @ Author: Your name
 # @ Create Time: 2023-04-25 15:27:33
 # @ Modified by: Your name
 # @ Modified time: 2023-04-25 15:27:41
 # @ Description:
 
 This script is responsible of creating the croco_bry files for the hindcast simulation. 
 For the moment this script just call the adapted matlab routines from crocotools. 
 Especially the "make_hindcast_mercator.m" routine in croco_tools/Forecast_tools/
 
 '''

 # ---------------------------------- IMPORTS --------------------------------- #
import os
import datetime
import numpy as np
import pandas as pd
from glob import glob


# ------------------------------- GENERAL STUFF ------------------------------ #
# crocotools_param.m static parameters
# Should be the same in ./HINDCAST/crocotools_param.m !!!
RUN_dir         = '/ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/'
DATADIR         = '/ceaza/lucas/CROCO-CEAZAMAR/DATA/' 
CROCO_files_dir = '/ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/CROCO_FILES/'
SCRATCH_dir     = RUN_dir+'/SCRATCH/'
Yorig           = 1950       
MERCATOR_delay  = 6
MERCATOR_offset = 1
 
 
# GLOBAL PARAMETERS
itolap_mercator = 2  
maindir         = '/ceaza/lucas/CROCO-CEAZAMAR/'
today           = datetime.datetime.utcnow()
fprefix         = 'crococeazah'
dates  = pd.date_range(
    (today-pd.Timedelta(days=MERCATOR_delay+MERCATOR_offset)).strftime('%F'),
    (today-pd.Timedelta(days=MERCATOR_delay)).strftime('%F')
    )

os.chdir(maindir)


# --------------------------------- FUNCTIONS -------------------------------- #
def make_hindcast_mercator(matlabexec='matlab -nodisplay -nosplash -nodesktop',
            execdirectory=RUN_dir):
    """
    This function is responsable of running adapted crocotools routine
    make_hindcast_mercator.m. The functions changes directory to the execution dir,
    then run the start.m file, and run the crocotools routine as usual.
    For this function to work is important to have a consistent
    crocotools_param.m file and configuration (crocotools paths, variables, etc).
    
    I have modified the standard make_GFS.m routine to include tides and to
    be consistent with the filenames used in this forecast package

    Args:
        matlabexec (str, optional):
            path to matlab executable with respective flags.
            Defaults to 'matlab -nodisplay -nosplash -nodesktop'.
        execdirectory (str, optional): 
            Directory where matlab is executed.
            Defaults to Run_dir (see global variable above).
    """
    os.chdir(execdirectory)
    routines     = ['start', 'make_hindcast_mercator']
    runarguments = 'clear all; '+'; '.join(routines)+'; exit;'
    command      = matlabexec+' -r "{}"'.format(runarguments)
    print(command)
    os.system(command)
    os.chdir(maindir)
    return 

def main_bry_hindcast():
    """
    This function just run all the previous routines and creates the file.
    """
    starttime = datetime.datetime.utcnow()
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Running crocotools make_hindcast_mercator.m, please wait...       ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    make_hindcast_mercator()
    print('Cleaning scratch directory...                                      ')
    for f in glob(SCRATCH_dir+'/*.nc'):
        os.remove(f)
    print('-------------------------------------------------------------------')
    print(datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'          ')
    endtime = datetime.datetime.utcnow()
    print('Elapsed time:',endtime-starttime)
    print('All good','                                                        ')
    return

if __name__=='__main__':
    main_bry_hindcast()