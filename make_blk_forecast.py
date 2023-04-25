#!/ceaza/lucas/miniconda3/envs/main/bin/python3
'''
 # @ Author: Your name
 # @ Create Time: 2023-04-24 15:51:58
 # @ Modified by: Your name
 # @ Modified time: 2023-04-24 15:52:10
 # @ Description:
 
 This script is responsible of creating the forecast surface forcing files
 
 '''
 
 # ---------------------------------- IMPORTS --------------------------------- #
import os
import datetime
import numpy as np
import pandas as pd

# ------------------------------- GENERAL STUFF ------------------------------ #
# crocotools_param.m static parameters
# Should be the same in crocotools_param.m !!!
RUN_dir         = '/ceaza/lucas/CROCO-CEAZAMAR/FORECAST/'
DATADIR         = '/ceaza/lucas/CROCO-CEAZAMAR/DATA/' 
CROCO_files_dir = '/ceaza/lucas/CROCO-CEAZAMAR/FORECAST/CROCO_FILES/'
Yorig           = 1950        
hdays           = 6
fdays           = 10

maindir         = '/ceaza/lucas/CROCO-CEAZAMAR/'
scratchdir      = RUN_dir+'/SCRATCH/'
today           = datetime.datetime.utcnow()
fprefix         = 'crococeazaf'

dates  = pd.date_range(
    (today-pd.Timedelta(days=hdays)).strftime('%F'),
    (today+pd.Timedelta(days=fdays)).strftime('%F')
    )
os.chdir(maindir)


# --------------------------------- FUNCTIONS -------------------------------- #
def make_forecast_GFS(matlabexec='matlab -nodisplay -nosplash -nodesktop',
            execdirectory=RUN_dir):
    """
    This function is responsable of running matlab crocotools routine
    make_GFS.m The functions changes directory to the execution dir,
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
    routines     = ['start', 'make_GFS']
    runarguments = 'clear all; '+'; '.join(routines)+'; exit;'
    command      = matlabexec+' -r "{}"'.format(runarguments)
    print(command)
    os.system(command)
    os.chdir(maindir)
    return 

def main_blk_forecast():
    """
    This function just run all the previous routines and creates the file.
    """
    starttime = datetime.datetime.utcnow()
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Running crocotools make_GFS.m, please wait...                     ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    blkname = CROCO_files_dir+fprefix+'_blk_'+today.strftime('%Y%m%d')+'.nc'
    if os.path.isfile(blkname):
            print('\t',blkname,' already exists!')
    else:
        make_forecast_GFS()
    print('-------------------------------------------------------------------')
    print(datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'          ')
    endtime = datetime.datetime.utcnow()
    print('Elapsed time:',endtime-starttime)
    print('All good','                                                        ')
    return
# ------------------------------- RUN ROUTINES ------------------------------- #
if __name__=='__main__':
    main_blk_forecast()
