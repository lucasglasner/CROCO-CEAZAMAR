#!/ceaza/lucas/miniconda3/envs/main/bin/python3
'''
 # @ Author: Your name
 # @ Create Time: 2023-05-12 11:27:35
 # @ Modified by: Your name
 # @ Modified time: 2023-05-12 11:27:46
 # @ Description:

 
 This script is responsible of creating the forecast lateral forcing files. 
 
 '''
 

 # ---------------------------------- IMPORTS --------------------------------- #
import os
import datetime
import xarray as xr
import numpy as np
import pandas as pd
import shutil
from utils import add_itolap_bry

# ------------------------------- GENERAL STUFF ------------------------------ #
# crocotools_param.m static parameters
# Should be the same in ./FORECAST/crocotools_param.m !!!
RUN_dir         = '/ceaza/lucas/CROCO-CEAZAMAR/FORECAST/'
DATADIR         = '/ceaza/lucas/CROCO-CEAZAMAR/DATA/' 
CROCO_files_dir = '/ceaza/lucas/CROCO-CEAZAMAR/FORECAST/CROCO_FILES/'
Yorig           = 1950        
hdays           = 6
fdays           = 10

# GLOBAL PARAMETERS
maindir         = '/ceaza/lucas/CROCO-CEAZAMAR/'
SCRATCH_dir     = RUN_dir+'/SCRATCH/'
today           = datetime.datetime.utcnow()
fprefix         = 'crococeazaf'
itolap_mercator = 2  
itolap_variables   = [
    'zeta_west','vbar_west','ubar_west','v_west','u_west','temp_west','salt_west',
    'zeta_north','vbar_north','ubar_north','v_north','u_north','temp_north','salt_north',
    'zeta_south','vbar_south','ubar_south','v_south','u_south','temp_south','salt_south',
    'bry_time']

dates  = pd.date_range(
    (today-pd.Timedelta(days=hdays)).strftime('%F'),
    (today+pd.Timedelta(days=fdays)).strftime('%F')
    )
os.chdir(maindir)


# --------------------------------- FUNCTIONS -------------------------------- #

def load_mercator_hindcast(dates, hindcastdirectory):
    """
    Simple function for loading mercator hindcast for desired dates.
    Hindcast files should be of the kind "yyyy-mm-dd.nc"

    Args:
        dates (array): times vector
        hindcastdirectory (str): path to hindcast directory

    Returns:
        hindcast (XDataset): xarray dataset
    """
    print('Loading mercator hindcast data...')
    paths = [hindcastdirectory+d.strftime('%F')+'.nc'
             for d in dates]
    hindcast = xr.open_mfdataset(paths,
                                engine='netcdf4')
    return hindcast

def load_mercator_forecast(date, forecastdirectory):
    """
    Simple function for loading mercator forecast for the desired date.
    Forecast file should be of the kind "yyyy-mm-dd.nc"

    Args:
        dates (datetime): datetime object
        forecastdirectory (str): path to forecast directory

    Returns:
        forecast (XDataset): xarray dataset
    """
    print('Loading mercator forecast data...')
    path = forecastdirectory+date.strftime('%F')+'.nc'
    forecast = xr.open_dataset(path)
    return forecast

def load_mercator(dates, hdays, fdays,
                  scratchdirectory,
                  hindcastdirectory,
                  forecastdirectory):
    """
    Simple function for loading the raw motu hindcast and forecast file. 
    The routine concatenates the data along the time dimension and write 
    the into the scratch directory for the crocotools interpolation 
    routine to grab. 

    Args:
        dates (): time vector of hindcast/forecast dates
        hdays (int): number of hindcast days
        fdays (int): number of forecast days
        hindcastdirectory (str): path to hindcast directory
        forecastdirectory (str): path to forecast directory

    Returns:
        mercator (XDataset): xarray dataset
    """
    
    encoding_opts = {
        'time': {'units':"hours since 1950-01-01",
                 'calendar':'gregorian',
                 'dtype':float}
        }
    
    hindcast = load_mercator_hindcast(dates[:hdays], hindcastdirectory)
    forecast = load_mercator_forecast(dates[hdays:][0], forecastdirectory)
    print('Merging data in a single dataset...')
    forecast = forecast.isel(time=slice(0, fdays))
    mercator = xr.concat([hindcast, forecast], 'time').drop_duplicates('time')
    print('Writing raw mercator file on disk...')
    mercator.load().to_netcdf(scratchdirectory+dates[hdays:][0].strftime('%F')+'.nc',
                              encoding=encoding_opts)
    return mercator

def make_forecast_mercator(
            matlabexec='/opt/MATLAB/R2017a/bin/matlab -nodisplay -singleCompThread',
            execdirectory=RUN_dir):
    """
    This function is responsable of running adapted crocotools routine
    make_forecast_mercator.m. The function changes directory to the execution dir,
    then run the start.m file, and run the crocotools routine as usual.
    For this function to work is important to have a consistent
    crocotools_param.m file and configuration (crocotools paths, variables, etc).
    
    Args:
        matlabexec (str, optional):
            path to matlab executable with respective flags.
            Defaults to 'matlab -nodisplay -singleCompThread'.
        execdirectory (str, optional): 
            Directory where matlab is executed.
            Defaults to Run_dir (see global variable above).
    """
    os.chdir(execdirectory)
    routines     = ['start', 'make_forecast_mercator']
    runarguments = 'clear all; '+'; '.join(routines)+'; exit;'
    command      = matlabexec+' -r "{}"'.format(runarguments)
    print(command)
    os.system(command)
    os.chdir(maindir)
    return 



def main_bry_forecast():
    """
    This function just run all the previous routines and creates the file.
    """
    starttime = datetime.datetime.utcnow()
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Creating mercator raw motu for croco_tools, please wait...        ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    if os.path.isfile(SCRATCH_dir+'/'+today.strftime('%F')+'.nc'):
        print(SCRATCH_dir+'/'+today.strftime('%F')+'.nc', 'already exists !!  ')
    else:
        load_mercator(dates, hdays, fdays, SCRATCH_dir, 
                    DATADIR+'/MERCATOR/HINDCAST/',
                    DATADIR+'/MERCATOR/')
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Running crocotools make_forecast_mercator.m, please wait...       ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    make_forecast_mercator()
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Adding overlap days to croco_bry files, please wait...            ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')

    add_itolap_bry(today, itolap=itolap_mercator,
                variables=itolap_variables,
                inputfiledir=CROCO_files_dir,
                outputfiledir=SCRATCH_dir,
                fprefix=fprefix+'_bry_')

    bryname = SCRATCH_dir+fprefix+'_bry_'+today.strftime('%Y%m%d')+'.nc'
    if os.path.isfile(bryname):
        print('Overwriting file...',bryname.replace(SCRATCH_dir,''),'     ')
        shutil.move(bryname,bryname.replace(SCRATCH_dir,CROCO_files_dir))     
    print('-------------------------------------------------------------------')
    print(datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'          ')
    endtime = datetime.datetime.utcnow()
    print('Elapsed time:',endtime-starttime)
    print('All good','                                                        ')
    return

# ------------------------------- RUN ROUTINES ------------------------------- #
if __name__=='__main__':
    todayfile = CROCO_files_dir+fprefix+'_bry_'+today.strftime('%Y%m%d')+'.nc'
    if os.path.isfile(todayfile):
        print('Today file:',todayfile,'already exists!')
    else:
        pass
    main_bry_forecast()


