#!/ceaza/lucas/miniconda3/envs/main/bin/python3
'''
# @ Author: lucas
# @ Create Time: 2023-04-12 10:44:56
# @ Modified by: lucas
# @ Modified time: 2023-04-12 10:45:01
# @ Description:

This script is responsible of creating the hindcast surface forcing files.
The algorithm consist of the following steps:
    (0)  Check if croco_blk files already exists. If not, then create the file,
         if yes but is corrupted somehow*, create it again and overwrite.
    (i)  Transform raw ERA5 NRT hourly data to a crocotools compatible format.
         Store this result in scratch directory.
    (ii) Run make_hindcast_ERA5.m matlab routine (croco_tools/Forecast_tools)
         and create the croco_blk.nc atmospheric forcing file and the 
         croco_frc.nc only with the tidal forcing variables. 
    (...) Add itolap to bulk files with a few hours
    (...) Clean scratch directory

In addition to the needed python packages, this script only works if a
matlab and nco (e.g ncks, ncatted, ncpdq, etc) executables exists on system path.


* No routine to check data corruption yet...
'''

 
# ---------------------------------- IMPORTS --------------------------------- #
import os
import json
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset as netcdf

# ------------------------------- GENERAL STUFF ------------------------------ #
# crocotools_param.m static parameters
# Should be the same in crocotools_param.m !!!
RUN_dir         = '/ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/'
DATADIR         = '/ceaza/lucas/CROCO-CEAZAMAR/DATA/' 
CROCO_files_dir = '/ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/CROCO_FILES/'
Yorig           = 1950              
ERA5_delay      = 6         
ERA5_offset     = 3          
itolap_era5     = 6   

maindir         = '/ceaza/lucas/CROCO-CEAZAMAR/'
scratchdir      = RUN_dir+'/SCRATCH/'
today           = datetime.datetime.utcnow()
fprefix         = 'crococeazah'

# ERA5 PARAMETERS
pathERA5raw = DATADIR+'/ERA5/'
dates  = pd.date_range(
    (today-pd.Timedelta(days=ERA5_delay+ERA5_offset)).strftime('%F'),
    (today-pd.Timedelta(days=ERA5_delay)).strftime('%F')
    )

os.chdir(maindir)

# --------------------------------- FUNCTIONS -------------------------------- #
def ERA5_coefficients():
    """
    This function returns a dictionary with 3 lists, the ER5 variable
    names, the corresponding unit conversion coefficients, and the
    final units as strings.
    
    Returns:
        dict: Dictionary with 'variables', 'conv_cff' and 'units' keys
    """
    with open('./croco_tools/Aforc_ERA5/ERA5_variables.json', 'r') as jf:
        era5 = json.load(jf)
    # TP: convert from accumlated m in a hour into   kg m-2 s-1
    cff_tp=1000./3600. # m in 1 hour -> kg m-2 s-1
    # Heat flux J m-2 in one hour into W m-2
    cff_heat=1./3600.   # J m-2 in 1 hour -> W m-2
    # Stress N m-2 accumulated in hour
    cff_stress=1./3600.

    variables = [ 'lsm'  , 'tp'        , 'e'         , 'u10'  , 'v10'  , 'ewss'    , 'nsss'    , 't2m', 'd2m', 'ssr'   , 'str'   , 'strd'  ]
    conv_cff  = [ 1.     , cff_tp      , cff_tp      , 1.     , 1.     , cff_stress, cff_stress, 1.   , 1.   , cff_heat, cff_heat, cff_heat]
    units     = [ '(0-1)', 'kg m-2 s-1', 'kg m-2 s-1', 'm s-1', 'm s-1', 'N m-2 s' , 'N m-2 s' , 'K'  , 'K'  ,'W m-2' , 'W m-2' , 'W m-2'  ]  
    long_name = [era5[n][0] for n in variables]
    return dict(variables=variables,
                conv_cff=conv_cff,
                units=units,
                long_name=long_name)
 
 
 
def ERA5_convert(fname, date, outdir=scratchdir, Yorig=Yorig):
    """
    This function is a translate of ./ERA5_convert.py
    croco_tools script for any kind of temporal data
    (not only interanual in Y????M?? format).
    Just give a raw ERA5 path in fname and the function 
    create some crocotools-compatible netcdf files in the
    outputdirectory.
    
    Args:
        fname     (str): path to raw ERA5 data.
                         File name should be something like
                         /path.../blabla_yyyymmdd.nc
                         or
                         /path.../blabla_yyyymm.nc
                         or
                         /path.../blabla_0001.nc
                         etc, etc, etc
        date      (datetime): time record in python format (datetime object)
        outdir    (str): path to crocotools-compatible ERA5 data 
                         directory
        Yorig (int, optional): Reference year for simulation.
        Defaults to Yorig=1950 see above.
    """
    fdate     = date.strftime('%F')
    coeffs    = ERA5_coefficients()
    for k in range(len(coeffs['variables'])):
        # Variable's name, long-name and level-type
        vname = coeffs['variables'][k]
        vlong = coeffs['long_name'][k]
        print('\t','Processing variable:',vname,vlong)

        # Read input filedate.toordinal(date(Yorig,1,1))
        nc = netcdf(fname,'r+',format='NETCDF4')
        time = nc.variables['time'][:]
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]+360
        data = nc.variables[vname][:,:,:]
        nc.close()
        
        # Flip latitudes (to have increasing latitudes...)
        lat = np.flip(lat, axis=0)
        data = np.flip(data, axis=1)
        
        # Missing values and multiply by cff to change unit
        try:
            mvalue=data.fill_value
        except AttributeError:
            print ('No fill value.. use nan')
            mvalue=np.nan

        data=np.array(data)
        data=coeffs['conv_cff'][k]*data
        data[np.where(data==mvalue)]=9999.
        
        # Convert time from
        # hours since 1900-01-01 00:00:00 into
        # days since Yorig-01-01 00:00:00
        time = time / 24.
        time = time - datetime.date.toordinal(datetime.date(Yorig,1,1))
        time = time + datetime.date.toordinal(datetime.date(1900,1,1))
        
        # Changes some names
        if vname=='u10':
            vname='u10m'

        if vname=='v10':
            vname='v10m'
            
        # Create and write output netcdf file
        fname_out = outdir+vname.upper()+'_'+fdate.replace('-','')+'.nc'
        nw = netcdf(fname_out,mode='w',format='NETCDF4')

        dimlon  = nw.createDimension('lon',  len(lon))
        dimlat  = nw.createDimension('lat',  len(lat))
        dimtime = nw.createDimension('time', None)

        varlon = nw.createVariable('lon'  , np.float64,('lon',))
        varlat = nw.createVariable('lat'  , np.float64,('lat',))
        vartime = nw.createVariable('time', np.float64,('time',))
        vardata = nw.createVariable(vname.upper(), 'f4',('time','lat','lon'))
        varlon.long_name = 'longitude of RHO-points'
        varlat.long_name = 'latitude of RHO-points'
        vartime.long_name = 'Time'
        varlon.units = 'degree_east'
        varlat.units = 'degree_north'
        vartime.units = 'days since 1950-01-01 00:00:00'
        vardata.missing_value = 9999.
        vardata.units = coeffs['units'][k]
        vardata.long_name = vlong
	
        varlon[:]=lon
        varlat[:]=lat
        vartime[:]=time
        vardata[:]=data
        nw.close()
    return 



def make_hindcast_ERA5(matlabexec='matlab -nodisplay -nosplash -nodesktop',
                       execdirectory=RUN_dir):
    """
    This function is responsable of running matlab crocotools routine
    make_hindcast_ERA5.m. The functions changes directory to the
    execution dir, then run the start.m file, and run the crocotools 
    routine as usual. For this function to work is important to have
    a consistent crocotools_param.m file and configuration
    (crocotools paths, variables, etc).

    Args:
        matlabexec (str, optional):
            path to matlab executable with respective flags.
            Defaults to 'matlab -nodisplay -nosplash -nodesktop'.
        execdirectory (str, optional): 
            Directory where matlab is executed.
            Defaults to Run_dir (see global variable above).
    """
    os.chdir(execdirectory)
    routines     = ['start', 'make_hindcast_ERA5']
    runarguments = 'clear all; '+'; '.join(routines)+'; exit;'
    command      = matlabexec+' -r "{}"'.format(runarguments)
    print(command)
    os.system(command)
    os.chdir(maindir)
    return 



def add_itolap_bulks(date, itolap=itolap_era5, bulkfreq=1):
    """
    This function grabs a croco blk file of a given date and add
    overlaps based on the previous and following days files. 
    If those files doesnt exist, the function uses the first/last
    record as the overlap value.
    
    Args:
        date (datetime): target date of file to fix
        itolap (int, optional): overlap records. Defaults to itolap_era5.
        bulkfreq (int, optional): era5 frequency in hours. Defaults to 1.
    """
    # Define previous and following dates
    previous   = date-pd.Timedelta(days=1)
    following  = date+pd.Timedelta(days=1)
    
    # Build file names
    fnamep = CROCO_files_dir+fprefix+'_blk_'+previous.strftime('%Y%m%d')+'.nc'
    fname  = CROCO_files_dir+fprefix+'_blk_'+date.strftime('%Y%m%d')+'.nc'
    fnamef = CROCO_files_dir+fprefix+'_blk_'+following.strftime('%Y%m%d')+'.nc'
    print('\n','Adding overlap to ',fname)
    # Load croco bulk associated with given date
    try:
        data    = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                                decode_coords=False, decode_timedelta=False,
                                use_cftime=False)
    except:
        print('\t',fname,' not found, doing nothing. !!')
        return
    # Define time frequency in bulk files
    dt      = bulkfreq/24
    # Build target previous and following timestamps
    ptimes  = np.array([data.bulk_time[0].item()-(n+1)*dt
                        for n in range(itolap)])
    ftimes  = np.array([data.bulk_time[0].item()-(n+1)*dt for n in range(itolap)])
    
    # If previous day data exists open and concat, else fill backwards with first record
    if os.path.isfile(fnamep):
        print('\t','Previous day file found: ',fnamep)
        pdata = xr.open_dataset(fnamep, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)
        pdata = pdata.isel(bulk_time=slice(-(itolap+1),-1))
        data  = xr.concat([pdata,data],'bulk_time')
    else:
        print('\t','Previous day file not found: filling backwards with the first record')
        ntime = np.concatenate([ptimes,data.bulk_time.values])
        data  = data.reindex({'bulk_time':ntime}).bfill('bulk_time')
        
    # If following day data exists open and concat, else forward fill with latest record
    if os.path.isfile(fnamef):
        print('\t','Following day file found: ',fnamef)
        fdata = xr.open_dataset(fnamef, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)
        fdata = fdata.isel(bulk_time=slice(0,itolap))
        data  = xr.concat([data,fdata],'bulk_time')
    else:
        print('\t','Following day file not found: filling forwards with the last record')
        ntime = np.concatenate([data.bulk_time.values, ftimes])
        data  = data.reindex({'bulk_time':ntime}).ffill('bulk_time')
    
    data.to_netcdf(fname+'.2', unlimited_dims=['bulk_time'])
    for var in data.keys():
        os.system('ncatted -a _FillValue,'+var+',d,, '+fname+'.2')
    os.system('ncatted -a _FillValue,bulk_time,d,, '+fname+'.2')
    return 
    
    

# ------------------------------- RUN ROUTINES ------------------------------- #
if __name__=='__main__':
    starttime = datetime.datetime.utcnow()
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Performing ERA5 data conversion, please wait...                   ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    for date in dates:
        fname = pathERA5raw+date.strftime('%F')+'.nc'
        print('Transforming to croco_tools format:',fname)
        blkname = CROCO_files_dir+fprefix+'_blk_'+date.strftime('%Y%m%d')+'.nc'
        if os.path.isfile(blkname):
            print('\t',blkname,' already exists!')
        else:
            ERA5_convert(fname,date)
        print('Done\n')
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Running crocotools make_hindcast_ERA5.m, please wait...           ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    make_hindcast_ERA5()
    print('\nCleaning scratch directory...','                                 ')
    os.system('rm -rf '+scratchdir+'/*.nc')
    print('-------------------------------------------------------------------')
    print('',datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'       ')
    print(' Adding overlap days to croco_blk files, please wait...            ')
    print(' Dates =',dates,'                                                  ')
    print('-------------------------------------------------------------------')
    for date in dates:
        add_itolap_bulks(date)
    print('Overwriting files...','                                            ')
    for date in dates:
        blkname = CROCO_files_dir+fprefix+'_blk_'+date.strftime('%Y%m%d')+'.nc'
        os.remove(blkname)
        os.rename(blkname+'.2',blkname)
            
    print('-------------------------------------------------------------------')
    print(datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),'          ')
    endtime = datetime.datetime.utcnow()
    print('Execution time:',endtime-starttime)
    print('All good','                                                        ')
