'''
 # @ Author: Your name
 # @ Create Time: 2023-04-27 14:30:37
 # @ Modified by: Your name
 # @ Modified time: 2023-04-27 14:42:56
 # @ Description:
 
 Utility functions for forecast/hindcast preprocessing
 '''
 
 
# ---------------------------------- IMPORTS --------------------------------- #
import os
import xarray as xr
import numpy as np
import pandas as pd
from glob import glob
 

# --------------------------------- FUNCTIONS -------------------------------- #
def cleandirectory(directory, filetypes=['*']):
    """
    This function just removes all files in the target directory
    with the file extension in filetypes (accepts glob like expressions)

    Args:
        directory (str): path to directory
        filetypes (list, optional): files extensions to remove.
        Defaults to ['*'].
    """
    for ftyp in filetypes:
        for f in glob(directory+'/'+ftyp):
            os.remove(f)
    return

def sliceandconcat_crocoblk(fname,data, itolap, timename, freq, timeshift):
    if timeshift=='previous':
        print('\t','Previous day file found: ',fname)
        pdata = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                                decode_coords=False, decode_timedelta=False,
                                use_cftime=False)
        if len(pdata[timename])==24/freq+itolap*2:
            # If pdata already has overlaps just grab the inner data
            pdata = pdata.isel({timename:slice(itolap,-itolap)})
        pdata = pdata.isel({timename:slice(-(itolap+1),-1)})
        ndata  = xr.concat([pdata,data],timename)
    elif timeshift=='following':
        print('\t','Following day file found: ',fname)
        fdata = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)
        if len(fdata[timename])==24/freq+itolap*2:
            # If fdata already has overlaps just grab the inner data
            fdata = fdata.isel({timename:slice(itolap,-itolap)})
        fdata = fdata.isel({timename:slice(0,itolap)})
        ndata  = xr.concat([data,fdata],timename)
    else:
        raise ValueError
    return ndata


def add_itolap(date, itolap,timename,
               inputfiledir, outputfiledir,
               fprefix='croco_blk_', freq=1):
    """
    This function grabs a croco forcing (blk or bry) file of a given date
    and add overlaps based on the previous and following days files. 
    If those files doesnt exist, the function uses the first/last
    record as the overlap value.
    
    Args:
        date (datetime): target date of file to fix
        itolap (int): number of overlap records.
        timename (str): dataset time dimension name.
        inputfiledir (str): input file directory
        outputfiledir (str): output file directory
        fprefix (str, optional): croco file prefix. Defaults to 'croco_blk_'
        freq (int, optional): data frequency in hours. Defaults to 1.
    """
    # Define previous and following dates
    previous   = date-pd.Timedelta(days=1)
    following  = date+pd.Timedelta(days=1)
    
    # Build file names
    fnamep = inputfiledir+fprefix+previous.strftime('%Y%m%d')+'.nc'
    fname  = inputfiledir+fprefix+date.strftime('%Y%m%d')+'.nc'
    fnamef = inputfiledir+fprefix+following.strftime('%Y%m%d')+'.nc'
    print('\n','Adding overlap to ',fname)
    # Load croco bulk associated with given date
    try:
        data    = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                                decode_coords=False, decode_timedelta=False,
                                use_cftime=False)
        if len(data[timename])==24/freq+itolap*2:
            print('\t',fname,' already has overlap times !!')
            return
    except Exception as e:
        print('\t',fname,e)
        return
    # Define time frequency in bulk files
    dt      = freq/24
    # Build target previous and following timestamps
    ptimes  = sorted(np.array([data[timename][0].item()-(n+1)*dt
                        for n in range(itolap)]))
    ftimes  = sorted(np.array([data[timename][-1].item()+(n+1)*dt
                        for n in range(itolap)]))

    # If previous day data exists open and concat,
    # else fill backwards with first record
    if os.path.isfile(fnamep):
        if 'blk' in fprefix:
            data = sliceandconcat_crocoblk(fnamep,data, itolap=itolap,
                                        timename=timename, freq=freq,
                                        timeshift='previous')
    else:
        print('\t','Previous day file not found:',
              'filling backwards with the first record')
        ntime = np.concatenate([ptimes,data[timename].values])
        data  = data.reindex({timename:ntime}).bfill(timename)
        
    # If following day data exists open and concat,
    # else forward fill with latest record
    
    if os.path.isfile(fnamef):
        if 'blk' in fprefix:
            data = sliceandconcat_crocoblk(fnamef,data, itolap=itolap,
                                        timename=timename, freq=freq,
                                        timeshift='following')
    else:
        print('\t','Following day file not found:',
              'filling forwards with the last record')
        ntime = np.concatenate([data[timename].values, ftimes])
        data  = data.reindex({timename:ntime}).ffill(timename)
    
    ofname = fname.replace(inputfiledir,outputfiledir)
    data.to_netcdf(ofname, mode='w', engine='netcdf4',
                   unlimited_dims=[timename])
    return 
