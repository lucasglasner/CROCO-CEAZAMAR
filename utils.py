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

def find_itolap_datesnfiles(date, inputfiledir, fprefix,
                            toffset=1):
    """
    Simple function to find the file names of the current date,
    previous and following.
    """
    # Define previous and following dates
    previous   = date-pd.Timedelta(days=toffset)
    following  = date+pd.Timedelta(days=toffset)
    
    # Build file names
    fnamep = inputfiledir+fprefix+previous.strftime('%Y%m%d')+'.nc'
    fname  = inputfiledir+fprefix+date.strftime('%Y%m%d')+'.nc'
    fnamef = inputfiledir+fprefix+following.strftime('%Y%m%d')+'.nc'
    return fnamep, fname, fnamef

def add_previous_itolap(pdata,data,itolap,timename):
    """
    Simple function to add previous overlaps based on 
    previous and current data
    """
    pdata = pdata.isel({timename:slice(-(itolap),None)})
    ndata  = xr.concat([pdata,data],timename)
    return ndata

def add_following_itolap(fdata,data,itolap,timename):
    """
    Simple function to add previous overlaps based on 
    following and current data
    """
    fdata = fdata.isel({timename:slice(None, itolap)})
    ndata  = xr.concat([data,fdata],timename)
    return ndata

def add_itolap_blk(date, itolap, variables, inputfiledir, outputfiledir,
                   fprefix='croco_blk_', timename='bulk_time',
                   freq=1):
    """
    This function grabs a croco blk file of a given date
    and add overlaps based on the previous and following days files. 
    If those files doesnt exist, the function uses the first/last
    record as the overlap value.
    
    Args:
        date (datetime): target date of file to fix
        variables (list): list of strings with the variables
        where to apply the overlap.
        itolap (int): number of overlap records.
        inputfiledir (str): input file directory
        outputfiledir (str): output file directory
        fprefix (str, optional): croco file prefix.
        Defaults to 'croco_blk_'
        freq (int, optional): data frequency in hours.
        Defaults to 1.
        timename (str, optional): dataset time dimension name.
        Defaults to 'bulk_time'
    """
    fnamep, fname, fnamef = find_itolap_datesnfiles(date,
                                                    inputfiledir,
                                                    fprefix)
    # Load croco forcing associated with given dates
    try:
        data    = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                                decode_coords=False, decode_timedelta=False,
                                use_cftime=False)
        if len(data[timename])==24/freq+itolap*2:
            print(fname,' already has overlap times !!')
            return
    except Exception as e:
        print('\t',fname,e)
        return
    keys       = list(data.keys())
    complement = np.array([k if k not in variables else None for k in keys])
    complement = complement[complement!=None]
    complement = data[complement]
    data       = data[variables]
    
    print('\n','Adding overlap to ',fname)
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
        print('\t','Previous day file found')
        pdata = xr.open_dataset(fnamep, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)[variables]
        if len(pdata[timename])==24/freq+itolap*2:
            pdata = pdata.isel({timename:slice(itolap,-itolap)})
        data = add_previous_itolap(pdata, data, itolap, timename)
    else:
        print('\t','Previous day file not found:',
            'filling backwards with the first record')
        ntime = np.hstack([ptimes,data[timename].values])
        data  = data.reindex({timename:ntime}).bfill(timename)
    
    # If following day data exists open and concat,
    # else forward fill with latest record
    if os.path.isfile(fnamef):
        print('\t','Following day file found')
        fdata = xr.open_dataset(fnamef, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)[variables]
        if len(fdata[timename])==24/freq+itolap*2:
            fdata = fdata.isel({timename:slice(itolap,-itolap)})
        data = add_following_itolap(fdata, data, itolap, timename)
    else:
        print('\t','Following day file not found:',
            'filling forwards with the last record')
        ntime = np.hstack([data[timename].values, ftimes])
        data  = data.reindex({timename:ntime}).ffill(timename)
    
    ofname = fname.replace(inputfiledir,outputfiledir)
    xr.merge([complement, data]).to_netcdf(ofname, mode='w', engine='netcdf4',
                   unlimited_dims=[timename])  
    return


def croco_bry_swapdims(data):
    """
    Simple function to remove repeated time variables
    and change everything to "bry_time"
    """
    try:
        timevars = ['temp_time','salt_time','zeta_time',
                    'v3d_time','v2d_time']
        swapdict = {keys:'bry_time' for keys in timevars}
        data = data.reset_index(list(timevars), drop=True)
        data = data.rename(swapdict)
    except:
        pass
    return data

def add_itolap_bry(date, itolap, variables, inputfiledir, outputfiledir,
                   fprefix='croco_bry_', timename='bry_time',
                   freq=6):
    """
    This function grabs a croco bry file of a given date
    and add overlaps based on the previous and following days files. 
    If those files doesnt exist, the function uses the first/last
    record as the overlap value.
    
    Args:
        date (datetime): target date of file to fix
        variables (list): list of strings with the variables
        where to apply the overlap.
        itolap (int): number of overlap records.
        inputfiledir (str): input file directory
        outputfiledir (str): output file directory
        fprefix (str, optional): croco file prefix.
        Defaults to 'croco_bry_'
        freq (int, optional): data frequency in hours.
        Defaults to 6.
        timename (str, optional): dataset time dimension name.
        Defaults to 'bry_time'
    """
    fnamep, fname, fnamef = find_itolap_datesnfiles(date, inputfiledir,fprefix)
    print('\n','Adding overlap to ',fname)
    # Load croco bry associated with given date
    try:
        data = xr.open_dataset(fname, decode_times=False, decode_cf=False,
                            decode_coords=False, decode_timedelta=False,
                            use_cftime=False)
        if len(data[timename])==24/freq+itolap*2:
            print(fname,' already has overlap times !!')
            return
        data = croco_bry_swapdims(data)
    except Exception as e:
        print('\t',fname,e)
        return
    keys       = list(data.keys())
    complement = np.array([k if k not in variables else None for k in keys])
    complement = complement[complement!=None]
    complement = data[complement]
    data       = data[variables]
    
    # Define time frequency in bulk files
    dt      = freq/24
    # Build target previous and following timestamps
    ptimes  = sorted(np.array([data[timename][0].item()-(n+1)*dt
                        for n in range(itolap)]))
    ftimes  = sorted(np.array([data[timename][-1].item()+(n+1)*dt
                        for n in range(itolap)]))
    #print(ptimes,data[timename].values,ftimes)
    # If previous day data exists open and concat,
    # else fill backwards with first record
    if os.path.isfile(fnamep):
        print('\t','Previous day file found')
        pdata = xr.open_dataset(fnamep, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)[variables]
        if len(pdata[timename])==24/freq+itolap*2:
            pdata = pdata.isel({timename:slice(itolap,-itolap)})
        pdata = croco_bry_swapdims(pdata).sortby(timename)
        data  = add_previous_itolap(pdata, data, itolap, timename)
    else:
        print('\t','Previous day file not found:',
            'filling backwards with the first record')
        ntime = np.hstack([ptimes,data[timename].values])
        data  = data.reindex({timename:ntime}).sortby(timename).bfill(timename)
    
    # If following day data exists open and concat,
    # else forward fill with latest record
    if os.path.isfile(fnamef):
        print('\t','Following day file found')
        fdata = xr.open_dataset(fnamef, decode_times=False, decode_cf=False,
                              decode_coords=False, decode_timedelta=False,
                              use_cftime=False)[variables]
        fdata = croco_bry_swapdims(fdata).sortby(timename)
        if len(fdata[timename])==24/freq+itolap*2:
            fdata = fdata.isel({timename:slice(itolap,-itolap)})
        data  = add_following_itolap(fdata, data, itolap, timename)
    else:
        print('\t','Following day file not found:',
            'filling forwards with the last record')
        ntime = np.hstack([data[timename].values, ftimes])
        data  = data.reindex({timename:ntime}).sortby(timename).ffill(timename)
        
    ofname = fname.replace(inputfiledir,outputfiledir)
    xr.merge([complement, data]).to_netcdf(ofname, mode='w', engine='netcdf4',
                   unlimited_dims=[timename])  
    return
