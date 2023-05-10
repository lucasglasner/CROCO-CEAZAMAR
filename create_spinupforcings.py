#!/ceaza/lucas/miniconda3/envs/main/bin/python3
'''
 # @ Author: Your name
 # @ Create Time: 2023-05-10 10:28:16
 # @ Modified by: Your name
 # @ Modified time: 2023-05-10 10:28:30
 # @ Description:
 
 This script just copies the ini, blk, frc and bry files of a single year
 and changes the name and timestamp for building the spin up years forcing.
 
 '''

import os
import pandas as pd
import time

CROCO_files_dir = '/ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/CROCO_FILES/'
TARGETYEAR      = 2022
SPINUPYEARS     = 2
SIMNAME         = 'crococeazah'


def modify_ini(infile,outfile, YEAR, NYEARS):
    # data      = xr.open_dataset(fname, decode_cf=False, decode_coords=False,
    #                             decode_timedelta=False, decode_times=False)
    timedelta = pd.to_datetime(str(YEAR)+'-01-01')
    timedelta = timedelta-pd.to_datetime(str(YEAR-NYEARS)+'-01-01')
    timedelta = timedelta.total_seconds()
    command = "ncap2 -O -s 'scrum_time=scrum_time-{}' {} /ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/SCRATCH/tmp.nc".format(timedelta,infile)
    os.system(command)
    command = "ncap2 -O -s 'ocean_time=ocean_time-{}' /ceaza/lucas/CROCO-CEAZAMAR/HINDCAST/SCRATCH/tmp.nc {}".format(timedelta,outfile)
    os.system(command)
    os.system('rm -f tmp.nc')
    return

def modify_blk(infile,outfile, date, ndate):
    # data      = xr.open_dataset(fname, decode_cf=False, decode_coords=False,
    #                             decode_timedelta=False, decode_times=False)
    timedelta = (date-ndate).days
    command = "ncap2 -O -s 'bulk_time=bulk_time-{}' {} {}".format(timedelta,infile,outfile)
    os.system(command)
    return

def modify_bry(infile,outfile, date, ndate):
    # data      = xr.open_dataset(fname, decode_cf=False, decode_coords=False,
    #                             decode_timedelta=False, decode_times=False)
    timedelta = (date-ndate).days
    command = "ncap2 -O -s 'bry_time=bry_time-{}' {} {}".format(timedelta,infile,outfile)
    os.system(command)
    return

def modify_frc(infile,outfile):
    command = "cp {} {}".format(infile,outfile)
    os.system(command)
    return


if __name__=='__main__':
    ininame = CROCO_files_dir+SIMNAME+'_ini_'+str(TARGETYEAR)+'0101.nc'
    print('Copy',ininame,'file')
    modify_ini(ininame, CROCO_files_dir+'spinup/'+SIMNAME+'_ini_'+str(TARGETYEAR-SPINUPYEARS)+'0101.nc',TARGETYEAR,SPINUPYEARS)
    for yr in range(1,SPINUPYEARS+1):
        YEAR = str(TARGETYEAR-yr)
        for i in range(365):
            date  = pd.to_datetime(str(TARGETYEAR)+'-01-01')+pd.Timedelta(days=i)
            ndate = pd.to_datetime(str(YEAR)+'-01-01')+pd.Timedelta(days=i)
            print(date, ndate)
            
            blkname = CROCO_files_dir+SIMNAME+'_blk_'+date.strftime('%Y%m%d')+'.nc'
            bryname = CROCO_files_dir+SIMNAME+'_bry_'+date.strftime('%Y%m%d')+'.nc'
            frcname = CROCO_files_dir+SIMNAME+'_frc_'+date.strftime('%Y%m%d')+'.nc'
            
            modify_blk(blkname, CROCO_files_dir+'spinup/'+SIMNAME+'_blk_'+ndate.strftime('%Y%m%d')+'.nc', date, ndate)
            # blk.to_netcdf(CROCO_files_dir+'spinup/'+SIMNAME+'_blk_'+ndate.strftime('%Y%m%d')+'.nc')

            modify_bry(bryname, CROCO_files_dir+'spinup/'+SIMNAME+'_bry_'+ndate.strftime('%Y%m%d')+'.nc', date, ndate)
            # bry.to_netcdf(CROCO_files_dir+'spinup/'+SIMNAME+'_bry_'+ndate.strftime('%Y%m%d')+'.nc')
            
            modify_frc(frcname, CROCO_files_dir+'spinup/'+SIMNAME+'_frc_'+ndate.strftime('%Y%m%d')+'.nc')
            # frc.to_netcdf(CROCO_files_dir+'spinup/'+SIMNAME+'_frc_'+ndate.strftime('%Y%m%d')+'.nc')