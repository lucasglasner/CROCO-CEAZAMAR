#
# For ERA5 python crocotools parameters list
#
# CAUTION IT MUST BE CONSISTENT with your MATLAB CROCOTOOLS_PARAM.m file in Run directory <- so make it consistent? lucas 2022
# *******************************************************************************
#                         U S E R  *  O P T I O N S
# *******************************************************************************
#
# General path
#
import os 
import re

# ---------------- go to find variables to crocotools_param.m ---------------- #
# ------------------------------ edit lucasglasner 2022 ---------------------- #
runfolder = '/ceaza/lucas/CROCO/run'
paramFile = open(runfolder+'/crocotools_param.m', 'r')
for line in paramFile:
    if re.search('makefrc', line):
        # ----------------- download and make frc required variables --------- #
        makefrc = int(line.split(";")[0].split("=")[1])
    if re.search('makeblk', line):
        # ----------------- download and make blk required variables --------- #
        makeblk = int(line.split(";")[0].split("=")[1])
    if re.search('RUN_dir=', line):
        # ----------------------------- config directory --------------------- #
        config_dir = line.split(";")[0].split("=")[1].replace("'","")
    if re.search('Yorig', line):
        # ---------------------------- year origin of time ------------------- #
        Yorig=int(line.split(";")[0].split("=")[1])
    if re.search('Ymin', line):
        # -------------------------------- year start ------------------------ #
        year_start=int(line.split(";")[0].split("=")[1])
    if re.search('Ymax', line):
        # --------------------------------- last year ------------------------ #
        year_end=int(line.split(";")[0].split("=")[1])
    if re.search('Mmin', line):
        # -------------------------------- month start ----------------------- #
        month_start=int(line.split(";")[0].split("=")[1])
    if re.search('Mmax', line):
        # --------------------------------- month end ------------------------ #
        month_end=int(line.split(";")[0].split("=")[1])
    if re.search('CROCO_config =', line):
        # --------------------------------- month end ------------------------ #
        config_name=line.split(";")[0].split("=")[1].replace("'","").replace(" ","")
# ---------------------------------------------------------------------------- #

# year_start=2012
# month_start=1

#year_start          = 1993; #         % first forcing year
#year_end          = 2020;  #        % last  forcing year
#month_start          = 1;      #       % first forcing month
#month_end          = 12;      #       % last  forcing month
# config_name = 'CHILE_CENTRAL'
#
# Original (raw) ERA5 directory
#
era5_dir_raw = config_dir + 'DATA/ERA5'
#
# Output ERA5 directory
#
era5_dir_processed = config_dir + 'DATA/ERA5_' + config_name



#
#
# Dates limits
#
# year_start = 2003
# month_start = 1
# year_end = 2003
# month_end = 12


# ---------------------------------------------------------------------------- #
#                                 OVERLAP DAYS                                 #
# ---------------------------------------------------------------------------- #
n_overlap = 8

# Request time (daily hours '00/01/.../23')
#
time = '00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23'

# Request area ([north, west, south, east])
#
ownArea = 0 	# 0 if area from a crocotools_param.m file
                # 1 if own area

if ownArea == 0: 
    # To complete if ownArea==0
    paramFile = config_dir + 'crocotools_param.m' # path the crocotools_param file of the simulation
    
else:
    # To complete if ownArea==1
    lonmin=7
    lonmax=23
    latmin=-45
    latmax=-20
#
# Variable names and conversion coefficients  
# TP: convert from accumlated m in a hour into   kg m-2 s-1
#
cff_tp=1000./3600. # m in 1 hour -> kg m-2 s-1
# Heat flux J m-2 in one hour into W m-2
#
cff_heat=1./3600.   # J m-2 in 1 hour -> W m-2
# Stress N m-2 accumulated in hour
#
cff_stress=1./3600.

# Names, conversion coefficients and new units
# ------------------- download blk and frc needed variables ------------------ #
# ------------------------------ edit lucas 2022 ----------------------------- #
variables = [ 'lsm'  , 'tp'        , 'e'         , 'q'      , 'u10'  , 'v10'  , 'ewss'    , 'nsss'    , 't2m', 'ssr'   , 'str'   , 'strd'  , 'slhf'  , 'sshf'  ]
conv_cff  = [ 1.     , cff_tp      , cff_tp      , 1.       , 1.     , 1.     , cff_stress, cff_stress, 1.   , cff_heat, cff_heat, cff_heat, cff_heat, cff_heat]
units     = [ '(0-1)', 'kg m-2 s-1', 'kg m-2 s-1', 'kg kg-1', 'm s-1', 'm s-1', 'N m-2 s' , 'N m-2 s' , 'K'  , 'W m-2' , 'W m-2' , 'W m-2' , 'W m-2' , 'W m-2' ]  

#variables = ['u10','v10','ewss','nsss']
#conv_cff = [1.,1.,cff_stress,cff_stress]
#units=['m s-1','m s-1','N m-2 s','N m-2 s']
# *******************************************************************************
#                         E N D     U S E R  *  O P T I O N S
# *******************************************************************************
