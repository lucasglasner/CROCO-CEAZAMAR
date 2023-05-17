%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% crocotools_param: common parameter file for the preprocessing
%                  of CROCO simulations using CROCOTOOLS
%
%                  This file is used by make_grid.m, make_forcing.m, 
%                  make_clim.m, make_biol.m, make_bry.m, make_tides.m,
%                  make_NCEP.m, make_OGCM.m, make_...
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Patrick Marchesiello and Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    6-Sep-2006 by Pierrick Penven
%  Updated    2006/10/05 by Pierrick Penven  (add tidegauge observations)
%  Updated    24-Oct-2006 by Pierrick Penven (diagnostics, chla etc...)
%  Updated    08-Apr-2009 by Gildas Cambon
%  Updated    23-Oct-2009 by Gildas Cambon
%  Updated    17-Nov-2011 by Pierrick Penven (CFSR)
%  Updated    07-Nov-2012 by Patrick Marchesiello (cleaning)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1  - Configuration parameters
%      used by make_grid.m (and others..)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isoctave=exist('octave_config_info');
%
%  CROCO title names and directories
%
CROCO_title  = 'CROCO CEAZA FORECAST';
CROCO_config = 'crococeazaf';
%
% Grid dimensions:
%
lonmin =  285.0;   % Minimum longitude [degree east]
lonmax =  289.2;   % Maximum longitude [degree east]
latmin =  -33.5;   % Minimum latitude  [degree north]
latmax =  -27.5;   % Maximum latitude  [degree north]
%
% Grid resolution [degree]
%
dl = 1/36;
%
% Number of vertical Levels (! should be the same in param.h !)
%
N = 50;
%
%  Vertical grid parameters (! should be the same in croco.in !)
%
theta_s    =  7.;
theta_b    =  2.;
hc         = 300.;
vtransform =  2.; % s-coordinate type (1: old- ; 2: new- coordinates)
                  % ! take care to define NEW_S_COORD cpp-key in cppdefs.h 
%
% Topography: choice of filter
%
topo_smooth =  2; % 1: old ; 2: new filter (better but slower)
%
% Minimum depth at the shore [m] (depends on the resolution,
% rule of thumb: dl=1, hmin=300, dl=1/4, hmin=150, ...)
% This affect the filtering since it works on grad(h)/h.
%
hmin = 15;
%
% Maximum depth at the shore [m] (to prevent the generation
% of too big walls along the coast)
%
hmax_coast = 25;
%
% Maximum depth [m] (cut the topography to prevent
% extrapolations below WOA data)
%
hmax = 7500;
%
% Slope parameter (r=grad(h)/h) maximum value for topography smoothing
%
rtarget = 0.2;
%
% Number of pass of a selective filter to reduce the isolated
% seamounts on the deep ocean.
%
n_filter_deep_topo=4;
%
% Number of pass of a single hanning filter at the end of the
% smooting procedure to ensure that there is no 2DX noise in the 
% topography.
%
n_filter_final=2;
%
%  GSHSS user defined coastline (see m_map) 
%  XXX_f.mat    Full resolution data
%  XXX_h.mat    High resolution data
%  XXX_i.mat    Intermediate resolution data
%  XXX_l.mat    Low resolution data
%  XXX_c.mat    Crude resolution data
%
coastfileplot = 'coastline_f.mat';
coastfilemask = 'coastline_f_mask.mat';
%
% Objective analysis decorrelation scale [m]
% (if Roa=0: nearest extrapolation method; crude but much cheaper)
%
%Roa=300e3;
Roa=0;
%
interp_method = 'spline'; % Interpolation method: 'linear' or 'spline'
%
makeplot     = 0;         % 1: create graphics after each preprocessing step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2 - Generic file and directory names 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  CROCOTOOLS directory
%
CROCOTOOLS_dir = '/ceaza/lucas/CROCO-CEAZAMAR/croco_tools/';
%
%  Run directory
%
RUN_dir='/ceaza/lucas/CROCO-CEAZAMAR/FORECAST/';
%
%  CROCO input netcdf files directory
%
CROCO_files_dir=[RUN_dir,'CROCO_FILES/'];
%
%  Global data directory (etopo, coads, datasets download from ftp, etc..)
%
DATADIR='/ceaza/lucas/CROCO-CEAZAMAR/DATA/'; 
%
%  Forcing data directory (ncep, quikscat, datasets download with opendap, etc..)
%
FORC_DATA_DIR = [RUN_dir,'DATA/'];
%
if (isoctave == 0)
	eval(['!mkdir ',CROCO_files_dir])
else
	system(['mkdir ',CROCO_files_dir])
end
%
% CROCO file names (grid, forcing, bulk, climatology, initial)
%
grdname  = [CROCO_files_dir,'crococeazaf_grd.nc'];
frcname  = [CROCO_files_dir,'crococeazaf_frc.nc'];
blkname  = [CROCO_files_dir,'crococeazaf_blk.nc'];
clmname  = [CROCO_files_dir,'crococeazaf_clm.nc'];

bryname  = [CROCO_files_dir,'crococeazaf_bry.nc'];
ininame  = [CROCO_files_dir,'crococeazaf_ini.nc'];
bioname  = [CROCO_files_dir,'crococeazaf_frcbio.nc']; % Iron Dust forcing for PISCES
rivname =  [CROCO_files_dir,'crococeazaf_runoff.nc'];
%
% intermediate z-level data files (not used in simulations)
%
oaname   = [CROCO_files_dir,'crococeazaf_oa.nc'];    % for climatology data processing
Zbryname = [CROCO_files_dir,'crococeazaf_bry_Z.nc']; % for boundary data processing
%
% Generic forcing file root names for interannual simulations (NCEP/GFS)
%
frc_prefix=[CROCO_files_dir,'crococeazaf_frc'];      % forcing file name 
blk_prefix=[CROCO_files_dir,'crococeazaf_blk'];      % bulk file name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Topography netcdf file name (ETOPO 2 or any other netcdf file
%  in the same format)
%
topofile = [DATADIR,'Topo/etopo2.nc'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3 - Surface forcing parameters
%     used by make_forcing.m and by make_bulk.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COADS directory (for climatology runs)
%
coads_dir=[DATADIR,'COADS05/'];
%
% COADS time (for climatology runs)
%
coads_time=(15:30:345); % days: middle of each month
coads_cycle=360;        % repetition of a typical year of 360 days  
%
%coads_time=(15.2188:30.4375:350.0313); % year of 365.25 days in case
%coads_cycle=365.25;                    % interannual QSCAT winds  
%                                       % are used with clim. heat flux
%
% Pathfinder SST data used by pathfinder_sst.m
%
pathfinder_sst_name=[DATADIR,...
                    'SST_pathfinder/climato_pathfinder.nc'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 4 - Open boundaries and initial conditions parameters
%     used by make_clim.m, make_biol.m, make_bry.m
%             make_OGCM.m and make_OGCM_frcst.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Open boundaries switches (! should be consistent with cppdefs.h !)
%
obc = [1 0 1 1]; % open boundaries (1=open , [S E N W])
%
%  Level of reference for geostrophy calculation
%
zref = -1000;
%
%  initial/boundary data options (1 = process)
%  (used in make_clim, make_biol, make_bry,
%   make_OGCM.m and make_OGCM_frcst.m)
%
makeini    = 1;   % initial data
makeclim   = 0;   % climatological data (for boundaries and nudging layers)
makebry    = 1;   % lateral boundary data
makenpzd   = 0;   % initial and boundary data for NChlPZD and N2ChlPZD2 models
makebioebus= 0;   % initial and boundary data for BioEBUS model
makepisces = 0;   % initial and boundary data for PISCES model
%
%
makeoa     = 1;   % oa data (intermediate file)
makeZbry   = 1;   % boundary data in Z coordinate (intermediate file)
insitu2pot = 1;   % transform in-situ temperature to potential temperature
%
%  Day of initialisation for climatology experiments (=0 : 1st january 0h)
%
tini=0;  
%
% Select Climatology Atlas (temp, salt and biological variables) from:
%    - World Ocean Atlas directory (WOA2009)  OR ...
%    - CARS2009 climatology directory (CARS2009)
%
woa_dir       = [DATADIR,'WOA2009/'];
cars2009_dir  = [DATADIR,'CARS2009/'];
climato_dir   = woa_dir;
%
% Pisces biogeochemical seasonal climatology
%
woapisces_dir = [DATADIR,'WOAPISCES/'];  % only compatible with woa_dir
%
% Surface chlorophyll seasonal climatology (SeaWifs)
%
chla_dir=[DATADIR,'SeaWifs/'];
%
% Runoff monthly seasonal climatology (Dai and Trenberth)
%
global_clim_riverdir=[DATADIR,'RUNOFF_DAI/'];
global_clim_rivername=[global_clim_riverdir,'Dai_Trenberth_runoff_global_clim.nc'];
%
%  Set times and cycles for the boundary conditions: 
%   monthly climatology 
%
woa_time=(15:30:345); % days: middle of each month
woa_cycle=360;        % repetition of a typical year of 360 days  
%
%woa_time=(15.2188:30.4375:350.0313); % year of 365.25 days in case
%woa_cycle=365.25;                    % interannual QSCAT winds are used 
%                                     % with clim. boundary conditions
%
%   For rivers setup : go in the routine Rivers/make_runoff.m to
%   setup your options
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 5 - Parameters for tidal forcing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TPXO file name (TPXO6 or TPXO7)
%
tidename=[DATADIR,'TPXO7/TPXO7.nc'];
%
% Self-Attraction and Loading GOT99.2 file name
%
sal_tides=1;
salname=[DATADIR,'GOT99.2/GOT99_SAL.nc'];
%
% Number of tides component to process
%
Ntides=10;
%
% Chose order from the rank in the TPXO file :
% "M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm"
% " 1  2  3  4  5  6  7  8  9 10"
%
tidalrank=[1 2 3 4 5 6 7 8 9 10];
%
% Compare with tidegauge observations
%
lon0 =  18.37;   % Example: 
lat0 = -33.91;   % Cape Town location
Z0   =  1;       % Mean depth of tide gauge
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 6 - Reference date and simulation times
%     (used for make_tides, make_CFSR (or make_NCEP), make_OGCM)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Yorig         = 1950;          % reference time for vector time
                               % in croco initial and forcing files
%
Ymin          = 2005;          % first forcing year
Ymax          = 2005;          % last  forcing year
Mmin          = 1;             % first forcing month
Mmax          = 3;             % last  forcing month
%
Dmin          = 1;             % Day of initialization
Hmin          = 0;             % Hour of initialization
Min_min       = 0;             % Minute of initialization
Smin          = 0;             % Second of initialization
%
SPIN_Long     = 0;             % SPIN-UP duration in Years
%
Mth_format    = '%02d';        % Number of digit for month on input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 7 - Parameters for Interannual forcing (SODA, ECCO, CFSR, NCEP, ...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Download_data = 1;   % Get data from OPENDAP sites  
level         = 0;   % AGRIF level; 0 = parent grid
%					  
NCEP_version  = 3;   % NCEP version: 
                     % [ CFSR up-to-date product are recommandated ]
                     %  1: NCEP/NCAR Reanalysis, 1/1/1948 - present
                     %  2: NCEP-DOE Reanalysis, 1/1/1979 - present
                     %  3: CFSR (Climate Forecast System Reanalysis), 
                     %           1/1/1979 - 31/3/2011
%
%--------------------------------------------
% Options for make_NCEP and make_QSCAT_daily
%--------------------------------------------
%
% NCEP data directory for files downloaded via opendap
%
if NCEP_version  == 1;
  NCEP_dir= [FORC_DATA_DIR,'NCEP1_',CROCO_config,'/']; 
elseif NCEP_version  == 2;
  NCEP_dir= [FORC_DATA_DIR,'NCEP2_',CROCO_config,'/'];
elseif NCEP_version  == 3;
  NCEP_dir= [FORC_DATA_DIR,'CFSR_',CROCO_config,'/']; % CFSR data dir. [croco format]
end
makefrc      = 0;       % 1: create forcing files
makeblk      = 1;       % 1: create bulk files
QSCAT_blk    = 0;       % 1: a) correct NCEP frc/bulk files with
                        %        u,v,wspd fields from daily QSCAT data
                        %    b) download u,v,wspd in QSCAT frc file
add_tides    = 0;       % 1: add tides
%
% Overlap parameters 
%
itolap_qscat = 0;       % 2 records for daily  QSCAT
itolap_ncep  = 0;       % 8 records for 4-daily NCEP
%
%--------------------------------------------------
% Options for make_QSCAT_daily and make_QSCAT_clim   
%--------------------------------------------------
%
QSCAT_dir        = [FORC_DATA_DIR,'QSCAT_',CROCO_config,'/']; % QSCAT data dir. [croco format]
QSCAT_frc_prefix = [frc_prefix,'_QSCAT_'];                    % QSCAT Generic file name for
                                                              % interannual simulations
QSCAT_clim_file  = [DATADIR,'QuikSCAT_clim/',...              % QSCAT climatology file
                    'roms_SCOW_month_clim_1999_2009.nc'];    % for make_QSCAT_clim.
%
%--------------------------------------------------
%  Options for make_ECMWF and make_ECMWF_daily  
%--------------------------------------------------
%
Reformat_ECMWF = 1;
ECMWF_dir= [FORC_DATA_DIR,'ECMWF_',CROCO_config,'/'];  % ERA-I data dir. [croco format]
My_ECMWF_dir=[FORC_DATA_DIR,'ERAI/'];                  % ERA-I native data downloaded 
                                                       % with python script
itolap_ecmwf = 0;                                      % 3 records for daily  ECMWF
%
%--------------------------------------------------
%  Options for make_ERA5 and make_FORECAST_ERA5
%--------------------------------------------------
%
ERA5_dir    = [RUN_dir,'/SCRATCH/'];          % ERA5 data dir. [croco format]
My_ERA5_dir = [FORC_DATA_DIR,'/ERA5/'];       % ERA5 native data downloaded with python script
itolap_era5 = 0;                              % 2 records = 2 hours
ERA5_delay  = 6;                              % Delay days of ERA5 NRT product
ERA5_offset = 10;                             % Days of ERA5 NRT to compute from the latest

%
%
%--------------------------------------------
% Options for make_OGCM or make_OGCM_mercator
%--------------------------------------------
%
OGCM        = 'mercator';        % Select OGCM: SODA, ECCO, mercator
%
OGCM_dir    = [FORC_DATA_DIR,OGCM,'_',CROCO_config,'/'];  % OGCM data dir. [croco format]
%
bry_prefix  = [CROCO_files_dir,'crococeazaf_bry_'];    % generic boundary file name
clm_prefix  = [CROCO_files_dir,'crococeazaf_clm_'];    % generic climatology file name
ini_prefix  = [CROCO_files_dir,'crococeazaf_ini_'];    % generic initial file name
OGCM_prefix = [OGCM,'_'];                                 % generic OGCM file name 

if strcmp(OGCM,'mercator')
    % For GLORYS 12 reanalysis extraction + download using python motuclient
    % ========================
    motu_url_reana='http://my.cmems-du.eu/motu-web/Motu';
    service_id_reana='GLOBAL_MULTIYEAR_PHY_001_030-TDS';
    product_id_reana='cmems_mod_glo_phy_my_0.083_P1D-m';
end
%
% Number of OGCM bottom levels to remove 
% (usefull if CROCO depth is shallower than OGCM depth)
%
rmdepth     = 2;
%
% Overlap parameters : nb of records around each monthly sequence
%
itolap_0    = 1;   % before
itolap_0    = 1;   % after
                   %
%--------------------------
% Options for make_bry_WKB 
%--------------------------
%
wkb_prefix=[CROCO_files_dir,'crococeazaf_wkb'];
wkb_obc= [1 1 1 1];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 8 - Parameters for the forecast system
%
%     --> select OGCM name above (mercator ...)
%     --> don't forget to define in cppdefs.h:
%                    - ROBUST_DIAG
%                    - CLIMATOLOGY
%                    - BULK_FLUX
%                    - TIDES if you choose so, but without TIDERAMP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
FRCST_dir = [FORC_DATA_DIR];  % path to local OGCM data directory
%
% Number of FORECAST/forecast days
%
if strcmp(OGCM,'ECCO')
  hdays=6;
  fdays=10;
elseif strcmp(OGCM,'mercator')
  hdays=6;
  fdays=11;
end
%
% Local time= UTC + timezone
%
timezone = 0;
%
% Add tides
%
add_tides_fcst = 1;       % 1: add tides
%
%  MERCATOR case: 
%  =============
%  To download data: set login/password (http://marine.copernicus.eu)
%  and path to croco's motuclient python package;
%  or set pathMotu='' (empty) to use your own motuclient
%
%  Various sets of data are proposed in the 
%  Copernicus web site (Mercator, UK Met Office ...)
%
if strcmp(OGCM,'mercator')
  SCRATCH_dir     = [RUN_dir,'/SCRATCH/'];
  user     = 'XXX';
  password = 'XXX';

  pathMotu =[CROCOTOOLS_dir,'Forecast_tools/'];

  mercator_type=1;   % 1 -->  1/12 deg Mercator forecast
                     % 2 -->  1/4  deg Met-Office forecast (GloSea5)
  if mercator_type==1
      motu_url_fcst='http://nrt.cmems-du.eu/motu-web/Motu';
      service_id_fcst='GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS';
      product_id_fcst='global-analysis-forecast-phy-001-024';
      
  elseif mercator_type==2
      motu_url_fcst='http://nrt.cmems-du.eu/motu-web/Motu';
      service_id_fcst='GLOBAL_ANALYSISFORECAST_PHY_CPL_001_015-TDS';
      product_id_fcst='MetO-GLO-PHY-CPL-dm-TEM'
  end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 9 Parameters for the diagnostic tools
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
DIAG_dir = [CROCOTOOLS_dir,'Diagnostic_tools/'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











