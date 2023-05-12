%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create and fill CROCO clim and bry files with OGCM data.
% for a forecast run
%
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
%  Copyright (c) 2006 by Pierrick Penven
%  e-mail:Pierrick.Penven@ird.fr
%
%  Updated    8-Sep-2006 by Pierrick Penven
%  Updated   20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Updated   12-Feb-2016 by P. Marchesiello
%  Updated   14-Oct-2020 by P. Marchesiello
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
tic
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%

%
% Common parameters
%
crocotools_param


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Transforming raw mercator download to a more crocotools compatible format...')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

vars = {'zos' ...
    'uo' ...
    'vo' ...
    'thetao' ...
    'so'};

now = datenum(datestr(datenum(datetime('now', 'TimeZone','Z')),'yyyy-mm-dd'));
raw_mercator_name=[SCRATCH_dir,'/',datestr(now,'yyyy-mm-dd'),'.nc'];
mercator_name=[SCRATCH_dir,'/mercator_',datestr(now,'yyyymmdd'),'.cdf'];
disp(['    Raw motu download: ',raw_mercator_name])
write_mercator_frcst(SCRATCH_dir,'',raw_mercator_name, ...
                      mercator_type,vars,now,Yorig,mercator_name); % write data
disp(' ')


%
%------------------------------------------------------------------------------------
%
% Get the OGCM grid 
% 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([' Get OGCM and model grids...'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
nc=netcdf([SCRATCH_dir,'/mercator_',datestr(now,'yyyymmdd'),'.cdf']);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=-nc{'depth'}(:);
NZ=length(Z);
NZ=NZ-rmdepth;
Z=Z(1:NZ);
close(nc)
disp(' ')


if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
rmask=nc{'mask_rho'}(:);
close(nc)


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['INTERPOLATION STEP'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Desired day: ',datestr(now,'yyyymmdd')])


bryname=[bry_prefix,datestr(now,'yyyymmdd'),nc_suffix];
clmname=[clm_prefix,datestr(now,'yyyymmdd'),nc_suffix];
mercator_name=[SCRATCH_dir,'/mercator_',datestr(now,'yyyymmdd'),'.cdf'];
%---------------------------------------------------------------
% Get time array 
%---------------------------------------------------------------
disp(' ')
disp(['Processing date: ',datestr(now,'yyyymmdd')])
disp(' ')

nc=netcdf(mercator_name);
OGCM_time=nc{'time'}(:);
time_cycle=0;
delta=1; % >1 if subsampling
trange=[1:delta:length(OGCM_time)];
time=zeros(length(trange),1);
for i=1:length(trange)
  time(i)=OGCM_time(trange(i));
end

if makeclim==1 | makebry==1
  if makebry==1
    create_bryfile(bryname,grdname,CROCO_title,obc,...
                  theta_s,theta_b,hc,N,...
                  time,time_cycle,'clobber',vtransform);
    nc_bry=netcdf(bryname,'write');
  else
    nc_bry=[];
  end
  if makeclim==1
    clmname=[clm_prefix,datestr(t,'yyyymmdd'),nc_suffix];
    create_climfile(clmname,grdname,CROCO_title,...
                    theta_s,theta_b,hc,N,...
                    time,time_cycle,'clobber',vtransform);
    nc_clm=netcdf(clmname,'write');
  else
    nc_clm=[];
  end
end

%---------------------------------------------------------------
% Perform interpolations for all selected records
%---------------------------------------------------------------
for tndx=1:length(time)
  disp([' Time step : ',num2str(tndx),' of ',num2str(length(time)),' :'])
  interp_OGCM_frcst(mercator_name,Roa,interp_method,...
                    lonU,latU,lonV,latV,lonT,latT,Z,trange(tndx),...
                    nc_clm,nc_bry,lon,lat,angle,h,pm,pn,rmask, ...
                    tndx,vtransform,obc)
end
%
% Close CROCO files
%
if ~isempty(nc_clm)
  close(nc_clm);
end
if ~isempty(nc_bry)
  close(nc_bry);
end
close(nc); 
disp(' ')

toc
