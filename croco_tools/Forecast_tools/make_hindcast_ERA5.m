%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  make_hindcast_ERA5.m
% 
%  Create and fill frc and bulk files with ERA5 data.
%  (ERA-5 Reanalysis)
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
%  Updated   L. Glasner, D. Donoso, G. Cambon. P. Penven (Apr 2023) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param                      
blk_prefix=[blk_prefix,'_'];
frc_prefix=[frc_prefix,'_'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if level==0
    nc_suffix='.nc';
else
    nc_suffix=['.nc.',num2str(level)];
    grdname=[grdname,'.',num2str(level)];
end

% Get the model grid
disp(' ')
disp(['Read the grid in ',grdname])
nc=netcdf(grdname);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
lonu=nc{'lon_u'}(:);
latu=nc{'lat_u'}(:);
lonv=nc{'lon_v'}(:);
latv=nc{'lat_v'}(:);
angle=nc{'angle'}(:);
close(nc);

% Define latest ERA5 data as today date minus the defined ERA5 delay
% Then try to create bulks for the last 10 days since the latest ERA5 file
now = datenum(datetime('now', 'TimeZone','Z'))-ERA5_delay;

%
%Loop on the years and the months
%
disp(['====================='])
disp(['INTERPOLATION STEP'])
disp(['====================='])
disp(['Loop on the desired days: ',datestr(now-ERA5_offset,'yyyymmdd'),' - ',datestr(now,'yyyymmdd')])
for t=(now-ERA5_offset):now
    blkname=[blk_prefix,datestr(t,'yyyymmdd'),nc_suffix];       
    frcname=[frc_prefix,datestr(t,'yyyymmdd'),nc_suffix];
    if and(exist([ERA5_dir,'LSM_',datestr(t,'yyyymmdd'),'.nc'],'file')==0,exist(blkname,'file')==2)
        disp([blkname,' already exists!!'])
    else
        % Get the ERA5 horizontal grids (it should be the same for every file)
        % Use ERA5_offset days ago data from the latest ERA5 file
        nc=netcdf([ERA5_dir,'LSM_',datestr(now-ERA5_offset,'yyyymmdd'),'.nc']);
        disp(['Use this land file :',char([ERA5_dir,'LSM_',datestr(now-ERA5_offset,'yyyymmdd'),'.nc'])])

        lon1=nc{'lon'}(:);
        lat1=nc{'lat'}(:);
        [lon1,lat1]=meshgrid(lon1,lat1);

        mask=squeeze(nc{'LSM'}(1,:,:));
        mask(mask ~=0 ) = 1; %we take the first record
        mask = 1-mask ;
        mask(mask ==0 ) = NaN ;
        close(nc);

        % for M=mo_min:mo_max
        disp(' ')
        disp(['Processing date: ',datestr(t,'yyyymmdd')])
        disp(' ')
        %-------------------------------------------------------------------%
        %
        % Process time (here in days), with LSM file (common for frc/blk)
        %
        %-------------------------------------------------------------------%
        nc=netcdf([ERA5_dir,'LSM_',datestr(t,'yyyymmdd'),'.nc']);
        ERA5_time=nc{'time'}(:);
        close(nc);
        dt=mean(gradient(ERA5_time));
        disp(['dt=',num2str(dt)])
        %-----------------------------------------------------------
        %Variable overlapping timesteps : 2 at the beginning and 2 at the end
        %------------------------------------------------------------
        tlen0=length(ERA5_time);
        disp(['tlen0=',num2str(tlen0)])
        freq=1; % hourly
        itolap=freq*itolap_era5;
        tlen=tlen0+2*itolap;
        disp(['tlen=',num2str(tlen)])
        disp(['Overlap is ',num2str(itolap_era5),' records before and after'])     
        time=0*(1:tlen);
        time(itolap+1:tlen0+itolap)=ERA5_time;   
        disp(['====================='])
        disp('Compute time for croco file')
        disp(['====================='])
        for aa=1:itolap
            time(aa)=time(itolap+1)-(itolap+1-aa)*dt;
        end
        for aa=1:itolap
            time(tlen0+itolap+aa)=time(tlen0+itolap)+aa*dt;
        end
        
        %-------------------------------------------------------------------%
        %
        % Create the CROCO bulk forcing files
        %
        % ------------------------------------------------------------------%
        %
        disp(['====================='])
        disp('Create the blk/frc netcdf file')
        disp(['====================='])
        %
        if makeblk==1 
            disp(['Create a new bulk file: ' blkname])
            create_bulk(blkname,grdname,CROCO_title,time,0);
            disp([' '])
        end
        if makefrc==1
            disp(['Create a new forcing file: ' frcname])
            create_forcing(frcname,grdname,CROCO_title,...
                            time,time,time,...
                            time,time,time,...
                            0,0,0,0,0,0);
            disp([' '])
        end
        
        %
        % Open the CROCO forcing files
        if makefrc==1
            nc_frc=netcdf(frcname,'write');
        else
            nc_frc=[];
        end
        if makeblk==1
            nc_blk=netcdf(blkname,'write');
        else
            nc_blk=[];
        end 
        disp(' ')
        disp('======================================================================')
        disp(['Perform interpolations for ',datestr(t,'yyyy-mm-dd'),'               '])
        disp('======================================================================')
        disp(' ')
        
        % Perform interpolations for the current month
        %
        for tndx=1:tlen0
            if mod(tndx,6)==0
                disp(['Step: ',num2str(tndx),' of ',num2str(tlen0)])
            end
            interp_hindcast_ERA5(ERA5_dir,datestr(t,'yyyymmdd'),Roa,interp_method,lon1,lat1,...
                                mask,tndx,nc_frc,nc_blk,lon,lat,angle,tndx+itolap)	
        end
        
        %
        % Add the tides
        %
        if add_tides_fcst==1
            disp(' ')
            disp(['Create a new only tide forcing file: ' frcname])
            create_forcing_tideonly(frcname,grdname,CROCO_title)
            disp(['Add tidal data ... '])
            [Y,M,d,h,mi,s] = datevec(datestr(t,'yyyy-mm-dd'));
            add_tidal_data(tidename,grdname,frcname,Ntides,tidalrank, ...
                        Yorig,Y,M,coastfileplot,sal_tides,salname)
        end
        %
        % Close the CROCO forcing files
        %
        if ~isempty(nc_frc)
            close(nc_frc);
            end
        if ~isempty(nc_blk)
            close(nc_blk);
        end
    end
end

