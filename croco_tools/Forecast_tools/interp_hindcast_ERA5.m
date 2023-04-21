function interp_ERA5(ERA5_dir,timestamp,Roa,interp_method,...
                     lon1,lat1,mask1,tin,...
		     nc_frc,nc_blk,lon,lat,angle,tout)
%
% Read the local daily ERA5 files and perform the space interpolations
%
%  Illig, 2010, from interp_NCEP
%  Updated    January 2016 (E. Cadier and S. Illig)
%  Updated    D.Donoso, G. Cambon. P. Penven (Oct 2021) 
%  Updated    L.Glasner (Apr 2023)
%---------------------------------------------------------------------------------

% % Fill the CROCO files
if ~isempty(nc_frc)
  %
  % 1: Sea surface temperature: Convert from Kelvin to Celsius
  %
  % vname='SST';
  % varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  % nc=netcdf(varpath);
  % sst=squeeze(nc{vname}(tin,:,:));
  % close(nc);
  % %sst=get_missing_val(lon1,lat1,mask1.*sst,nan,Roa,nan);
  % sst=sst-273.15;
  % sst=interp2(lon1,lat1,sst,lon,lat,interp_method);
  % nc_frc{'SST'}(tout,:,:)=sst;


  %
  % 2: Water fluxes: Precipitation - Evaporation 
  % Convert from [kg/m^2/s] to cm/day
  %
  vname='TP';
  tp_varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(tp_varpath);
  prate=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %prate=get_missing_val(lon1,lat1,mask1.*prate,nan,Roa,nan);
  prate=prate*0.1*(24.*60.*60.0);
  prate=interp2(lon1,lat1,prate,lon,lat,interp_method);
  prate(prate<1.e-4)=0;

  vname='E';
  e_varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(e_varpath);
  erate=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %erate=get_missing_val(lon1,lat1,mask1.*erate,nan,Roa,nan);
  erate=erate*0.1*(24.*60.*60.0);
  erate=interp2(lon1,lat1,erate,lon,lat,interp_method);
  erate(erate<1.e-4)=0;
  
  swflux=erate-prate;
  nc_frc{'swflux'}(tout,:,:)=swflux; % <- WRITE

  %
  % 3: Momentum fluxes (wind surface stress)
  %
  cosa=cos(angle);
  sina=sin(angle);
  vname='EWSS';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  tx=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %tx=get_missing_val(lon1,lat1,mask1.*tx,nan,0,nan);
  tx=interp2(lon1,lat1,tx,lon,lat,interp_method);
  
  vname='NSSS';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  ty=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %ty=get_missing_val(lon1,lat1,mask1.*ty,nan,0,nan);
  ty=interp2(lon1,lat1,ty,lon,lat,interp_method);

  sustr=rho2u_2d(tx.*cosa+ty.*sina);
  svstr=rho2v_2d(ty.*cosa-tx.*sina);
  nc_frc{'sustr'}(tout,:,:)=sustr;% => ERA5 downloaded
  nc_frc{'svstr'}(tout,:,:)=svstr;% => ERA5 downloaded


  %
  % 4: Heat fluxes
  %
  % 4.1: Shortwave flux: [W/m^2]
  %      CROCO convention: downward = positive
  %
  %       Net solar shortwave radiation
  %
  vname='SSR';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  radsw=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %radsw=get_missing_val(lon1,lat1,mask1.*dswrf,nan,Roa,nan);
  radsw=interp2(lon1,lat1,radsw,lon,lat,interp_method);
  radsw(radsw<1.e-10)=0;
  nc_frc{'swrad'}(tout,:,:)=radsw;


  %
  %  4.2 Get the net longwave flux [W/m^2]
  %
  vname='STR';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  radlw=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %radlw=get_missing_val(lon1,lat1,mask1.*dlwrf,nan,Roa,nan);
  radlw=interp2(lon1,lat1,radlw,lon,lat,interp_method);

  %
  %  4.3 Get the sensible heat flux [W/m^2]
  %
  vname='SSHF';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  sshf=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %sshf=get_missing_val(lon1,lat1,mask1.*sshf,nan,Roa,nan);
  sshf=interp2(lon1,lat1,sshf,lon,lat,interp_method);


  %
  %  4.4 Get the latent heat flux [W/m^2]
  %
  vname='SLHF';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  slhf=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %slhf=get_missing_val(lon1,lat1,mask1.*slhf,nan,Roa,nan);
  slhf=interp2(lon1,lat1,slhf,lon,lat,interp_method);

  %
  % 4.5 STORE NET HEAT FLUX
  %
  netheatflux=radsw+radlw+sshf+slhf;
  nc_frc{'shflux'}(tout,:,:)=netheatflux;

  %
  % 5: COMPUTE HEAT SENSITIVITY TO SST (dQdSST)
  %

  %
  % 5.1: Get air temperature (°C)
  %
  vname='T2M';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  tair=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %tair=get_missing_val(lon1,lat1,mask1.*tair,nan,Roa,nan);
  tair=tair-273.15;
  tair=interp2(lon1,lat1,tair,lon,lat,interp_method);

  %
  % 5.2: Relative humidity: Get Specific Humidity [Kg/Kg]
  %
  vname='Q';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  shum=squeeze(nc{vname}(tin,:,:));
  close(nc);
  %shum=get_missing_val(lon1,lat1,mask1.*shum,nan,Roa,nan);
  shum=interp2(lon1,lat1,shum,lon,lat,interp_method);
  % computes specific humidity at saturation (Tetens  formula)
  % (see air_sea tools, fonction qsat)
  rhum=shum./qsat(tair);

  %
  % 5.3: Air density at sea level: [kg m-3]
  %
  rhoa=air_dens(tair,rhum*100);

  %
  % 5.4: Wind speed (m s-1)
  %
  % 5.4.1: 
  vname='U10M';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  uwnd=squeeze(nc{vname}(tin,:,:));
  close(nc)
  %uwnd=get_missing_val(lon1,lat1,mask1.*uwnd,nan,Roa,nan);
  uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
  % 5.4.2:
  vname='V10M';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  vwnd=squeeze(nc{vname}(tin,:,:));
  close(nc)
  %vwnd=get_missing_val(lon1,lat1,mask1.*vwnd,nan,Roa,nan);
  vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
  % 5.4.3:
  wspd=sqrt(uwnd.^2+vwnd.^2);

  %
  % 5.5: Store dQdSST sensitivity
  %
  dqdsst=get_dqdsst(sst,tair,rhoa,wspd,shum);
  nc_frc{'dQdSST'}(tout,:,:)=dqdsst;


%   nc_frc{'Awave'}(tout,:,:)=awave; => ERA5 not downloaded
%   nc_frc{'Dwave'}(tout,:,:)=dwave; => ERA5 not downloaded
%   nc_frc{'Pwave'}(tout,:,:)=pwave; => ERA5 not downloaded
end


if ~isempty(nc_blk)
  %
  % 1: Air temperature: Convert from Kelvin to Celsius
  %
  vname='T2M';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  tair=squeeze(nc{vname}(tin,:,:));
  close(nc);
  tair=tair-273.15;
  tair=interp2(lon1,lat1,tair,lon,lat,interp_method);
  nc_blk{'tair'}(tout,:,:)=tair; %<- WRITE

  %
  % 2: Relative humidity: Get Dew Point [K]
  %
  vname='D2M';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  dewp=squeeze(nc{vname}(tin,:,:));
  dewp=dewp-273.15;
  close(nc);
  dewp=interp2(lon1,lat1,dewp,lon,lat,interp_method);
  % computes vapor pressure at saturation (Clausius-Clapeyron formula)
  % (see air_sea tools, fonction satvap)
  rhum=satvap(dewp)./satvap(tair);
  nc_blk{'rhum'}(tout,:,:)=rhum; % <- WRITE


  %
  % 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
  %
  vname='TP';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  prate=squeeze(nc{vname}(tin,:,:));
  close(nc);
  prate=prate*0.1*(24.*60.*60.0);
  prate=interp2(lon1,lat1,prate,lon,lat,interp_method);
  prate(prate<1.e-4)=0;
  nc_blk{'prate'}(tout,:,:)=prate; % <- WRITE

  %
  % 4: Shortwave flux: [W/m^2]
  %      CROCO convention: downward = positive
  %
  %  Net solar shortwave radiation
  %
  vname='SSR';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  radsw=squeeze(nc{vname}(tin,:,:));
  close(nc);
  radsw=interp2(lon1,lat1,radsw,lon,lat,interp_method);
  radsw(radsw<1.e-10)=0;
  nc_blk{'radsw'}(tout,:,:)=radsw;
  % 5: Longwave flux:  [W/m^2]
  %      CROCO convention: positive upward.

  %
  %  5.1 Get the net longwave flux [W/m^2]
  %
  vname='STR';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf(varpath);
  radlw=squeeze(nc{vname}(tin,:,:));
  close(nc);
  radlw=interp2(lon1,lat1,radlw,lon,lat,interp_method);
  nc_blk{'radlw'}(tout,:,:)=-1*radlw; % <- WRITE

  %
  %  5.2 get the downward longwave flux [W/m^2]
  %
  vname='STRD';
  varpath=[ERA5_dir,vname,'_',timestamp,'.nc'];
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  radlw_in=squeeze(nc{vname}(tin,:,:));
  close(nc);
  radlw_in=interp2(lon1,lat1,radlw_in,lon,lat,interp_method);
  nc_blk{'radlw_in'}(tout,:,:)=radlw_in; % <- WRITE
  
  %
  % 6: Wind  [m/s]
  %
  % Rotation matrix for the croco grid
  cosa=cos(angle);
  sina=sin(angle);
  upath= [ERA5_dir,'U10M','_',timestamp,'.nc'];
  vpath= [ERA5_dir,'V10M','_',timestamp,'.nc'];

  vname='U10M';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  uwnd=squeeze(nc{vname}(tin,:,:));
  close(nc)
  uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
  %
  %
  vname='V10M';
  nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  vwnd=squeeze(nc{vname}(tin,:,:));
  close(nc)
  vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);

  u10=rho2u_2d(uwnd.*cosa+vwnd.*sina);
  v10=rho2v_2d(vwnd.*cosa-uwnd.*sina);
  wspd=sqrt(uwnd.^2+vwnd.^2);

  nc_blk{'uwnd'}(tout,:,:)=u10;
  nc_blk{'vwnd'}(tout,:,:)=v10;
  nc_blk{'wspd'}(tout,:,:)=wspd;
  %   %
  %   % Compute the stress
  %   %
  % ewss_path=[ERA5_dir,'EWSS','_',timestamp,'.nc'];
  % nwss_path=[ERA5_dir,'NSSS','_',timestamp,'.nc'];

  % [Cd,uu]=cdnlp(wspd,10.);
  % rhoa=air_dens(tair,rhum*100);
  % tx=Cd.*rhoa.*uwnd.*wspd;
  % ty=Cd.*rhoa.*vwnd.*wspd;
  % sustr=rho2u_2d(tx.*cosa+ty.*sina);
  % svstr=rho2v_2d(ty.*cosa-tx.*sina);
  % nc_blk{'sustr'}(tout,:,:)=sustr;% => ERA5 not downloaded
  % nc_blk{'svstr'}(tout,:,:)=svstr;% => ERA5 not downloaded
  % vname='EWSS';
  % nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  % tx=squeeze(nc{vname}(tin,:,:));
  % close(nc)
  % tx=get_missing_val(lon1,lat1,mask1.*tx,nan,0,nan);
  % tx=interp2(lon1,lat1,tx,lon,lat,interp_method);
  % %
  % vname='NSSS';
  % nc=netcdf([ERA5_dir,vname,'_',timestamp,'.nc']);
  % ty=squeeze(nc{vname}(tin,:,:));
  % close(nc)
  % ty=get_missing_val(lon1,lat1,mask1.*ty,nan,0,nan);
  % ty=interp2(lon1,lat1,ty,lon,lat,interp_method);

  % sustr=rho2u_2d(tx.*cosa+ty.*sina);
  % svstr=rho2v_2d(ty.*cosa-tx.*sina);
  % nc_blk{'sustr'}(tout,:,:)=sustr;% => ERA5 downloaded
  % nc_blk{'svstr'}(tout,:,:)=svstr;% => ERA5 downloaded

end
end

