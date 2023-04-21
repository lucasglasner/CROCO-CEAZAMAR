# -*- coding: UTF-8 -*-
#
# generated by wxGlade 0.8.0b3 on Tue Jan 30 13:49:27 2018
#

import numpy as np


###############################################################
# Ertel Potential vorticity


def get_pv(croco, tindex, depth=None, minlev=None, maxlev=None,
           lonindex=None, latindex=None, typ='ijk'):

    # mask = croco.wrapper.masks['mask_r']
    # pv = np.full_like(mask,np.nan)

    # pv from minlev to maxlev
    if depth is None:
        if lonindex is not None:
            var = np.full((maxlev - minlev, croco.wrapper.M), np.nan)
        else:
            var = np.full((maxlev - minlev, croco.wrapper.L), np.nan)
        var[:, 1:-1] = calc_ertel(croco, tindex, minlev=minlev, maxlev=maxlev,
                                  lonindex=lonindex, latindex=latindex, typ=typ)

    # pv on a depth
    elif depth <= 0:
        var = np.full((maxlev - minlev, croco.wrapper.M, croco.wrapper.L), np.nan)
        # pv = np.tile(pv,(maxlev-minlev,1,1))
        var[:, 1:-1, 1:-1] = calc_ertel(croco, tindex, minlev=minlev, maxlev=maxlev,
                                        lonindex=lonindex, latindex=latindex, typ=typ)

    # pv on a level
    elif depth > 0:
        var = np.full((croco.wrapper.M, croco.wrapper.L), np.nan)
        var[1:-1, 1:-1] = calc_ertel(croco, tindex,
                                     minlev=int(depth) - 2, maxlev=int(depth - 1),
                                     lonindex=lonindex, latindex=latindex, typ=typ)
    return var


def calc_ertel(croco, tindex, minlev=None, maxlev=None, lonindex=None, latindex=None, typ='ijk'):
    """
    #
    #   epv    - The ertel potential vorticity with respect to property 'lambda'
    #
    #                                       [ curl(u) + f ]
    #   -  epv is given by:           EPV = --------------- . del(lambda)
    #                                            rho
    #
    #   -  pvi,pvj,pvk - the x, y, and z components of the potential vorticity.
    #
    #   -  Ertel PV is calculated on horizontal rho-points, vertical w-points.
    #
    #
    #   tindex   - The time index at which to calculate the potential vorticity.
    #   depth    - depth
    #
    # Adapted from rob hetland.
    #
    """

    # Grid parameters

    minlon = 0 if lonindex is None else lonindex - 1
    maxlon = croco.wrapper.L - 1 if lonindex is None else lonindex + 1
    minlat = 0 if latindex is None else latindex - 1
    maxlat = croco.wrapper.M - 1 if latindex is None else latindex + 1

    pm = croco.wrapper.metrics['dx_r'][minlat:maxlat + 1, minlon:maxlon + 1]
    pm = np.tile(pm, (maxlev - minlev + 1, 1, 1))
    pn = croco.wrapper.metrics['dy_r'][minlat:maxlat + 1, minlon:maxlon + 1]
    pn = np.tile(pn, (maxlev - minlev + 1, 1, 1))
    f = croco.wrapper.metrics['f'][minlat:maxlat + 1, minlon:maxlon + 1]
    f = np.tile(f, (maxlev - minlev + 1, 1, 1))
    rho0 = croco.rho0
    #
    # 3D variables
    #

    ssh = croco.variables['ssh'].isel(t=tindex, y_r=slice(minlat, maxlat + 1),
                                      x_r=slice(minlon, maxlon + 1)).values
    dz = croco.wrapper.scoord2dz_r(ssh, alpha=0., beta=0, lonindex=lonindex,
                                   latindex=latindex)[minlev:maxlev + 1, :]
    u = croco.variables['u'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                  y_r=slice(minlat, maxlat + 1), x_u=slice(minlon, maxlon)).values
    v = croco.variables['v'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                  y_v=slice(minlat, maxlat), x_r=slice(minlon, maxlon + 1)).values
    w = croco.variables['w'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                  y_r=slice(minlat, maxlat + 1), x_r=slice(minlon, maxlon + 1))

    try:
        rho = croco.variables['rho'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                          y_r=slice(minlat, maxlat + 1), x_r=slice(minlon, maxlon + 1))
    except Exception:
        # temp = croco.variables['temp'].isel(t=tindex, z_r=slice(minlev,maxlev+1),\
                    # y_r=slice(minlat,maxlat+1), x_r=slice(minlon,m axlon + 1))
        # salt = croco.variables['salt'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),\
                    # y_r=slice(minlat,maxlat+1),x_r=slice(minlon,maxlon+1))
        # rho = croco.rho_eos(temp,salt,0)
        print('rho not in history file')
        return

    if 'k' in typ:
        #
        #
        #  Ertel potential vorticity, term 1: [f + (dv/dx - du/dy)]*drho/dz
        #
        # Compute d(v)/d(xi) at PSI-points.
        #
        dxm1 = 0.25 * (pm[:, :-1, 1:] + pm[:, 1:, 1:] + pm[:, :-1, :-1] + pm[:, 1:, :-1])
        dvdxi = np.diff(v, n=1, axis=2) * dxm1
        #
        #  Compute d(u)/d(eta) at PSI-points.
        #
        dym1 = 0.25 * (pn[:, :-1, 1:] + pn[:, 1:, 1:] + pn[:, :-1, :-1] + pn[:, 1:, :-1])
        dudeta = np.diff(u, n=1, axis=1) * dym1
        #
        #  Compute d(rho)/d(z) at horizontal RHO-points and vertical W-points
        #
        dz_w = 0.5 * (dz[:-1, :, :] + dz[1:, :, :])
        drhodz = np.diff(rho, n=1, axis=0) / dz_w
        #
        #  Compute Ertel potential vorticity <k hat> at horizontal RHO-points and
        #  vertical W-points.
        omega = dvdxi - dudeta
        omega = f[:, 1:-1, 1:-1] + 0.25 * (omega[:, :-1, 1:] + omega[:, 1:, 1:] +
                                           omega[:, :-1, :-1] + omega[:, 1:, :-1])
        pvk = 0.5 * (omega[:-1, :, :] + omega[1:, :, :]) * drhodz[:, 1:-1, 1:-1]
    else:
        pvk = 0.

    if 'i' in typ:
        #
        #
        #  Ertel potential vorticity, term 2: (dw/dy - dv/dz)*(drho/dx)
        #
        #  Compute d(w)/d(y) at horizontal V-points and vertical RHO-points
        #
        dym1 = 0.5 * (pn[:, :-1, :] + pn[:, 1:, :])
        dwdy = np.diff(w, axis=1) * dym1
        #
        #  Compute d(v)/d(z) at horizontal V-points and vertical W-points
        #
        dz_v = 0.5 * (dz[:, 1:, :] + dz[:, :-1, :])
        # dz_v = croco.crocoGrid.scoord2dz_v(zeta, alpha=0., beta=0.)
        dvdz = np.diff(v, axis=0) / (0.5 * (dz_v[:-1, :, :] + dz_v[1:, :, :]))
        #
        #  Compute d(rho)/d(xi) at horizontal U-points and vertical RHO-points
        #
        dxm1 = 0.5 * (pm[:, :, 1:] + pm[:, :, :-1])
        drhodx = np.diff(rho, axis=2) * dxm1
        #
        #  Add in term 2 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points.
        #
        pvi = (0.25 * (dwdy[1:, :-1, 1:-1] + dwdy[1:, 1:, 1:-1] + dwdy[:-1, :-1, 1:-1] + dwdy[:-1, 1:, 1:-1]) -
               0.5 * (dvdz[:, :-1, 1:-1] + dvdz[:, 1:, 1:-1])) *\
            0.25 * (drhodx[1:, 1:-1, :-1] + drhodx[1:, 1:-1, 1:] +
                    drhodx[:-1, 1:-1, :-1] + drhodx[:-1, 1:-1, 1:])
    else:
        pvi = 0.

    if 'j' in typ:
        #
        #
        #  Ertel potential vorticity, term 3: (du/dz - dw/dx)*(drho/dy)
        #
        #  Compute d(u)/d(z) at horizontal U-points and vertical W-points
        #
        dz_u = 0.5 * (dz[:, :, 1:] + dz[:, :, :-1])
        dudz = np.diff(u, axis=0) / (0.5 * (dz_u[:-1, :, :] + dz_u[1:, :, :]))
        #
        #  Compute d(w)/d(x) at horizontal U-points and vertical RHO-points
        #
        dxm1 = 0.5 * (pm[:, :, 1:] + pm[:, :, :-1])
        dwdx = np.diff(w, axis=2) * dxm1
        #
        #  Compute d(rho)/d(eta) at horizontal V-points and vertical RHO-points
        #
        dym1 = 0.5 * (pn[:, 1:, :] + pn[:, :-1, :])
        drhodeta = np.diff(rho, axis=1) * dym1
        #
        #  Add in term 3 contribution to Ertel potential vorticity at horizontal RHO-points and
        #  vertical W-points..
        #
        pvj = (0.5 * (dudz[:, 1:-1, 1:] + dudz[:, 1:-1, :-1]) -
               0.25 * (dwdx[1:, 1:-1, 1:] + dwdx[1:, 1:-1, :-1] + dwdx[:-1, 1:-1, 1:] + dwdx[:-1, 1:-1, :-1])) * \
            0.25 * (drhodeta[1:, :-1, 1:-1] + drhodeta[1:, 1:, 1:-1] +
                    drhodeta[:-1, :-1, 1:-1] + drhodeta[:-1, 1:, 1:-1])
    else:
        pvj = 0.

    #
    #
    # Sum potential vorticity components, and divide by rho0
    #
    pvi = pvi / rho0
    pvj = pvj / rho0
    pvk = pvk / rho0
    #
    return(np.squeeze(pvi + pvj + pvk))
    #
    #
    ####################################################################################


###############################################################
# Zeta_k term

def get_zetak(croco, tindex, depth=None, minlev=None, maxlev=None,
              lonindex=None, latindex=None):

    # mask = croco.wrapper.masks['mask_r']
    # pv = np.full_like(mask,np.nan)

    # pv from level 0 to N at latitude or longitude index
    if depth is None:
        if lonindex is not None:
            var = np.full((maxlev - minlev, croco.wrapper.M), np.nan)
        else:
            var = np.full((maxlev - minlev, croco.wrapper.L), np.nan)
        var[:, 1:-1] = calc_zetak(croco, tindex, minlev=minlev, maxlev=maxlev - 1,
                                  lonindex=lonindex, latindex=latindex)
    elif depth <= 0:
        var = np.full((maxlev - minlev, croco.wrapper.M, croco.wrapper.L), np.nan)
        var[:, 1:-1, 1:-1] = calc_zetak(croco, tindex, minlev=minlev, maxlev=maxlev - 1,
                                        lonindex=lonindex, latindex=latindex)
    # pv on a level
    elif depth > 0:
        var = np.full((croco.wrapper.M, croco.wrapper.L), np.nan)
        var[1:-1, 1:-1] = calc_zetak(croco, tindex,
                                     minlev=int(depth) - 1, maxlev=int(depth - 1),
                                     lonindex=lonindex, latindex=latindex)
    return var


def calc_zetak(croco, tindex, minlev=None, maxlev=None, lonindex=None, latindex=None):
    """
    #   -  zetak is given by:      (dv/dx - du/dy)/f
    #
    #   -  zetak is calculated at RHO-points
    #
    #
    #   tindex   - The time index at which to calculate the potential vorticity.
    #   depth    - depth
    #
    # Adapted from rob hetland.
    #
    """
    #
    # Grid parameters
    #
    minlon = 0 if lonindex is None else lonindex - 1
    maxlon = croco.wrapper.L - 1 if lonindex is None else lonindex + 1
    minlat = 0 if latindex is None else latindex - 1
    maxlat = croco.wrapper.M - 1 if latindex is None else latindex + 1

    pm = croco.wrapper.metrics['dx_r'][minlat:maxlat + 1, minlon:maxlon + 1]
    pm = np.tile(pm, (maxlev - minlev + 1, 1, 1))
    pn = croco.wrapper.metrics['dy_r'][minlat:maxlat + 1, minlon:maxlon + 1]
    pn = np.tile(pn, (maxlev - minlev + 1, 1, 1))
    f = croco.wrapper.metrics['f'][minlat:maxlat + 1, minlon:maxlon + 1]
    f = np.tile(f, (maxlev - minlev + 1, 1, 1))
    #
    # 3D variables
    #

    u = croco.variables['u'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                  y_r=slice(minlat, maxlat + 1), x_u=slice(minlon, maxlon)).values
    v = croco.variables['v'].isel(t=tindex, z_r=slice(minlev, maxlev + 1),
                                  y_v=slice(minlat, maxlat), x_r=slice(minlon, maxlon + 1)).values

    #
    #
    #  Ertel potential vorticity, term 1: (dv/dx - du/dy)/f
    #
    # Compute d(v)/d(xi) at PSI-points.
    #
    dxm1 = 0.25 * (pm[:, :-1, 1:] + pm[:, 1:, 1:] + pm[:, :-1, :-1] + pm[:, 1:, :-1])
    dvdxi = np.diff(v, n=1, axis=2) * dxm1
    #
    #  Compute d(u)/d(eta) at PSI-points.
    #
    dym1 = 0.25 * (pn[:, :-1, 1:] + pn[:, 1:, 1:] + pn[:, :-1, :-1] + pn[:, 1:, :-1])
    dudeta = np.diff(u, n=1, axis=1) * dym1
    #
    #  Compute Ertel potential vorticity <k hat> at horizontal RHO-points and
    #  vertical RHO-points.
    omega = dvdxi - dudeta
    return(np.squeeze(0.25 * (omega[:, :-1, 1:] + omega[:, 1:, 1:] +
                      omega[:, :-1, :-1] + omega[:, 1:, :-1]) / f[:, 1:-1, 1:-1]))

###############################################################
# dtdz term

def get_dtdz(croco, tindex, depth=None, minlev=None, maxlev=None, lonindex=None, latindex=None):

    # dtdz from levels 1 to N at a given longitude or latitude
    if depth is None:
        if lonindex is not None:
            dtdz = np.full((maxlev - minlev, croco.wrapper.M), np.nan)
        else:
            dtdz = np.full((maxlev - minlev, croco.wrapper.L), np.nan)
        dtdz[:] = calc_dtdz(croco, tindex, minlev=minlev, maxlev=maxlev,
                            lonindex=lonindex, latindex=latindex)
    # dtdz from levels minlev to maxlev
    elif depth <= 0:
        dtdz = np.full((maxlev - minlev, croco.wrapper.M, croco.wrapper.L), np.nan)
        dtdz[:] = calc_dtdz(croco, tindex, minlev=minlev, maxlev=maxlev)
    # dtdz on a level
    elif depth > 0:
        dtdz = np.full((croco.wrapper.M, croco.wrapper.L), np.nan)
        dtdz[:] = calc_dtdz(croco, tindex, minlev=int(depth) - 2, maxlev=int(depth) - 1)
    return dtdz


def calc_dtdz(croco, tindex, minlev=None, maxlev=None, lonindex=None, latindex=None):

    #
    # 3D variables
    #
    if lonindex is not None:
        ssh = croco.variables['ssh'].isel(t=tindex,
                                          x_r=slice(lonindex - 1, lonindex + 2)).values
        dz = np.squeeze(croco.wrapper.scoord2dz_r(ssh, alpha=0., beta=0, lonindex=lonindex)
                        [minlev:maxlev + 1, :, 1:-1])
        t = croco.variables['temp'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), x_r=lonindex)
    elif latindex is not None:
        ssh = croco.variables['ssh'].isel(t=tindex,
                                          y_r=slice(latindex - 1, latindex + 2)).values
        dz = np.squeeze(croco.wrapper.scoord2dz_r(ssh, alpha=0., beta=0, latindex=latindex)
                        [minlev:maxlev + 1, 1:-1, :])
        t = croco.variables['temp'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), y_r=latindex)
    else:
        ssh = croco.variables['ssh'].isel(t=tindex).values
        dz = croco.wrapper.scoord2dz_r(ssh, alpha=0., beta=0)[minlev:maxlev + 1, :, :]
        t = croco.variables['temp'].isel(t=tindex, z_r=slice(minlev, maxlev + 1))
    dtdz = np.diff(t, axis=0) / (0.5 * (dz[:-1, :] + dz[1:, :]))
    return(dtdz)


###############################################################
# Richardson Number

def get_richardson(croco, tindex, depth=None, minlev=None, maxlev=None,
                   lonindex=None, latindex=None):

    # mask = croco.wrapper.masks['mask_r']
    # pv = np.full_like(mask,np.nan)

    # Ri from level 0 to N at latitude or longitude index
    if depth is None:
        if lonindex is not None:
            var = np.full((maxlev - minlev, croco.wrapper.M), np.nan)
        else:
            var = np.full((maxlev - minlev, croco.wrapper.L), np.nan)
        var[:, 1:-1] = calc_richardson(croco, tindex, minlev=minlev, maxlev=maxlev,
                                       lonindex=lonindex, latindex=latindex)
    # pv at depth
    elif depth <= 0:
        var = np.full((maxlev - minlev, croco.wrapper.M, croco.wrapper.L), np.nan)
        var[:, 1:-1, 1:-1] = calc_richardson(croco, tindex, minlev=minlev, maxlev=maxlev)
    # pv on a level
    elif depth > 0:
        var = np.full((croco.wrapper.M, croco.wrapper.L), np.nan)
        var[1:-1, 1:-1] = calc_richardson(croco, tindex,
                                          minlev=int(depth) - 2, maxlev=int(depth) - 1)
    return var


def calc_richardson(croco, tindex, minlev=None, maxlev=None, lonindex=None, latindex=None):
    """
      -  Ri is given by:      N²/((du/dz)² - (dv/dz)²)
         with N = sqrt(-g/rho0 * drho/dz)
    
      -  Ri is calculated at RHO-points and w level
    
    
      tindex   - The time index at which to calculate the potential vorticity.
      depth    - depth
    
    """

    # If u or v or rho not in netcdf file, abort
    try:
        croco.variables['rho']
        croco.variables['u']
        croco.variables['v']
    except Exception:
        print("Variable rho, u or v missing in the netcdf file")
        return

    # Longitude section
    if lonindex is not None:
        ssh = croco.variables['ssh'].isel(t=tindex,
                                          x_r=slice(lonindex - 1, lonindex + 2)).values
        z_r = np.squeeze(croco. get_coord("rho", direction='z', timeIndex=tindex)
                         [:, :, lonindex])
        z_u = np.squeeze(croco. get_coord("u", direction='z', timeIndex=tindex)
                         [:, :, lonindex])
        z_v = np.squeeze(croco. get_coord("v", direction='z', timeIndex=tindex)
                         [:, :, lonindex])
        rho = croco.variables['rho'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), x_r=lonindex)
        u = croco.variables['u'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), x_u=lonindex)
        v = croco.variables['v'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), x_r=lonindex)
        drhodz = np.diff(rho, axis=0) / np.diff(z_r, axis=0)
        N2 = (-croco.g / croco.rho0) * (drhodz)
        dudz = np.diff(u, axis=0) / np.diff(z_u, axis=0)
        dvdz = np.diff(v, axis=0) / np.diff(z_v, axis=0)
        Ri = np.squeeze(np.log10(N2[:, 1:-1] / (
                        (0.5 * (dudz[:, 1:-1] + dudz[:, 1:-1]))**2 +
                        (0.5 * (dvdz[:, :-1] + dvdz[:, 1:]))**2)))

    # Latitude section
    elif latindex is not None:
        ssh = croco.variables['ssh'].isel(t=tindex,
                                          y_r=slice(latindex - 1, latindex + 2)).values
        z_r = np.squeeze(croco. get_coord("rho", direction='z', timeIndex=tindex)
                         [:, latindex, :])
        z_u = np.squeeze(croco. get_coord("u", direction='z', timeIndex=tindex)
                         [:, latindex, :])
        z_v = np.squeeze(croco. get_coord("v", direction='z', timeIndex=tindex)
                         [:, latindex, :])
        rho = croco.variables['rho'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), y_r=latindex)
        u = croco.variables['u'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), y_r=latindex)
        v = croco.variables['v'].isel(t=tindex, z_r=slice(minlev, maxlev + 1), y_v=latindex)
        drhodz = np.diff(rho, axis=0) / np.diff(z_r, axis=0)
        N2 = (-croco.g / croco.rho0) * (drhodz)
        dudz = np.diff(u, axis=0) / np.diff(z_u, axis=0)
        dvdz = np.diff(v, axis=0) / np.diff(z_v, axis=0)
        Ri = np.squeeze(np.log10(N2[:, 1:-1] / (
                        (0.5 * (dudz[:, 1:] + dudz[:, :-1]))**2 +
                        (0.5 * (dvdz[:, 1:-1] + dvdz[:, 1:-1]))**2)))

    # Level or depth section 
    else:
        ssh = croco.variables['ssh'].isel(t=tindex).values
        dz = croco.wrapper.scoord2dz_r(ssh, alpha=0., beta=0)[minlev:maxlev + 1, :, :]
        z_r = np.squeeze(croco. get_coord("rho", direction='z', timeIndex=tindex)
                         [minlev:maxlev + 1, :, :])
        z_u = np.squeeze(croco. get_coord("u", direction='z', timeIndex=tindex)
                         [minlev:maxlev + 1, :, :])
        z_v = np.squeeze(croco. get_coord("v", direction='z', timeIndex=tindex)
                         [minlev:maxlev + 1, :, :])
        rho = croco.variables['rho'].isel(t=tindex, z_r=slice(minlev, maxlev + 1))
        u = croco.variables['u'].isel(t=tindex, z_r=slice(minlev, maxlev + 1))
        v = croco.variables['v'].isel(t=tindex, z_r=slice(minlev, maxlev + 1))
        drhodz = np.diff(rho, axis=0) / np.diff(z_r, axis=0)
        N2 = (-croco.g / croco.rho0) * (drhodz)
        dudz = np.diff(u, axis=0) / np.diff(z_u, axis=0)
        dvdz = np.diff(v, axis=0) / np.diff(z_v, axis=0)
        Ri = np.squeeze(np.log10(N2[:, 1:-1, 1:-1] / (
                        (0.5 * (dudz[:, 1:-1, 1:] + dudz[:, 1:-1, :-1]))**2 +
                        (0.5 * (dvdz[:, :-1, 1:-1] + dvdz[:, 1:, 1:-1]))**2)))

    return(Ri)