import numpy as np
from datetime import datetime, timedelta
from pyproj import Proj
import xarray
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os, sys, glob
import sys

os.chdir('/media/adrian/runs/g16pedro/CONUS')

for ff in sorted(glob.glob('OR_ABI-L2-MCMIPC-M6_G16_*.nc')):

    # Using xarray, I assign the opened file to the variable C for the CONUS domain.
    #
    FILE = './'+ff

    C = xarray.open_dataset(FILE)


    # Scan's start time, converted to datetime object
    scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

    # Scan's end time, converted to datetime object
    scan_end = datetime.strptime(C.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')

    # File creation time, convert to datetime object
    file_created = datetime.strptime(C.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')

    # The 't' variable is the scan's midpoint time
    # I'm not a fan of numpy datetime, so I convert it to a regular datetime object
    midpoint = str(C['t'].data)[:-8]
    scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

    print('Scan Start    : %s' % scan_start)
    print('Scan midpoint : %s' % scan_mid)
    print('Scan End      : %s' % scan_end)
    print('File Created  : %s' % file_created)
    print('Scan Duration : %.2f minutes' % ((scan_end-scan_start).seconds/60))


    # Confirm that each band is the wavelength we are interested in
    for band in [3, 5, 7, 13]:
        print("%s is %.2f %s (units = %s)" % (C['band_wavelength_C%02d' % band].long_name,
                                              C['band_wavelength_C%02d' % band][0],
                                              C['band_wavelength_C%02d' % band].units,
                                              C['CMI_C%02d' % band].units))
        
    # NOTE: Units K is in Kelvin and 1 is in Radiance between 0-1


    # Load the three channels into appropriate R, G, and B variables
    R = C['CMI_C15'].data - C['CMI_C13'].data  # 
    G = C['CMI_C13'].data - C['CMI_C07'].data
    B = C['CMI_C13'].data

    # Normalize each channel by the appropriate range of values  e.g. R = (R-minimum)/(maximum-minimum)
    # Normalizing the channels forces the values to range from 0-1, what we need to create the RGB.
#    R = (R-0)/(1-0)
#    G = (G-0)/(0.7-0)
#    B = (B-0)/(30-0)

#    # Apply range limits for each channel. RGB values must be between 0 and 1
#    R = np.clip(R, 0, 1)
#    G = np.clip(G, 0, 1)
#    B = np.clip(B, 0, 1)

#    # Apply a gamma correction to the image
    gamma = 1.0
    R = np.power(R, 1/gamma)
    G = np.power(G, 1/gamma)
    B = np.power(B, 1/gamma)


    #fig, ([ax1, ax2, ax3]) = plt.subplots(1, 3, figsize=(16, 3))

    #ax1.imshow(R, cmap='Reds', vmax=1, vmin=0)
    #ax1.set_title('Red', fontweight='semibold')
    #ax1.axis('off')

    #ax2.imshow(G, cmap='Greens', vmax=1, vmin=0)
    #ax2.set_title('Green', fontweight='semibold')
    #ax2.axis('off')

    #ax3.imshow(B, cmap='Blues', vmax=1, vmin=0)
    #ax3.set_title('Blue', fontweight='semibold')
    #ax3.axis('off')

    #plt.subplots_adjust(wspace=.02)


    # For fun, lets get the True Color RGB
    def get_TrueColor_RGB():
        R = C['CMI_C02'].data
        G = C['CMI_C03'].data
        B = C['CMI_C01'].data

        # Apply range limits for each channel. RGB values must be between 0 and 1
        R = np.clip(R, 0, 1)
        G = np.clip(G, 0, 1)
        B = np.clip(B, 0, 1)

        # Apply a gamma correction to the image
        gamma = 2.2
        R = np.power(R, 1/gamma)
        G = np.power(G, 1/gamma)
        B = np.power(B, 1/gamma)

        # Calculate the "True" Green
        G_true = 0.45 * R + 0.1 * G + 0.45 * B
        G_true = np.maximum(G_true, 0)
        G_true = np.minimum(G_true, 1)

        return np.dstack([R, G_true, B])

    RGB_true = get_TrueColor_RGB()

    # The RGB array for the true color image
    RGB = np.dstack([R, G, B])

    #fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16,6))

    #ax1.imshow(RGB)
    #ax1.set_title('GOES-16 RGB Day snow-fog', fontweight='semibold', loc='left', fontsize=12);
    #ax1.set_title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');
    #ax1.axis('off');

    #ax2.imshow(RGB_true)
    #ax2.set_title('GOES-16 RGB True Color', fontweight='semibold', loc='left', fontsize=12);
    #ax2.set_title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');
    #ax2.axis('off');




    #####################




    C['goes_imager_projection']

    ######################################################################
    #

    # Satellite height
    sat_h = C['goes_imager_projection'].perspective_point_height

    # Satellite longitude
    sat_lon = C['goes_imager_projection'].longitude_of_projection_origin

    # Satellite sweep
    sat_sweep = C['goes_imager_projection'].sweep_angle_axis

    semi_major = C['goes_imager_projection'].semi_major_axis
    semi_minor = C['goes_imager_projection'].semi_minor_axis

    # The projection x and y coordinates equals the scanning angle (in radians)
    # multiplied by the satellite height See details here:
    # https://proj4.org/operations/projections/geos.html?highlight=geostationary
    x = C['x'][:] * sat_h
    y = C['y'][:] * sat_h



    #fig = plt.figure(figsize=(15, 12))

    #globe = ccrs.Globe(semimajor_axis=semi_major, semiminor_axis=semi_minor)
    #geos = ccrs.Geostationary(central_longitude=sat_lon, 
    #                         satellite_height=sat_h, globe=globe)

    #ax = fig.add_subplot(1, 1, 1, projection=geos)

    #ax.imshow(RGB, origin='upper',
    #                   extent=(x.min(), x.max(), y.min(), y.max()),
    #                   transform=geos,
    #                   interpolation='nearest', vmin=162., vmax=330.)
    #ax.coastlines(resolution='50m', color='black', linewidth=2)
    #ax.add_feature(ccrs.cartopy.feature.STATES)

    #plt.title('GOES-16 Day Snow-Fog', loc='left', fontweight='semibold', fontsize=15)
    #plt.title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');

    globe = ccrs.Globe(semimajor_axis=semi_major, semiminor_axis=semi_minor)
    geos = ccrs.Geostationary(central_longitude=sat_lon, 
                             satellite_height=sat_h, globe=globe)

    dtext = ('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '))

    fig = plt.figure(1, figsize=(15, 12), dpi=300)

    ax1 = fig.add_subplot(1, 1, 1, projection=geos)

    ax1.imshow(RGB, origin='upper',
                       extent=(x.min(), x.max(), y.min(), y.max()),
                       transform=geos,
                       interpolation='nearest', vmin=162., vmax=330.)
    ax1.coastlines(resolution='10m', color='black', linewidth=1)
    #ax1.add_feature(ccrs.cartopy.feature.STATES)

    ax1.set_title('GOES-16 Nighttime Microphysics RGB product for Fog Detection', loc='left', fontweight='semibold', fontsize=15)
    ax1.set_title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');

    ax1.set_extent([-85, -80, 21, 24], crs=ccrs.PlateCarree())

    plt.savefig('fog_rgb1_'+dtext+'.png',bbox_inches='tight',dpi=300)
    plt.clf()
    plt.cla()
    plt.close('all')


#    fig = plt.figure(2, figsize=(15, 12), dpi=300)

#    ax2 = fig.add_subplot(1, 1, 1, projection=geos)

#    ax2.imshow(RGB_true, origin='upper',
#                       extent=(x.min(), x.max(), y.min(), y.max()),
#                       transform=geos,
#                       interpolation='nearest', vmin=162., vmax=330.)
#    ax2.coastlines(resolution='10m', color='black', linewidth=1)
#    #ax2.add_feature(ccrs.cartopy.feature.STATES)

#    ax2.set_title('GOES-16 True Color RGB product', loc='left', fontweight='semibold', fontsize=15)
#    ax2.set_title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');

#    ax2.set_extent([-85, -80, 21, 24], crs=ccrs.PlateCarree())

#    plt.savefig('truecolor_rgb_'+dtext+'.png',bbox_inches='tight',dpi=300)
#    plt.clf()
#    plt.cla()
#    plt.close('all')

    #plt.show()



