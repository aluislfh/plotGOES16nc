from datetime import datetime

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy  # noqa: F401
import numpy as np
import xarray
import os, sys, glob


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
    midpoint = str(C['t'].data)[:-8]
    scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

    print('Scan Start    : {}'.format(scan_start))
    print('Scan midpoint : {}'.format(scan_mid))
    print('Scan End      : {}'.format(scan_end))
    print('File Created  : {}'.format(file_created))
    print('Scan Duration : {:.2f} minutes'.format((scan_end-scan_start).seconds/60))

    # We'll use the `CMI_C02` variable as a 'hook' to get the CF metadata.
    dat = C.metpy.parse_cf('CMI_C02')

    geos = dat.metpy.cartopy_crs

    # We also need the x (north/south) and y (east/west) axis sweep of the ABI data
    x = dat.x
    y = dat.y
    

    # Create the RGB like we did before

    # Load the three channels into appropriate R, G, and B
    R = C['CMI_C02'].data
    G = C['CMI_C03'].data
    B = C['CMI_C01'].data

    # Apply range limits for each channel. RGB values must be between 0 and 1
    R = np.clip(R, 0, 1)
    G = np.clip(G, 0, 1)
    B = np.clip(B, 0, 1)

    # Apply the gamma correction
    gamma = 2.2
    R = np.power(R, 1/gamma)
    G = np.power(G, 1/gamma)
    B = np.power(B, 1/gamma)

    # Calculate the "True" Green
    G_true = 0.45 * R + 0.1 * G + 0.45 * B
    G_true = np.clip(G_true, 0, 1)

    # The final RGB array :)
    RGB = np.dstack([R, G_true, B])





    # Apply the normalization...
    cleanIR = C['CMI_C13'].data

    # Normalize the channel between a range.
    #       cleanIR = (cleanIR-minimumValue)/(maximumValue-minimumValue)
    cleanIR = (cleanIR-90)/(313-90)

    # Apply range limits to make sure values are between 0 and 1
    cleanIR = np.clip(cleanIR, 0, 1)

    # Invert colors so that cold clouds are white
    cleanIR = 1 - cleanIR

    # Lessen the brightness of the coldest clouds so they don't appear so bright
    # when we overlay it on the true color image.
    cleanIR = cleanIR/1.4

    # Yes, we still need 3 channels as RGB values. This will be a grey image.
    RGB_cleanIR = np.dstack([cleanIR, cleanIR, cleanIR])



#    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

#    ax1.set_title('True Color', fontweight='bold')
#    ax1.imshow(RGB)
#    ax1.axis('off')

#    ax2.set_title('Clean IR', fontweight='bold')
#    ax2.imshow(RGB_cleanIR)
#    ax2.axis('off')



    # Maximize the RGB values between the True Color Image and Clean IR image
    RGB_ColorIR = np.dstack([np.maximum(R, cleanIR), np.maximum(G_true, cleanIR),
                             np.maximum(B, cleanIR)])


#    fig = plt.figure(figsize=(15, 12))

#    ax = fig.add_subplot(1, 1, 1, projection=geos)

#    ax.imshow(RGB_ColorIR, origin='upper',
#              extent=(x.min(), x.max(), y.min(), y.max()),
#              transform=geos)

#    ax.coastlines(resolution='50m', color='black', linewidth=2)
#    ax.add_feature(ccrs.cartopy.feature.STATES)

#    plt.title('GOES-16 True Color and Night IR', loc='left', fontweight='bold',
#              fontsize=15)
#    plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y'), loc='right'),
#              loc='right')


    def contrast_correction(color, contrast):
        """Modify the contrast of an RGB.
        See:
        https://www.dfstudios.co.uk/articles/programming/image-programming-algorithms/image-processing-algorithms-part-5-contrast-adjustment/

        Input:
            color    - an array representing the R, G, and/or B channel
            contrast - contrast correction level
        """
        F = (259*(contrast + 255))/(255.*259-contrast)
        COLOR = F*(color-.5)+.5
        COLOR = np.clip(COLOR, 0, 1)  # Force value limits 0 through 1.
        return COLOR


    # Amount of contrast
    contrast_amount = 105

    # Apply contrast correction
    RGB_contrast = contrast_correction(RGB, contrast_amount)

    # Add in clean IR to the contrast-corrected True Color image
    RGB_contrast_IR = np.dstack([np.maximum(RGB_contrast[:, :, 0], cleanIR),
                                 np.maximum(RGB_contrast[:, :, 1], cleanIR),
                                 np.maximum(RGB_contrast[:, :, 2], cleanIR)])


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


    globe = ccrs.Globe(semimajor_axis=semi_major, semiminor_axis=semi_minor)
    geos = ccrs.Geostationary(central_longitude=sat_lon, 
                             satellite_height=sat_h, globe=globe)

    dtext = ('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '))


    # Plot on map with Cartopy

    fig = plt.figure(1, figsize=(15, 12), dpi=300)

    ax1 = fig.add_subplot(1, 1, 1, projection=geos)

    ax1.imshow(RGB_ColorIR, origin='upper',
               extent=(x.min(), x.max(), y.min(), y.max()),
               transform=geos)
    ax1.coastlines(resolution='10m', color='black', linewidth=2)
    ax1.add_feature(ccrs.cartopy.feature.BORDERS)
    ax1.set_title('True Color and Night IR')

    ax1.set_extent([-85, -80, 21, 24], crs=ccrs.PlateCarree())

    plt.savefig('truecolor1_rgb_'+dtext+'.png',bbox_inches='tight',dpi=300)
    plt.clf()
    plt.cla()
    plt.close('all')
    
    fig = plt.figure(2, figsize=(15, 12), dpi=300)

    ax2 = fig.add_subplot(1, 1, 1, projection=geos)

    ax2.imshow(RGB_contrast_IR, origin='upper',
               extent=(x.min(), x.max(), y.min(), y.max()),
               transform=geos)
    ax2.coastlines(resolution='10m', color='black', linewidth=2)
    ax2.add_feature(ccrs.cartopy.feature.BORDERS)
    ax2.set_title('Contrast Correction = {}'.format(contrast_amount))

    ax2.set_extent([-85, -80, 21, 24], crs=ccrs.PlateCarree())

    plt.subplots_adjust(wspace=.02)

    plt.savefig('IR_13_corrected_'+dtext+'.png',bbox_inches='tight',dpi=300)
    plt.clf()
    plt.cla()
    plt.close('all')
    
#    plt.show()













