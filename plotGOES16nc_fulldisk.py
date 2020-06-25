import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit
import numpy as np # Import the Numpy package
from datetime import datetime, date
from pyproj import Proj
from multiprocessing import Pool
from cpt_convert import loadCPT # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap # Linear interpolation for color maps
import os,glob,sys
import multiprocessing
import time
import datetime

def main():
    ruta = "/home/adrian/Desktop/python3"
    os.chdir(ruta)

    PROCESS_LIMIT=8

#-------------------------------------------------------------------------------

    # GOES 16 - Channel 14
    lista = glob.glob("OR_ABI-L2-CMIPF-M3C14_G16_*.nc")
    lista.sort()

    for i in lista:

        process=multiprocessing.Process(target=dataproc,args=(ruta,i,'IR4AVHRR6.cpt'))
        while(len(multiprocessing.active_children()) == PROCESS_LIMIT):
            time.sleep(1)
        process.start()


#-------------------------------------------------------------------------------

    # GOES 16 - Channel 14
    lista = glob.glob("OR_ABI-L2-CMIPF-M3C09_G16_*.nc")
    lista.sort()

    for i in lista:

        process=multiprocessing.Process(target=dataproc,args=(ruta,i,'WVCOLOR35.cpt'))
        while(len(multiprocessing.active_children()) == PROCESS_LIMIT):
            time.sleep(1)
        process.start()


#-------------------------------------------------------------------------------

    # GOES 16 - Channel 14
    lista = glob.glob("OR_ABI-L2-CMIPF-M3C08_G16_*.nc")
    lista.sort()

    for i in lista:

        process=multiprocessing.Process(target=dataproc,args=(ruta,i,'WVCOLOR35.cpt'))
        while(len(multiprocessing.active_children()) == PROCESS_LIMIT):
            time.sleep(1)
        process.start()


def dataproc(ruta,path,paleta):

    # Open NC file
    nc = Dataset(path)

    # Load data
    data_subset = nc.variables['CMI'][:][:,:]

    # Create the projection variables
    ori_proj = nc.variables['goes_imager_projection']
    sat_h = ori_proj.perspective_point_height
    sat_lon = ori_proj.longitude_of_projection_origin
    sat_sweep = ori_proj.sweep_angle_axis
    # The projection x and y coordinates equals
    # the scanning angle (in radians) multiplied by the satellite height (http://proj4.org/projections/geos.html)
    X = nc.variables['x'][:] * sat_h
    Y = nc.variables['y'][:] * sat_h
    p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
    # Convert map points to latitude and longitude with the magic provided by Pyproj
    XX, YY = np.meshgrid(X, Y)
    lons, lats = p(XX, YY, inverse=True)

    ii=[]
    jj=[]
    for i in range(lons.shape[1]):
        for j in range(lats.shape[0]):
            if lons[1515,i]>=-100.0 and lats[j,2325]>=15.0 and lons[1515,i]<=-65.0 and lats[j,2325]<=40.0:
                ii.append(i)
                jj.append(j)

    jmin,jmax = (int(np.min(ii)),int(np.max(ii)))
    imin,imax = (int(np.min(jj)),int(np.max(jj)))


    # Create map
    DPI = 150
    ax = plt.figure(figsize=(2000/float(DPI), 2000/float(DPI)), frameon=True, dpi=DPI)
    bmap = Basemap(projection='cyl', llcrnrlon=lons[1515,jmin], llcrnrlat=lats[imax,2325], urcrnrlon=lons[1515,jmax], urcrnrlat=lats[imin,2325],  resolution='l')

    bmap.drawparallels(np.arange(-90.0, 90.0, 10.0), linewidth=0.1, color='black', labels=[True, False, False, True])
    bmap.drawmeridians(np.arange(0.0, 360.0, 10.0), linewidth=0.1, color='black', labels=[True, False, False, True])

    # Load Shapefile
#    bmap.readshapefile('Nueva_DPA','Nueva_DPA',linewidth=0.90,color='darkslategray')
#    bmap.readshapefile('ne_10m_admin_0_countries','ne_10m_admin_0_countries',linewidth=0.90,color='darkslategray')
    bmap.drawstates(color='gray', linewidth=0.25)
    bmap.drawcoastlines(color='k', linewidth=0.9)
    bmap.drawcountries(color='k', linewidth=0.9)

    # Converts a CPT file to be used in Python
    cpt = loadCPT(paleta)
    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt)
     
    # Plot the GOES-16 channel with the converted CPT colors
    #bmap.pcolormesh(lons[imin:imax,jmin:jmax],lats[imin:imax,jmin:jmax],data_subset[imin:imax,jmin:jmax], cmap=cpt_convert, vmin=170, vmax=378)
    bmap.pcolormesh(lons[imin:imax,jmin:jmax],lats[imin:imax,jmin:jmax],data_subset[imin:imax,jmin:jmax]-273.15, cmap=cpt_convert, vmin=-103, vmax=104)

    # Search for the GOES-16 channel in the file name
    Band = (path[path.find("M3C")+3:path.find("_G16")])
    # Search for the Scan start in the file name
    Start = (path[path.find("_s")+2:path.find("_e")])
    Start_Formatted = Start[0:4] + " Day " + Start[4:7] + " - " + Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + "." + Start [13:14] + " UTC"
    # Search for the Scan end in the file name
    End = (path[path.find("_e")+2:path.find("_c")])
    End_Formatted = End[0:4] + " Day " + End[4:7] + " - " + End [7:9] + ":" + End [9:11] + ":" + End [11:13] + "." + End [13:14] + " UTC"
     
    # Add a title to the plot
    plt.title("GOES-16 ABI Band " + Band + "\n Scan from " + Start_Formatted + " to " + End_Formatted,fontsize=10)

    #bmap.colorbar(location='right', label='Brightness Temperature [K]')
    bmap.colorbar(location='right', label='Brightness Temperature [C]')


    # Save the result
    # Converting from julian day to dd-mm-yyyy
    year = int(Start[0:4])
    dayjulian = int(Start[4:7]) - 1 # Subtract 1 because the year starts at "0"
    dayconventional = datetime.datetime(year,1,1) + datetime.timedelta(dayjulian) # Convert from julian to conventional
    date = dayconventional.strftime('%d-%b-%Y') # Format the date according to the strftime directives
    date_save = dayconventional.strftime('%d%m%Y')
    time_save = Start [7:9] + Start [9:11] + Start [11:13]
    plt.savefig('G16_C' + str(Band) + '_' + date_save + '_' + time_save + '_fulldisk.png', dpi=DPI, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.cla()


main()




























