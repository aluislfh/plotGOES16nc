from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit
import numpy as np # Import the Numpy package
from datetime import datetime
from pyproj import Proj
import os,glob
import multiprocessing
import time
import datetime
import netCDF4

def main():

    odir = "/home/adrian/Desktop/G16_crop/nc_output"

    wdir = "/home/adrian/Desktop/G16_crop"
    os.chdir(wdir)

    PROCESS_LIMIT=6

#-------------------------------------------------------------------------------

    # GOES 16 - Channel 14
    glist = glob.glob("OR_ABI-L2-CMIPC-M6C14_G16_*.nc")
    glist.sort()

    for gfile in glist:

        process=multiprocessing.Process(target=dataproc,args=(wdir,odir,gfile))
        while(len(multiprocessing.active_children()) == PROCESS_LIMIT):
            time.sleep(1)
        process.start()


def dataproc(wdir,odir,path):

    # Open NC file
    nc = netCDF4.Dataset(path)

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

    # Select points of the grid
    ii=[]
    jj=[]
    for i in range(lons.shape[1]):
        for j in range(lats.shape[0]):
            if lons[-1,i]>=-86.0 and lats[j,-1]>=19.0 and lons[-1,i]<=-73.0 and lats[j,-1]<=25.0:
                ii.append(i)
                jj.append(j)

    # d03 -85.71704  -73.765396 19.340546 24.264412 1 28 183 411

    jmin,jmax = (int(np.min(ii)),int(np.max(ii)))
    imin,imax = (int(np.min(jj)),int(np.max(jj)))
     
    # Search for the GOES-16 channel in the file name
    Band = (path[path.find("M3C")+3:path.find("_G16")])
    # Search for the Scan start in the file name
    Start = (path[path.find("_s")+2:path.find("_e")])
    Start_Formatted = Start[0:4] + " Day " + Start[4:7] + " - " + Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + "." + Start [13:14] + " UTC"
    # Search for the Scan end in the file name
    End = (path[path.find("_e")+2:path.find("_c")])
    End_Formatted = End[0:4] + " Day " + End[4:7] + " - " + End [7:9] + ":" + End [9:11] + ":" + End [11:13] + "." + End [13:14] + " UTC"
     
    # Save the result
    # Converting from julian day to dd-mm-yyyy
    year = int(Start[0:4])
    dayjulian = int(Start[4:7]) - 1 # Subtract 1 because the year starts at "0"
    dayconventional = datetime.datetime(year,1,1) + datetime.timedelta(dayjulian) # Convert from julian to conventional
    date = dayconventional.strftime('%d-%b-%Y') # Format the date according to the strftime directives
    date_save = dayconventional.strftime('%d%m%Y')
    time_save = Start [7:9] + Start [9:11] + Start [11:13]
    otime = date_save[6:8] + date_save[2:4] + date_save[:2] + Start[7:9] + Start[9:11] + Start[11:13]
    print(otime)

    # Escribir el archivo NetCDF
    write_netcdf(lons[imin:imax,jmin:jmax], lats[imin:imax,jmin:jmax], data_subset[imin:imax,jmin:jmax]-273.15, otime, odir+'/G16_C' + str(Band) + '_' + date_save + '_' + time_save + '_subregion.nc')

    nc.close()

    os.system('rm -f ' + wdir + '/' + path)


def write_netcdf(xlon, xlat, idata, rdate, outfilename):

    iyy, imm, idd, ihh, imn, iss = ('20'+rdate[:2], rdate[2:4], rdate[4:6], rdate[6:8] , rdate[8:10], rdate[10:12])
    print(iyy+'-'+imm+'-'+idd+' '+ihh+':'+imn+':'+iss)
    
    dates = [datetime.datetime(int(iyy), int(imm), int(idd), int(ihh), int(imn), int(iss))]


    # Create the new netCDF file
    fid = netCDF4.Dataset(outfilename,'w', format='NETCDF4_CLASSIC')

    # Define the dimensions
    time = fid.createDimension('time', len(dates))
    lon  = fid.createDimension('lon', idata.shape[1])
    lat  = fid.createDimension('lat', idata.shape[0])

    # Create global attributes
    fid.title          = 'GOES-16 CH14 --> '+str(dates[0])
    fid.description    = 'GOES-16 level 2 data subregion'
    fid.institution    = "CFA/INSMET"
    fid.acknowledgment = "AWS NASA SERVICE"

    # Create variable TIMES
    times = fid.createVariable('time', np.float64, ('time',))
    times.calendar = 'standard'
    times.units = 'hours since '+iyy+'-'+imm+'-'+idd+' '+ihh+':'+imn+':'+iss

    times[:] = netCDF4.date2num(dates, units = times.units, calendar = times.calendar)

    fid.variables['time'] = times[:]
    fid.variables['time'].standard_name='time'
    fid.variables['time'].units = 'hours since '+iyy+'-'+imm+'-'+idd+' '+ihh+':'+imn+':'+iss
    fid.variables['time'].comment='Time dimenssion'

    # Create variable LONGITUDES
    longitudes = fid.createVariable('lon', np.float32, ('lon'))
    fid.variables['lon'].standard_name='longitude'
    fid.variables['lon'].long_name='longitude'
    fid.variables['lon'].units='degrees_east'
    fid.variables['lon'].comment='Longitude'
    fid.variables['lon'].axis='X'
    fid.variables['lon'][:] = xlon[0,:]

    # Create variable LATITUDES
    latitudes = fid.createVariable('lat', np.float32, ('lat'))
    fid.variables['lat'].standard_name='latitude'
    fid.variables['lat'].long_name='latitude'
    fid.variables['lat'].units='degrees_north'
    fid.variables['lat'].comment='Latitude'
    fid.variables['lat'].axis='Y'
    fid.variables['lat'][:] = xlat[:,0]

    # Create variable VAR
    nc_var = fid.createVariable('xlon', np.float32,('time','lat', 'lon')) 

    fid.variables['xlon'][0,:,:] = xlon[:,:]
    fid.variables['xlon'].units='degrees_east'
    fid.variables['xlon'].missing_value=-999.0
    fid.variables['xlon'].description='Longitudes'

    # Create variable VAR
    nc_var = fid.createVariable('xlat', np.float32,('time','lat', 'lon')) 

    fid.variables['xlat'][0,:,:] = xlat[:,:]
    fid.variables['xlat'].units='degrees_north'
    fid.variables['xlat'].missing_value=-999.0
    fid.variables['xlat'].description='Latitudes'

    # Create variable VAR
    nc_var = fid.createVariable('g16ch14', np.float32,('time','lat', 'lon')) 

    fid.variables['g16ch14'][0,:,:] = idata[:,:]
    fid.variables['g16ch14'].units='Celsius'
    fid.variables['g16ch14'].missing_value=-999.0
    fid.variables['g16ch14'].description='Brightness Temperature'

    # Closing file
    fid.close()



main()




























