#!/usr/bin/env python
# coding: utf-8

from osgeo import gdal
import argparse, os
import numpy as np
from pyproj import Proj
gdal.UseExceptions()

def calculate_Sx(raw_y_coords_2d,raw_x_coords_2d,raw_ele,raw_mask,buf_y_coords_2d,buf_x_coords_2d,buf_ele,nodata):
    
    # raw data are in shape (ny_subset,nx_subset) 
    # buf data are in shape (ny_buf_subset, nx_buf_subset)

    # define upwind window
    R = 6371.0 #km
#     cell_size = 9.259259259300000407e-05 #0.000833333 # grid cell size (degree) obtained from raw DEM
#     dmax = 0.5 #km 0.1, 0.2, 0.5
#     azims_wind = range(235, 340+1, 5)

    cell_size = 0.0008333333333333467 # grid cell size (degree) obtained from raw DEM
    dmax = 0.1 #km 0.1
    azims_wind = range(190, 280+1, 5)

    (ny,nx) = np.shape(raw_y_coords_2d)
    raw_num = ny*nx
    (ny_buf,nx_buf) = np.shape(buf_y_coords_2d)
    buf_num = ny_buf*nx_buf

    # reshape buf array (raw_num, buf_num). Repeat on rows.
    print('(1) process buf data')
#     p = Proj(proj='utm',zone=13,ellps='GRS80',datum='NAD83', preserve_units=False)
    p = Proj("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")
    buf_lon, buf_lat = p(buf_x_coords_2d, buf_y_coords_2d, inverse=True) # inverse transform
    buf_lat = np.deg2rad(buf_lat.reshape(1,buf_num).repeat(raw_num, axis=0))
    buf_lon = np.deg2rad(buf_lon.reshape(1,buf_num).repeat(raw_num, axis=0))
    buf_ele = (buf_ele.reshape(1,buf_num).repeat(raw_num, axis=0))*0.001 #m to km

    # reshape raw array (raw_num, buf_num). Repeat on columns.
    print('(2) process raw data')
#     p = Proj(proj='utm',zone=13,ellps='GRS80',datum='NAD83', preserve_units=False)
    p = Proj("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")
    raw_lon, raw_lat = p(raw_x_coords_2d, raw_y_coords_2d, inverse=True) # inverse transform
    raw_lat = np.deg2rad(raw_lat.reshape(raw_num,1).repeat(buf_num, axis=1))
    raw_lon = np.deg2rad(raw_lon.reshape(raw_num,1).repeat(buf_num, axis=1))
    raw_ele = (raw_ele.reshape(raw_num,1).repeat(buf_num, axis=1))*0.001
    delta_ele = buf_ele - raw_ele
    raw_mask = raw_mask.reshape(raw_num,1)

    # Haversine distance (km)
    print('(3) calculate azim')
    a = (np.sin((raw_lat - buf_lat)/2))**2 + np.cos(raw_lat)*np.cos(buf_lat) * (np.sin((raw_lon - buf_lon)/2))**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    d = R * c 

    # 4-corner Azimuth (radiance) from the interested cell to a buffered cell
    buf_lat_top = buf_lat + 0.5*np.deg2rad(cell_size) 
    buf_lat_bottom = buf_lat - 0.5*np.deg2rad(cell_size) 
    buf_lon_left = buf_lon - 0.5*np.deg2rad(cell_size) 
    buf_lon_right = buf_lon + 0.5*np.deg2rad(cell_size) 

    theta_top_left = azimuth(theta1 = raw_lat, lamda1 = raw_lon, theta2 = buf_lat_top, lamda2 = buf_lon_left)
    theta_top_right = azimuth(theta1 = raw_lat, lamda1 = raw_lon, theta2 = buf_lat_top, lamda2 = buf_lon_right)
    theta_bottom_left = azimuth(theta1 = raw_lat, lamda1 = raw_lon, theta2 = buf_lat_bottom, lamda2 = buf_lon_left)
    theta_bottom_right = azimuth(theta1 = raw_lat, lamda1 = raw_lon, theta2 = buf_lat_bottom, lamda2 = buf_lon_right)
    del raw_lat,raw_ele,buf_lat,buf_ele

    azim = np.zeros((raw_num,buf_num,4),dtype=float)
    azim[:,:,0] = theta_top_left
    azim[:,:,1] = theta_top_right
    azim[:,:,2] = theta_bottom_left
    azim[:,:,3] = theta_bottom_right

    azim_min = np.nanmin(azim, axis=2) # shape(raw_num,buf_num)
    azim_max = np.nanmax(azim, axis=2)

    # Calculate SX along every increment wind vector
    print('(4) calculate Sx')
    wind_num=len(azims_wind)
    Sx_azim = np.zeros((raw_num,wind_num)) # each grid's Max Sx per wind direction

    for i in range(raw_num):
        if raw_mask[i]==False:
            for j in range(wind_num):
                azim_wind = azims_wind[j]
                azim_wind_rad = azim_wind*np.pi/180.0

                # Select buffered cells containing the wind vector
                condition1 = d[i,:] <= dmax # max distance from cell of interest
                condition2 = np.logical_and(azim_min[i,:] <= azim_wind_rad, azim_max[i,:] >= azim_wind_rad) # cover wind vector
                condition3 = buf_lon_right[i,:] <= raw_lon[i,0] # on the west of cell of interest
                select_col = np.where(np.logical_and(np.logical_and(condition1, condition2), condition3) == True)

                if len(select_col[0])==0:
                    Sx_ij=np.nan
                else:
                    Sx_ij = np.max(np.arctan(np.tan(np.divide(delta_ele[i,select_col], d[i,select_col])))*180/np.pi)
                    if Sx_ij < 0.0: #[-pi, pi] -> [0, 2pi]
                        Sx_ij = Sx_ij + 2*np.pi 
                Sx_azim[i,j]=Sx_ij
        else:
            Sx_azim[i,:]=np.nan
    Sx = np.nanmean(Sx_azim, axis=1) # each grid's mean Max Sx of all wind directions.(raw_um,)
    Sx = np.reshape(Sx, (ny,nx))
    Sx_fill = np.where(~np.isnan(Sx),Sx,nodata)
    return Sx_fill

def azimuth(theta1, lamda1, theta2, lamda2):
    # Azimuth (in radiance)
    # reference: https://www.omnicalculator.com/other/azimuth#how-to-calculate-the-azimuth-an-example
    delta_lamda = lamda2 - lamda1
    y = np.multiply(np.sin(delta_lamda),np.cos(theta2))
    x1 = np.multiply(np.cos(theta1), np.sin(theta2))
    x2 = np.multiply(np.multiply(np.sin(theta1), np.cos(theta2)), np.cos(delta_lamda))
    x = x1-x2
    theta = np.arctan2(y, x) # NOTE: In EXCEL, ATAN2(x,y). The places of x and y are tricky.    

    theta_2pi = np.where(theta >= 0.0, theta, 2*np.pi+theta) #[-pi, pi] -> [0, 2pi]
    return theta_2pi 

def process_command_line():
    '''Parse the commandline'''
    parser = argparse.ArgumentParser(description='Script to subset a netcdf file based on a list of IDs.')
    parser.add_argument('dem_raster', help='path of file with watershed DEM grid cells.')
    parser.add_argument('dem_buf_raster',help='path of file with buffered DEM grid cells.')
    parser.add_argument('row_start',help='start row id of rawfile.')
    parser.add_argument('row_end',help='end row id of rawfile.')
    parser.add_argument('col_start',help='start column id of rawfile.')
    parser.add_argument('col_end',help='end column id of rawfile.')
    parser.add_argument('buf_window',help='number of buffered grids for buffered dem extraction.')    
    parser.add_argument('opath',help='output folder path.')
    args = parser.parse_args()
    return(args)


# main script
# process command line
args = process_command_line()
dem_raster=args.dem_raster #'dem.tif'
dem_buf_raster=args.dem_buf_raster #'dem_buf.tif'
row_start=int(args.row_start)
row_end=int(args.row_end)
col_start=int(args.col_start)
col_end=int(args.col_end)
buf_window=int(args.buf_window)

outfolder = args.opath
if not os.path.exists(outfolder):
    os.makedirs(outfolder)
dem_nodata=-9999
Sx_nodata=-9999

#====================================================
# read raw raster [ELEVATION]
print('read raw data')
r = gdal.Open(dem_raster)
band = r.GetRasterBand(1) #bands start at one
elev = band.ReadAsArray().astype(np.float)
mask = (elev==dem_nodata)
elev = np.where(elev==dem_nodata,np.nan,elev)

elev_subset = elev[row_start:row_end,col_start:col_end]
mask_subset = mask[row_start:row_end,col_start:col_end]
del elev,mask

#====================================================
# calculate raw lat/lon coordinates
# reference: https://gis.stackexchange.com/questions/42790/gdal-and-python-how-to-get-coordinates-for-all-cells-having-a-specific-value
print('calculate raw lat/lon')
(upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = r.GetGeoTransform()
x_index=np.arange(col_start,col_end)
y_index=np.arange(row_start,row_end)
x_coords = x_index * x_size + upper_left_x + (x_size / 2.0) #add half the cell size to centre the point
y_coords = y_index * y_size + upper_left_y + (y_size / 2.0) 

ny = row_end-row_start
nx = col_end-col_start
x_coords_2d = np.repeat(np.reshape(x_coords,(1,nx)), ny, axis=0)
y_coords_2d = np.repeat(np.reshape(y_coords,(ny,1)), nx, axis=1)
del x_coords,y_coords,r

#====================================================
print('read buf data')
# read buffered raster [ELEVATION]
r = gdal.Open(dem_buf_raster)
band = r.GetRasterBand(1) #bands start at one
elev_buf = band.ReadAsArray().astype(np.float)
elev_buf = np.where(elev_buf==dem_nodata,np.nan,elev_buf)
(ny_buf,nx_buf) = np.shape(elev_buf)

# calculate buffered lat/lon coordinates
print('calculate buf lat/lon')
(upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = r.GetGeoTransform()
x_index_buf=np.arange(nx_buf)
y_index_buf=np.arange(ny_buf)
x_coords_buf = x_index_buf * x_size + upper_left_x + (x_size / 2.0) #add half the cell size to centre the point
y_coords_buf = y_index_buf * y_size + upper_left_y + (y_size / 2.0) 

x_coords_2d_buf = np.repeat(np.reshape(x_coords_buf,(1,nx_buf)), ny_buf, axis=0)
y_coords_2d_buf = np.repeat(np.reshape(y_coords_buf,(ny_buf,1)), nx_buf, axis=1)
del x_coords_buf,y_coords_buf,r

# subset
print('subset buf data')
buf_row_start = np.argmin(abs(y_coords_2d_buf[:,0] - (np.amax(y_coords_2d[:,0])+buf_window)))
buf_row_end = np.argmin(abs(y_coords_2d_buf[:,0] - (np.amin(y_coords_2d[:,0])-buf_window)))
buf_col_start = np.argmin(abs(x_coords_2d_buf[0,:] - (np.amin(x_coords_2d[0,:])-buf_window)))
buf_col_end = np.argmin(abs(x_coords_2d_buf[0,:] - (np.amax(x_coords_2d[0,:])+ buf_window)))

x_coords_2d_buf_subset = x_coords_2d_buf[buf_row_start:buf_row_end,buf_col_start:buf_col_end]
y_coords_2d_buf_subset = y_coords_2d_buf[buf_row_start:buf_row_end,buf_col_start:buf_col_end]
elev_buf_subset = elev_buf[buf_row_start:buf_row_end,buf_col_start:buf_col_end]   
del x_coords_2d_buf,y_coords_2d_buf,elev_buf

# #====================================================
# calcualte Sx
print('calcualte Sx')
Sx=calculate_Sx(y_coords_2d,x_coords_2d,elev_subset,mask_subset,
                y_coords_2d_buf_subset,x_coords_2d_buf_subset,elev_buf_subset,Sx_nodata)

# write Sx output
ofile='Sx_Row'+str(row_start)+'_'+str(row_end)+'_Col'+str(col_start)+'_'+str(col_end)+'.txt'
np.savetxt(os.path.join(outfolder,ofile),Sx,delimiter=',',fmt='%f')

print('Done')