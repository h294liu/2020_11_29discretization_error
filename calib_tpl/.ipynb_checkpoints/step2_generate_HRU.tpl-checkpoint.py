#!/usr/bin/env python
# coding: utf-8
"""
Created on Fri Oct 16 09:35:11 2020
@author: hongli
"""

# ### 0. Import libraries ###
import os
import numpy as np
import rasterio as rio
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt 
import geospatial_functions.geospatial as gs

def read_raster(raster_file, var_name, area_name):
    with rio.open(raster_file) as ff:
        data_value  = ff.read(1)
        data_mask = ff.read_masks(1)
        gt = ff.transform

    dx = gt[0]
    dy = -gt[4]
    pixelArea = dx*dy

    (unique, counts) = np.unique(data_value[data_mask!=0], return_counts=True)
    frequencies = np.asarray((unique, counts)).T
    frequencies[:,1]=frequencies[:,1]*pixelArea

    df = pd.DataFrame(frequencies, columns = [var_name, area_name]) 
    return df 

def get_cdf_error(raw_df, dis_df, var_column, area_column):

    # get raw grid cell and hru numbers
    raw_num = len(raw_df)
    hru_num = len(dis_df)

    # sort dataframe based on interested variable values
    raw_df_sort = raw_df.sort_values(by=[var_column])
    dis_df_sort = dis_df.sort_values(by=[var_column])
    
    # create evenly spaced variable values 
    line_var_num = 200
    line_var = np.reshape(np.linspace(0.9*np.minimum(np.amin(raw_df_sort[var_column].values), np.amin(dis_df_sort[var_column].values)),
                    1.1*np.maximum(np.amax(raw_df_sort[var_column].values), np.amax(dis_df_sort[var_column].values)),
                    num=line_var_num),[line_var_num,1])
    
    # calculate variable cdf in raw_cdf
    line_var_rep = np.repeat(line_var, raw_num, axis=1) # shape [line_var_num, raw_num]
    raw_var_rep = np.repeat(np.reshape(raw_df_sort[var_column].values, [1, raw_num]), line_var_num, axis=0) 
    raw_area_rep =  np.repeat(np.reshape(raw_df_sort[area_column].values, [1, raw_num]), line_var_num, axis=0)
    
    condition_area = np.where(line_var_rep>raw_var_rep, raw_area_rep, 0.0)
    condition_area_cumsum = np.sum(condition_area, axis=1)
    total_area = np.sum(raw_area_rep, axis=1)
    line_var_cdf_in_raw = np.divide(condition_area_cumsum,total_area)
    del line_var_rep, raw_var_rep, raw_area_rep, condition_area, condition_area_cumsum, total_area
    
    # calculate variable in hru_cdf
    line_var_rep = np.repeat(line_var, hru_num, axis=1) # shape [line_var_num, hru_num]
    dis_var_rep = np.repeat(np.reshape(dis_df_sort[var_column].values, [1, hru_num]), line_var_num, axis=0) 
    hru_area_rep =  np.repeat(np.reshape(dis_df_sort[area_column].values, [1, hru_num]), line_var_num, axis=0)
    
    condition_area = np.where(line_var_rep>dis_var_rep, hru_area_rep, 0.0) # Area array with elements from hru_area_rep where condition is True, and 0 elsewhere.
    condition_area_cumsum = np.sum(condition_area, axis=1)
    total_area = np.sum(hru_area_rep, axis=1)
    line_var_cdf_in_hru = np.divide(condition_area_cumsum,total_area)
    del line_var_rep, dis_var_rep, hru_area_rep, condition_area, condition_area_cumsum, total_area
    
    # calucalte CDF difference
    cdf_dif = np.abs(line_var_cdf_in_raw-line_var_cdf_in_hru)
    line_var = np.reshape(line_var, np.shape(cdf_dif))
    var_error = np.trapz(cdf_dif, x=line_var)

    # relative error (normalized Wasserstein distance)     
    var_base = np.trapz(line_var_cdf_in_raw, x=line_var)    
    var_error_relative = var_error/var_base
    
    #  Kolmogorov–Smirnov statistic
    cdf_dif_max = max(cdf_dif)
    return var_error_relative, cdf_dif_max


# ### 0. Specify file paths ###
# --- common source file paths ---
root_dir = 'ROOT_DIR'

# --- case study depedent file paths ----
case_dir = root_dir #os.path.join(root_dir, case)

buf_dist = 600 # meter
    
dem_crop = os.path.join(case_dir, 'dem_crop.tif')
dem_crop_buf = os.path.join(case_dir, 'dem_crop_buf.tif')
dem_class_raster = os.path.join(case_dir, 'dem_class.tif')
dem_value_raster = os.path.join(case_dir, 'dem_value.tif')

slp_crop_buf = os.path.join(case_dir, 'slope_crop_buf.tif')
slp_crop = os.path.join(case_dir, 'slope_crop.tif')
slp_class_raster = os.path.join(case_dir, 'slope_class.tif')
slp_value_raster = os.path.join(case_dir, 'slope_value.tif')

asp_crop_buf = os.path.join(case_dir, 'aspect_crop_buf.tif')
asp_crop = os.path.join(case_dir, 'aspect_crop.tif')
asp_crop_180 = os.path.join(case_dir, 'aspect_crop_180.tif')
asp_class_raster = os.path.join(case_dir, 'aspect_class.tif')
asp_value_raster = os.path.join(case_dir, 'aspect_value.tif')

lc_crop = os.path.join(case_dir, 'landcover_crop.tif')
lc_crop_resample = os.path.join(case_dir, 'landcover_crop_resample.tif')
lc_class_raster = os.path.join(case_dir, 'landcover_class.tif')

sx_raster = os.path.join(case_dir, 'sx.tif')
sw_raster = os.path.join(case_dir, 'sw.tif')

sub_shp_prj = os.path.join(case_dir, 'subbasin_prj.shp')
sub_shp_prj_buf = os.path.join(case_dir, 'subbasin_prj_buf.shp')
sub_raster = os.path.join(case_dir, 'subbasin.tif')
sub_corr_txt=os.path.join(case_dir, 'subNo_HUC12_corr.txt') #correspondence between HUC12 and subbasin number.
 
hru_raster = os.path.join(case_dir, 'hru.tif')
hru_vector = os.path.join(case_dir, 'hru.shp')
hru_raster_diss = os.path.join(case_dir, 'hru_diss.tif')
hru_vector_diss = os.path.join(case_dir, 'hru_diss.shp')

hru_attrb_elev = os.path.join(case_dir, 'hru_attrb_elevation.tif') # attribute raster of HRU
hru_attrb_lc = os.path.join(case_dir, 'hru_attrb_landcover.tif')
hru_attrb_slp = os.path.join(case_dir, 'hru_attrb_slope.tif')
hru_attrb_asp = os.path.join(case_dir, 'hru_attrb_aspect.tif')
hru_attrb_sw_mean = os.path.join(case_dir, 'hru_attrb_sw_mean.tif')
hru_attrb_sw_basename = os.path.join(case_dir, 'hru_attrb_sw_DOY') # e.g., hru_attrb_sw_DOY100.tif

error_file = os.path.join(case_dir, 'Diagnostics.txt')
hist_ofile = os.path.join(case_dir, 'Sw_hist.png')
cdf_ofile = os.path.join(case_dir, 'Sw_cdf.png')

# -- define GRU and HRU field names and data types --
subNo_field = 'GRUNo'
subNo_field_type = 'int32'
subName_field = 'GRUId'

hruNo_field = 'HRUNo'
hruNo_field_type = 'int32'
hruName_field = 'HRUId'

area_field = 'areaSqm'

# --- define common projection, nodata value, reference raster----
proj4="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs"
dst_crs = rio.crs.CRS.from_string(proj4)
## Albers Equal Area Conic Projection. 
## reference 1: https://gisgeography.com/conic-projection-lambert-albers-polyconic/
## reference 2: https://epsg.io/102008
nodatavalue = -9999  
refraster = dem_crop # reference raster to rasterize vector and resample

# HRU elimination threshold
hru_area_thld = 10**6 #1km2


print('PART 1. classify slope raster')
bins = [0,slp_thrsh1_value,slp_thrsh2_value,slp_thrsh3_value,90] # including the rightmost edge
bins.sort()
slp_classes = len(bins)-1
gs.classify_elevation(slp_crop, sub_raster, bins, slp_class_raster, slp_value_raster, nodatavalue)

print('PART 2. classify aspect raster')
bins = [-1,0,asp_thrsh1_value,asp_thrsh2_value,asp_thrsh3_value,181]
bins.sort()
asp_classes = len(bins)-1
gs.classify_elevation(asp_crop_180, sub_raster, bins, asp_class_raster, asp_value_raster, nodatavalue)

print('PART 3. genearte HRU')
raster_list = [sub_raster, dem_class_raster, lc_class_raster, slp_class_raster, asp_class_raster]
raster_fieldname_list = [subNo_field, 'elevClass', 'lcClass', 'slpClass', 'aspClass']

gs.define_hru(raster_list, raster_fieldname_list, sub_raster, sub_corr_txt, subNo_field, subName_field,
              nodatavalue, hru_raster, hru_vector, hruNo_field, hruNo_field_type, hruName_field)

print('PART 4. zonal area')
in_gpd = gpd.read_file(hru_vector)
in_gpd[area_field] = in_gpd.area
in_gpd.to_file(hru_vector)

print('PART 5. Sw zonal statistics')
gs.zonal_statistic(sw_raster, hru_vector, hruNo_field, hruNo_field_type, refraster, 'mean', 
                   hru_attrb_sw_mean, nodatavalue, output_column_prefix='sw')

print('PART 6. Calculate errors')
sw_var = 'Sw'
sx_var = 'Sx'
area_var = 'area_m'

# read sw
raw_df_sw = read_raster(sw_raster, sw_var, area_var)
dis_df_sw = read_raster(hru_attrb_sw_mean, sw_var, area_var)

# calculate cdf error
sw_error, sw_cdf_dif_max = get_cdf_error(raw_df_sw, dis_df_sw, sw_var, area_var)

# output Diagnostics
f = open(error_file, 'w')
f.write('HRU_number, Sw_error, SW_cdf_dif_max\n')
f.write('%d, ' % len(dis_df_sw))
f.write('%.6f, ' % sw_error)
f.write('%.6f, ' % sw_cdf_dif_max)
f.close() 

print('Done')



