#!/usr/bin/env python
# coding: utf-8

"""
Created on Fri Oct 16 09:35:11 2020

@author: hongli

Note: in python, there are two types of nodata mask for raster files.
The first is in the format of GDAL. 0 is invalid region, 255 is valid region. 
The second is in the foramt of Numpy masked array. True is invalid region, False is valid region.
When read mask use ff.read_masks(1), it by default returns mask in the first format.
(reference: https://rasterio.readthedocs.io/en/latest/topics/masks.html)
"""

import os, fiona
import numpy as np
import rasterio as rio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import geopandas as gpd
import fiona.crs 
import rasterio.mask
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt 

def reproject_raster(inraster, outraster, dst_crs):
    # reference: https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/reproject-raster/    
    with rio.open(inraster) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rio.open(outraster, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
    return
                
def reproject_vector(invector, outvector, new_crs):
    # reference:https://www.earthdatascience.org/courses/use-data-open-source-python/intro-vector-data-python/vector-data-processing/reproject-vector-data-in-python/
    in_gdf = gpd.read_file(invector)            # read vector file
    in_gdf_prj = in_gdf.to_crs(new_crs)         # convert projection   
    in_gdf_prj.to_file(outvector)               # save projected geodataframe
    return 

def buffer_vector(invector, outvector, buf_dist):
    # reference: https://stackoverflow.com/questions/51263138/how-to-create-an-accurate-buffer-of-5-miles-around-a-coordinate-in-python
    in_gdf = gpd.read_file(invector)                       # read vector file
    in_gdf['geometry'] = in_gdf.geometry.buffer(buf_dist)  # buffer (unit:meter)
    in_gdf.to_file(outvector)                              # save geodataframe    
    return 

def calculate_slope_aspect(dem_file, out_slope_raster, out_aspect_raster):
    # read dem file
    with rio.open(dem_file) as ff:
        dem  = ff.read(1)
        dem_mask = ff.read_masks(1)
        out_meta = ff.meta.copy()
        gt = ff.transform

    (nx,ny) = np.shape(dem)
    dx = gt[0]*np.ones_like(dem)
    dy = -gt[4]*np.ones_like(dem)
    slope_arr = dem.copy()
    aspect_arr = dem.copy()

    # calculate slope and aspect
    for i in range(nx):
        for j in range(ny):            
            # filter out masked grid cell
            if (dem_mask[i,j] != 0):
                imin = i-1
                imax = i+1
                jmin = j-1
                jmax = j+1

                # adjust boundary neighbors 
                if (imin < 0): imin = 0
                if (imax > nx-1): imax = nx-1
                if (jmin < 0): jmin = 0
                if (jmax > ny-1): jmax = ny-1

                # calculate dz/dy
                dzdy = ((dem[imin,jmin] + 2*dem[i,jmin] + dem[imax,jmin]) \
                      - (dem[imin,jmax] + 2*dem[i,jmax] + dem[imax,jmax])) \
                      /((dy[imin,jmin] + 2*dy[i,jmin] + dy[imax,jmin]) \
                      + (dy[imin,jmax] + 2*dy[i,jmax] + dy[imax,jmax]))
                # calcualte dz/dx
                dzdx = ((dem[imin,jmin] + 2*dem[imin,j] + dem[imin,jmax]) \
                      - (dem[imax,jmin] + 2*dem[imax,j] + dem[imax,jmax])) \
                      /((dx[imin,jmin] + 2*dx[imin,j] + dx[imin,jmax]) \
                      + (dx[imax,jmin] + 2*dx[imax,j] + dx[imax,jmax])) 

                # calculate slope
                rise_run = (dzdx**2 + dzdy**2)**0.5
                slope_degrees = np.arctan(rise_run) * 180/float(np.pi)

                # calculate aspect
                aspect = 180/float(np.pi) * np.arctan2(dzdy,-dzdx)
                # note special case where slope=0
                if slope_degrees == 0.0:
                    aspect_degrees = -1
                else:
                    if aspect < 0:
                        aspect_degrees = 90.0 - aspect  
                    elif aspect > 90.0:
                        aspect_degrees = 360.0 - aspect + 90.0
                    else:
                        aspect_degrees = 90.0 - aspect
                # assign
                slope_arr[i,j] = slope_degrees
                aspect_arr[i,j]= aspect_degrees

    # save slope and asepct rasters
    with rio.open(out_slope_raster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(slope_arr, dem_mask==0)
        outf.write(out_arr_ma,1)             
    with rio.open(out_aspect_raster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(aspect_arr, dem_mask==0)
        outf.write(out_arr_ma,1)       
    return

def convert_aspect_180(in_aspect_raster, out_aspect_raster):
    # read dem file
    with rio.open(in_aspect_raster) as ff:
        asp  = ff.read(1)
        asp_mask = ff.read_masks(1)
        out_meta = ff.meta.copy()
    # convert
    asp = np.where(asp>180, 360-asp, asp)
    # save
    with rio.open(out_aspect_raster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(asp, asp_mask==0)
        outf.write(out_arr_ma,1)       
    return
    
def rasterize_subbasin_vector(invector,infield,infield_rename,subNo_field,subNo_field_dtype,
                              refraster,outraster,sub_corr_txt,nodatavalue):
    # reference: https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
    '''
    invector: input, vector, the vector to be rasterized.
    infield: input, str, attribute field of input vector to use as a burn-in value.
    infield_rename: input, str, rename the attribute infield to be easy understand or commonly acceptable.
    subNo_field: input, str, subbasin number field name to be added to the updated invector.
    subNo_field_type: input, str, data type for subNo_field.
    refraster: input, raster, reference raster to get meta. 
    outraster: output, raster, path of output raster.
    sub_corr_txt: output, txt, subNo-infield_rename correspondence relationship.
    nodatavalue: input, int, use to fill no data grids.
    '''
    # open input vector
    in_gdf = gpd.read_file(invector)    

    # add a subNo column to gdf
    in_gdf[subNo_field] = np.arange(1,len(in_gdf)+1)

    # rename infield with HUC12
    in_gdf = in_gdf.rename(columns = {infield: infield_rename}) 

    # save the correspondence between field (HUC12) and subNo if field exists
    if infield_rename in in_gdf.columns:
        in_gdf[[subNo_field,infield_rename]].to_csv(sub_corr_txt, sep=',', index=False)    

    # copy and update the metadata from the input raster for the output
    with rio.open(refraster) as src:                 
        ref_mask = src.read_masks(1)
        meta = src.meta.copy()
    meta.update(count=1, dtype=subNo_field_dtype, compress='lzw', nodata=nodatavalue)

    with rio.open(outraster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # burn the features into the raster and write it out
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(in_gdf['geometry'], in_gdf[subNo_field]))  
        # Areas not covered by input geometries are replaced with an optional fill value, which defaults to 0.
        burned = rio.features.rasterize(shapes=shapes, out=out_arr, fill=nodatavalue, transform=out.transform)
        burned_ma = np.ma.masked_array(burned, ref_mask==0)
        out.write(burned_ma,1) 
    
    # save updated subbasin shapefile
    in_gdf.to_file(invector) 
    return

def rasterize_vector(invector,infield,infield_dtype,refraster,outraster,nodatavalue):
    # reference: https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
    '''
    invector: input, vector, the vector to be rasterized.
    infield: input, string. field of input vector to use as a burn-in value.
    infield_dtype: input, string. data type of infield, used to rasterize invector.
    refraster: input, raster, reference raster to get meta. 
    outraster: output, raster, path of output raster.
    nodatavalue: input, int, use to fill no data grids.
    '''
    # open input vector
    in_gdf = gpd.read_file(invector)    
        
    # copy and update the metadata from the refraster for the output
    with rio.open(refraster) as src:                 
        ref_mask = src.read_masks(1)
        meta = src.meta.copy()
    meta.update(count=1, dtype=infield_dtype, compress='lzw', nodata=nodatavalue) 
    
    with rio.open(outraster, 'w+', **meta) as out:
        out_arr = out.read(1)
    
        # burn the features into the raster and write it out
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(in_gdf['geometry'], in_gdf[infield]))  
        # Areas not covered by input geometries are replaced with an fill value.
        # Note that this process may induce pixels missing if the pixels' center is not within the polygon.
        # reference: https://rasterio.readthedocs.io/en/latest/api/rasterio.features.html
        burned = rio.features.rasterize(shapes=shapes, out=out_arr, fill=nodatavalue, transform=out.transform)
        burned_ma = np.ma.masked_array(burned, ref_mask==0)
        out.write(burned_ma,1)    
    return

def crop_raster(inraster,invector,outraster,nodatavalue):
    '''
    inraster: input, raster, rater to be cropped.
    invector: input, vector, provide crop extent.
    outraster: output, raster, output raster after crop.
    nodatavalue'''
    
    # open crop extent
    in_gdf = gpd.read_file(invector) 
    shapes = in_gdf['geometry']    
    
    # crop raster with the extent
    with rio.open(inraster) as src:
        if src.crs == in_gdf.crs:     # check projection consistency
            out_arr, out_transform = rio.mask.mask(src, shapes, crop=True, nodata=nodatavalue)
            out_meta = src.meta.copy()
     
    # update the new cropped meta info
    out_meta.update({"driver": "GTiff",
                     "height": out_arr.shape[1],
                     "width": out_arr.shape[2],
                     "transform": out_transform})
    
    # write cropped raster data
    with rio.open(outraster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(out_arr, out_arr==nodatavalue)
        outf.write(out_arr_ma)            
    return

def resample_raster(inraster,refraster,outraster):
    '''
    inraster: input, raster, rater to be cropped.
    refraster: input, raster, provide reference of array shape.
    outraster: output, raster, output raster after crop.'''
    
    # read refraster to get band size
    with rio.open(refraster) as ff:
        ref_crs = ff.crs
        ref_mask = ff.read_masks(1)
        ref_height, ref_width = ff.height, ff.width
        
    # method 1. use rasterio
    with rio.open(inraster) as ff:
        data  = ff.read(1)
        out_meta = ff.meta.copy()

        # reference: https://rasterio.readthedocs.io/en/stable/topics/resampling.html
        with rio.open(inraster) as ff:
            # resample data to target shape
            data = ff.read(
                out_shape=(int(ref_height),int(ref_width)),
                resampling=rio.enums.Resampling.nearest)

            # scale image transform
            transform = ff.transform * ff.transform.scale(
                (ff.width / data.shape[-1]),
                (ff.height / data.shape[-2]))

    # write resampled raster data
    out_meta.update(transform=transform, driver='GTiff', 
                    height=ref_height, width=ref_width, crs=ref_crs)

    with rio.open(outraster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(data, ref_mask==0)
        outf.write(out_arr_ma) 

    # method 2. use gdal 
    # reference: https://gis.stackexchange.com/questions/271226/up-sampling-increasing-resolution-raster-image-using-gdal
    return

def classify_elevation(dem_raster, base_raster, bins, dem_class_raster, dem_value_raster, nodatavalue):
    '''
    base_raster: input, raster, base as mask to classify elevation per base unit.
    dem_raster: input, raster, elevation raster file.
    bins: input, Three types:
        - If bins is an int, it defines the number of equal-width bins in the given range. 
        - If bins is a sequence, it defines a monotonically increasing array of bin edges, 
        including the rightmost edge, allowing for non-uniform bin widths.
        - If bins is a string, it defines the method used to calculate the bin edge sequence.
    dem_class_raster: output, raster, elevation class.
    dem_value_raster: output, raster, new elevation value of its class. 
    nodatavalue'''
    
    # read dem data
    with rio.open(dem_raster) as ff:
        dem  = ff.read(1)
        dem_mask = ff.read_masks(1)
        out_meta = ff.meta.copy()

    # read base and identify unique base
    with rio.open(base_raster) as ff:
        bases  = ff.read(1)
        bases_mask = ff.read_masks(1)
    unique_bases = np.unique(bases[bases_mask!=0])

    # define array in the same shape of dem, and specify dtype!
    elev_class = np.ones(np.shape(dem), dtype=np.int32)*int(nodatavalue)
    elev_value = np.ones(np.shape(dem), dtype=np.float32)*int(nodatavalue)

    # reclassify elevation 
    # (1) if 'bins' is a string
    if isinstance(bins, str):
        bin_name = bins

        if bin_name == 'median':
            # loop through bases
            for base in unique_bases:
                smask = bases == base            
                #Bin the elevation data based on the dem median
                bins = [np.min(dem[smask]), np.median(dem[smask]), np.max(dem[smask])]
                (hist,bin_edges) = np.histogram(dem[smask],bins=bins)            
                #Place the mean elev and elev band
                for ibin in np.arange(len(bin_edges)-1):
                    smask = (bases == base) & (dem >= bin_edges[ibin]) & (dem < bin_edges[ibin+1])
                    elev_class[smask] = ibin+1  # ele_band starts from one, not zero.
                    elev_value[smask] = np.mean(dem[smask])                

    # (2) if 'bins' is an integer or a sequence
    elif (np.ndim(bins) == 0) or (np.ndim(bins) == 1): 
        # loop through bases
        for base in unique_bases:
            smask = bases == base
            #Bin the elevation data based on bins
            (hist,bin_edges) = np.histogram(dem[smask],bins=bins)
            #Place the mean elev and elev band
            for ibin in np.arange(len(bin_edges)-1):
                smask = (bases == base) & (dem >= bin_edges[ibin]) & (dem < bin_edges[ibin+1])
                elev_class[smask] = ibin+1  # ele_band starts from one, not zero.
                elev_value[smask] = np.mean(dem[smask])                

    # mask elev_class and elev_value
    elev_class_ma = np.ma.masked_array(elev_class,dem_mask==0)
    elev_value_ma = np.ma.masked_array(elev_value,dem_mask==0) 

    # save elev_class_ma and elev_value_ma into rasters
    # update data type of meta and save raster
    out_meta.update(count=1, dtype='int32', compress='lzw', nodata=nodatavalue)
    with rio.open(dem_class_raster, 'w', **out_meta) as outf:
        outf.write(elev_class_ma, 1)    
    # update data type of meta and save raster
    out_meta.update(count=1, dtype='float32', compress='lzw', nodata=nodatavalue)
    with rio.open(dem_value_raster, 'w', **out_meta) as outf:
        outf.write(elev_value_ma, 1)
    return

def classify_landcover(lc_raster, lc_class_raster, nodatavalue):
    '''
    lc_raster: input, raster, landcover raster file.
    lc_class_raster: output, raster, elevation class.
    nodatavalue'''
    
    # read landcover data
    with rio.open(lc_raster) as ff:
        lc  = ff.read(1)
        lc_mask = ff.read_masks(1)
        out_meta = ff.meta.copy()

    # define array and specify dtype!
    lc_class = np.ones(np.shape(lc), dtype=np.int32)*int(nodatavalue)

    # specify landcover class
    # lc_class 1: crop. 2: non-crop. 
    ## Crop class includes: 0 Evergreen needleleaf forests, 1 Evergreen broadleaf forests, 
    ## 2 Deciduous needleleaf forests, 3 Deciduous broadleaf forests, 4 Mixed forests, 7 Woody savannas.  
    lc_class[(lc<5) | (lc==7)]=1  # crop
    lc_class[(lc>=5) & (lc!=7) & (lc!=255)]=2 # non-crop    

    # mask lc_class
    lc_class_ma = np.ma.masked_array(lc_class,lc_mask==0)

    # save lc_class_ma into rasters
    out_meta.update(count=1, dtype='int32', compress='lzw', nodata=nodatavalue)
    with rio.open(lc_class_raster, 'w', **out_meta) as outf:
        outf.write(lc_class_ma, 1)    
    return

def polygonize_raster(inraster, nodatavalue, outvector, attrb_field, attrb_field_typ):
    '''
    inraster: input, raster, source to convert to vector.
    nodatavalue: input, int, no data value of the inraster.
    outvector: output, vector, output vector.
    attrb_field: input, str, attribute field name of vector attribtue table to save raster values.
    attrb_field_typ: input, str, type of the output attribute field. For example, int'''
    # reference: https://programtalk.com/vs2/?source=python/7910/rasterio/examples/rasterio_polygonize.py
    
    driver='ESRI Shapefile'     
    with rio.Env():
        
        # read raster
        with rio.open(inraster) as src:
            image = src.read(1)
            mask = src.read_masks(1)
            
        # extract shapes of raster features
        results = (
            {'properties': {attrb_field: v}, 'geometry': s}
            for i, (s, v) in enumerate(rio.features.shapes(image, mask=mask, transform=src.transform)))
        
        # write GeoJSON style geometry dict to shapefile
        with fiona.open(
                outvector, 'w',
                driver=driver,
                crs=fiona.crs.to_string(src.crs),
                schema={'properties': [(attrb_field, attrb_field_typ)],
                        'geometry': 'Polygon'}) as dst:
            dst.writerecords(results)      
    return

def define_hru(raster_list, raster_fieldname_list, basin_raster, sub_corr_txt, subNo_field, subName_field,
               nodatavalue, outraster, outvector, hruNo_field, hruNo_field_type, hruName_field):
    '''
    raster_list: input, list, a list raster inputs that are used to define HRU.
    raster_fieldname_list: input, str, a list of field names corresponding to raster_list.
    basin_raster: input, raster, subbasin input raster.
    sub_corr_txt: input, text, subNo-HUC12 correspondence txt file.
    nodatavalue: no data value for HRU polygonization.
    outraster: output, raster, HRU raster with integer value.
    outvector: output, raster, HRU vector with HRU_full and HRU_int fields. '''
    
    # remove output files if exist
    for file in [outraster,outvector]:
        if os.path.exists(file): 
            os.remove(file)
        
    # --- PART 1. rasterize HRU via raster overlay (i.e., array concatenation) ---
    raster_str_max_len_list = [] # recrod the max str lenght of input class str
    for iraster, raster in enumerate(raster_list):
        # read raster
        with rio.open(raster) as ff:
            raster_value  = ff.read(1)
            raster_mask = ff.read_masks(1)

        # convert raster value to str, and format str to the same length (i.e., the max str length) 
        raster_str = raster_value.astype(str)
        max_len = max(map(len,raster_str[raster_mask!=0]))
        raster_str_max_len_list.append(max_len)
        raster_str_fmt = np.char.zfill(raster_str,max_len)

        # concatenate element-wise two arrays to generate HRU
        if iraster == 0:
            hru_str_fmt = raster_str_fmt
            hru_mask = (raster_mask!=0)
        else:
            hru_str_fmt = np.char.add(hru_str_fmt, raster_str_fmt)
            hru_mask = (hru_mask & raster_mask)

    # assign the unique integer to unique HRU (HRU_str_fmt)
    hru_int = np.zeros(np.shape(hru_str_fmt), dtype=np.int32)
    unique_hrus_str = np.unique(hru_str_fmt[hru_mask!=0]) 
    for ihru, hru in enumerate(unique_hrus_str):
        hru_mask = hru_str_fmt == hru
        hru_int[hru_mask] = int(ihru)+1
    hru_int_ma = np.ma.masked_array(hru_int,hru_mask==0)

    # save hru_int_ma into raster based on basin_raster 
    with rio.open(basin_raster) as src:
        dst_meta = src.meta.copy()
    with rio.open(outraster, 'w', **dst_meta) as dst:
        dst.write(hru_int_ma, 1)

    # --- PART 2. polgonize HRU and add attributes ---
    # polygonize hru rater
    polygonize_raster(outraster, nodatavalue, outvector, hruNo_field, hruNo_field_type)

    # dissolve and exclude HRU=0 blank area
    hru_gpd = gpd.read_file(outvector)
    hru_gpd_disv = hru_gpd.dissolve(by=hruNo_field)
    hru_gpd_disv = hru_gpd_disv.reset_index()
    hru_gpd_disv = hru_gpd_disv.loc[hru_gpd_disv.loc[:,hruNo_field] > 0]
    hru_gpd_disv = hru_gpd_disv.reset_index()

    # add other fields: HRU and detailed input raster classes
    # Note: since we want to name HRU as HUC12+hru# (xxxx01, xxxx02), 
    # we need to find out the string location of subbasin among the compelte hru str.
    basin_index = raster_list.index(basin_raster)
    basin_str_start = sum(raster_str_max_len_list[0:basin_index])
    basin_str_len = raster_str_max_len_list[basin_index]

    # read subNo-HUC12 corresponding text file to get subNo <-> HUC12.
    corr_df = pd.read_csv(sub_corr_txt, sep=",")

    # loop through HRUs to add attributes
    for irow, row in hru_gpd_disv.iterrows():

        # (a) add new field columns, define other needed variables.
        if irow == 0:
            hru_gpd_disv[hruName_field] = ""      # simplified HRU name, e.g., hruNo + HRU#.
            hru_gpd_disv[subName_field] = ""    # HUC12 int
            for jfield, field in enumerate(raster_fieldname_list): # other fields
                hru_gpd_disv[field] = ""          
            hru_num_per_sub = 0           # count the number of HRUs per subbasin
            subNo_before = ''             # subNo before this irow to update hru_num_per_sub
        else:
            subNo_before = subNo_current

        # (b) identify current subbasin ID (subNo_current: 1,2,3...)
        hru_str = unique_hrus_str[row[hruNo_field]-1]
        subNo_current = hru_str[basin_str_start:basin_str_start+basin_str_len] 

        # (c) update the numebr of HRUs whtin subNo_current (hru_num_per_sub: 1,2)
        if subNo_current == subNo_before:  # if the current subNo is the same as before, cumulate hru num.
            hru_num_per_sub = hru_num_per_sub+1
        else: # otherwise, reset hru num as zero.
            hru_num_per_sub = 0

        # (d) fill HRU field = HUC12 + hru_count in subNo_current, and HUC12 field.   
        sub_HUC = corr_df[corr_df[subNo_field]==int(subNo_current)].reset_index().loc[0,subName_field]
        hru_gpd_disv.loc[irow,hruName_field] = str(sub_HUC) + str(hru_num_per_sub+1).zfill(2) # two-digit HRU#
        hru_gpd_disv.loc[irow,subName_field] = str(sub_HUC) 

        # (e) fill other fields for class checks
        for jfield, field in enumerate(raster_fieldname_list):
            field_str_start = sum(raster_str_max_len_list[0:jfield])
            field_str_len = raster_str_max_len_list[jfield]
            hru_gpd_disv.loc[irow,field] = hru_str[field_str_start:field_str_start+field_str_len]          

    # re-order two columns (HUC12, subNo -> subNo, HUC12)
    cols = list(hru_gpd_disv.columns)
    a, b = cols.index(subName_field), cols.index(subNo_field)
    cols[b], cols[a] = cols[a], cols[b]
    hru_gpd_disv = hru_gpd_disv[cols]
    
    # save gpd to shapefile
    hru_gpd_disv.to_file(outvector, index=False)
    return

def zonal_statistic(attr_raster, invector, infield, infield_dtype, refraster, metric, out_raster, nodatavalue, *args, **kwargs):
    '''
    attr_raster: input, raster. Raster to analyze.
    invector: input, vector. Vector with zones boundaries.
    infield: input, string. Field name of invector that is used to identify zone of invector.
    infield_dtype: input, string. data type of infield, used to rasterize invector.
    metric: input, string. metric to calcualte zonal statistics. Choose from 'mean', 'median', 'max', 'min'.
    nodatavalue: no data value.   
    out_raster: output, raster. Raster to save the zonal attribtue value of invector. 
    raster_band: optional input, integer. Number of raster band to analyze.
    output_column_prefix: optional input, str. Prefix for output fields.'''

    output_column_prefix = kwargs.get('output_column_prefix', 'DN')
    raster_band = kwargs.get('raster_band', 1)
        
    # --- PART 1. rasterize invector based on refraster ---
    in_raster = os.path.join(os.path.dirname(invector), os.path.basename(invector).split('.')[0]+'tmp.shp')
    rasterize_vector(invector,infield,infield_dtype,refraster,in_raster,nodatavalue)

    # --- PART 2. loop through HRU to get attribute value ---
    # read invector, in_raster and attr_raster
    in_gpd = gpd.read_file(invector)
    with rio.open(in_raster) as ff:
        in_raster_value  = ff.read(1)
        in_raster_mask = ff.read_masks(1)
    with rio.open(attr_raster) as ff:
        attr_raster_value  = ff.read(raster_band)
        out_meta = ff.meta.copy()

    # initialize the output column prefix 
    if output_column_prefix in in_gpd.columns:
        in_gpd = in_gpd.drop(columns=[output_column_prefix])
    in_gpd[output_column_prefix] = "" 

    # initialize the output array
    output_array = (np.ones_like(attr_raster_value))*nodatavalue
    
    # loop through invector polygons
    unique_field_value = np.unique(in_raster_value[in_raster_mask!=0])
    for field_value in unique_field_value:
        smask = in_raster_value == field_value         
        if isinstance(metric, str):
            if metric == 'mean':
                zonal_value = np.mean(attr_raster_value[smask])
            elif metric == 'mode':
                zonal_value = stats.mode(attr_raster_value[smask])[0]
            in_gpd.loc[in_gpd[infield]==field_value, output_column_prefix] = zonal_value
            output_array[smask] = zonal_value
        else: 
            print("Invalid metric string. Please choose from 'mean', 'median', 'max', 'min'.")
    
    # save vector
    in_gpd.to_file(invector)
    os.remove(in_raster)        

    # save outraster
    with rio.open(out_raster, 'w', **out_meta) as outf:
        out_arr_ma = np.ma.masked_array(output_array, in_raster_mask==0)
        outf.write(out_arr_ma,1)             
    return

def plot_vector(invector, column):    
    in_gpd = gpd.read_file(invector)
    fig, ax = plt.subplots(figsize=(10, 6))
    in_gpd.plot(column=column, ax=ax)
    ax.set_axis_off()
    # plt.axis('equal')
    plt.show()   
    return
