#!/bin/bash
set -e

# H. liu, April 28, 2020.
# Convert original NLDAS netcdf to the GMET output format.
# This is for ensemble bias correction.
  
source_code_dir=/glade/u/home/hongli/github/2020_11_29discretization_error
ostIn_tpl=$source_code_dir/calib_tpl/ostIn.tpl.txt
run_model_tpl=$source_code_dir/calib_tpl/run_model.sh
discretize_code_tpl=$source_code_dir/calib_tpl/step2_generate_HRU.tpl.py
discretize_functions=$source_code_dir/geospatial_functions

case=yampa
source_data_dir=/glade/u/home/hongli/scratch/2020_11_29discretization_error/$case
dem_crop=$source_data_dir/dem_crop.tif
dem_class_raster=$source_data_dir/dem_class.tif
lc_class_raster=$source_data_dir/landcover_class.tif
slp_crop=$source_data_dir/slope_crop.tif
asp_crop_180=$source_data_dir/aspect_crop_180.tif
sw_raster=$source_data_dir/step9_merge_raw_Sw/sw.tif
sub_raster=$source_data_dir/subbasin.tif
sub_corr_txt=$source_data_dir/subNo_HUC12_corr.txt

optimize_dir=/glade/u/home/hongli/scratch/2020_11_29discretization_error/${case}_optimize
if [ ! -d $optimize_dir ]; then mkdir -p $optimize_dir; fi
model_dir=${optimize_dir}/model
Ostrich_exe=/glade/u/home/hongli/tools/ostrich/OstrichGCC
MaxIterations=30

# --- Part 1. Optimize --- 
cd $optimize_dir
# (1) Ostrich.exe
rm -f Ostrich.exe
ln -s $Ostrich_exe Ostrich.exe
# (2) ostIn.txt
ostIn_basename=${ostIn_tpl##*/}         # get basename of filename
ostIn="${ostIn_basename/.tpl.txt/.txt}" # remove ".tpl"
if [ -e ${ostIn} ]; then rm -rf ${ostIn}; fi
sed "s,ROOT_DIR,$model_dir,g" $ostIn_tpl |\
sed "s,CASE,$case,g" |\
sed "s,MAXITERATIONS,$MaxIterations,g" > $ostIn
chmod 740 ${ostIn}
# (3) run_model.sh
run_model_file=${run_model_tpl##*/}
if [ -e ${run_model_file} ]; then rm -rf ${run_model_file}; fi
cp $run_model_tpl $run_model_file
# (4) discretization code
discretize_code=${discretize_code_tpl##*/}       # get basename of filename
if [ -e ${discretize_code} ]; then rm -rf ${discretize_code}; fi
sed "s,ROOT_DIR,$model_dir,g" $discretize_code_tpl > $discretize_code
chmod 740 ${discretize_code}

# --- Part 2. Model --- 
if [ -e model ]; then rm -rf model; fi
mkdir model
cd model
# (1) discretization function library
cp -r $discretize_functions .
# (2) copy needed raster and vector files
cp $dem_crop .
cp $dem_class_raster .
cp $lc_class_raster .
cp $slp_crop .
cp $asp_crop_180 .
cp $sw_raster .
cp $sub_raster .
cp $sub_corr_txt .
cd ../

# --- Part 3. Submit job --- 
./Ostrich.exe ostIn.txt