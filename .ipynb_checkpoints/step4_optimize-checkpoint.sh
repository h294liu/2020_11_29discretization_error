#!/bin/bash
set -e

# H. liu, April 28, 2020.
# Convert original NLDAS netcdf to the GMET output format.
# This is for ensemble bias correction.
  
source_dir=/glade/u/home/hongli/github/2020_11_29discretization_error
ostIn_tpl=${source_dir}/calib_tpl/ostIn.txt.tpl
discretize_code_tpl=${source_dir}/calib_tpl/step2_generate_HRU.py.tpl
discretize_functions=${source_dir}/geospatial_functions

work_dir=/glade/u/home/hongli/scratch/2020_11_29discretization_error
case=yampa
case_discretize_dir=${work_dir}/${case}
case_optimize_dir=${work_dir}/${case}_optimize
Ostrich_exe=/glade/u/home/hongli/tools/ostrich/OstrichGCC
MaxIterations=3

if [ ! -d $case_discretize_dir ]; then mkdir -p $case_discretize_dir; fi
if [ ! -d $case_optimize_dir ]; then mkdir -p $case_optimize_dir; fi

# Create optimize dictionary and files
cd $case_optimize_dir

rm -f Ostrich.exe
ln -s $Ostrich_exe Ostrich.exe

ostIn_basename=${ostIn_tpl##*/}         # get basename of filename
ostIn="${ostIn_basename/.tpl/}" # remove suffix ".tpl"
if [ -e ${ostIn} ]; then rm -rf ${ostIn}; fi
echo $ostIn

sed "s,ROOT_DIR,$work_dir,g" $ostIn_tpl |\
sed "s,CASE,$case,g" |\
sed "s,MAXITERATIONS,$MaxIterations,g" > $ostIn
chmod 740 ${ostIn}

# Create model dictionary and files
if [ -e ${discretize_functions} ]; then rm -rf ${discretize_functions}; fi
mkdir model
cd model
cp -r $discretize_functions .

discretize_code="${discretize_code_tpl/.tpl/}"
discretize_code_basename=${ostIn_tpl##*/}         # get basename of filename
discretize_code="${discretize_code_basename/.tpl/}" # remove suffix ".tpl"
if [ -e ${discretize_code} ]; then rm -rf ${discretize_code}; fi
echo $discretize_code

sed "s,ROOT_DIR,$work_dir,g" $discretize_code_tpl |\
sed "s,CASE,$case,g" > $discretize_code
chmod 740 ${discretize_code}

cd ../
# ./Ostrich.exe ostIn.txt