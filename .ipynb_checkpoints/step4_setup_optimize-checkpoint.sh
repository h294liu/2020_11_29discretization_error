#!/bin/bash
set -e

# H. liu, April 28, 2020.
# Convert original NLDAS netcdf to the GMET output format.
# This is for ensemble bias correction.
  
source_code_dir=/glade/u/home/hongli/github/2020_11_29discretization_error
ostIn_tpl_path=$source_code_dir/calib_tpl/ostIn.tpl.txt
run_model_tpl_path=$source_code_dir/calib_tpl/run_model.tpl.sh
dis_code_tpl_path=$source_code_dir/calib_tpl/step2_generate_HRU.tpl.py
# dis_code_tpl_path=$source_code_dir/calib_tpl/step2_generate_HRU_eliminate.tpl.py
dis_functions=$source_code_dir/geospatial_functions

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

Ostrich_exe=/glade/u/home/hongli/tools/ostrich/OstrichGCC
MaxIterations=5000
opt_num=3

for i in $(seq 1 $opt_num); do

    echo $i

    optimize_dir=/glade/u/home/hongli/scratch/2020_11_29discretization_error/${case}_optimize${i}
    # optimize_dir=/glade/u/home/hongli/scratch/2020_11_29discretization_error/${case}_optimize_eliminate
    if [ -d $optimize_dir ]; then rm -rf $optimize_dir; fi
    mkdir -p $optimize_dir
    model_dir=${optimize_dir}/model

    # --- Part 1. Optimize --- 
    cd $optimize_dir
    # (1) discretization code
    dis_code_tpl=${dis_code_tpl_path##*/}  # get basename of filename
    dis_code="${dis_code_tpl/.tpl.py/.py}" # remove ".tpl"
#     if [ -e ${dis_code_tpl} ]; then rm -rf ${dis_code_tpl}; fi
    sed "s,ROOT_DIR,$model_dir,g" $dis_code_tpl_path > $dis_code_tpl
    chmod 740 ${dis_code_tpl}

    # (2) run_model.sh
    run_model_tpl=${run_model_tpl_path##*/}
    run_model="${run_model_tpl/.tpl.sh/.sh}" # remove ".tpl"
#     if [ -e ${run_model} ]; then rm -rf ${run_model}; fi
    sed "s,DISCRETIZE_CODE,$dis_code,g" $run_model_tpl_path > $run_model
    chmod 740 ${run_model}

    # (3) Ostrich.exe
    rm -f Ostrich.exe
    ln -s $Ostrich_exe Ostrich.exe

    # (4) ostIn.txt
    ostIn_tpl=${ostIn_tpl_path##*/}     # get basename of filename
    ostIn="${ostIn_tpl/.tpl.txt/.txt}"  # remove ".tpl"
#     if [ -e ${ostIn} ]; then rm -rf ${ostIn}; fi
    sed "s,ROOT_DIR,$model_dir,g" $ostIn_tpl_path |\
    sed "s,DISCRETIZE_CODE_TPL,$dis_code_tpl,g" |\
    sed "s,DISCRETIZE_CODE,$dis_code,g" |\
    sed "s,MAXITERATIONS,$MaxIterations,g" > $ostIn
    chmod 740 ${ostIn}

    # --- Part 2. Model --- 
    mkdir model
    cd model
    # (1) discretization function library
    cp -r $dis_functions .

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
    # ./Ostrich.exe ostIn.txt
    CommandFile=${optimize_dir}/qsub.sh       
    LogFile=${optimize_dir}/log
#     rm -f $CommandFile
#     rm -f $LogFile.*
    jobName=${optimize_dir##*/}

    echo '#!/bin/bash' > ${CommandFile}
    echo "#PBS -N ${jobName}" >> ${CommandFile}
    echo '#PBS -A P48500028' >> ${CommandFile}
    echo '#PBS -q regular' >> ${CommandFile}
    echo '#PBS -l select=1:ncpus=1:mpiprocs=1' >> ${CommandFile}
    echo '#PBS -l walltime=12:00:00' >> ${CommandFile}
    echo "#PBS -o ${LogFile}.o" >> ${CommandFile}
    echo "#PBS -e ${LogFile}.e" >> ${CommandFile}

    echo "mkdir -p /glade/scratch/hongli/temp" >> ${CommandFile}
    echo "export TMPDIR=/glade/scratch/hongli/temp" >> ${CommandFile}

    echo "cd ${optimize_dir}" >> ${CommandFile}
    echo "./Ostrich.exe ostIn.txt" >> ${CommandFile}
    # echo "${Program} ${ConfigFile}" >> ${CommandFile}
    chmod 740 ${CommandFile}

    qsub ${CommandFile}
done