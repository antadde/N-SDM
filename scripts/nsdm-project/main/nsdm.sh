#!/bin/bash

##########################
## nsdm.sh
## Core N-SDM script
## Date: 20-05-2022
## Author: Antoine Adde 
##########################

# Load required modules [cluster specific]
module load gcc/9.3.0
module load r/4.0.5
module load proj/5.2.0
module load perl
module load curl
module load geos/3.8.1
module load gdal/2.4.4-proj-5.2.0

# Retrieve main paths
wp=$(awk -F ";" '$1 == "w_path" { print $2}' ./settings/settings.csv)    # working path
sop=$(awk -F ";" '$1 == "scr_path" { print $2}' ./settings/settings.csv) # scratch output path
svp=$(awk -F ";" '$1 == "svp_path" { print $2}' ./settings/settings.csv) # saving output path

# Retrieve project name
project=$(awk -F ";" '$1 == "project" { print $2}' ./settings/settings.csv)

# General definitions
own=$(awk -F ";" '$1 == "sess_own" { print $2}' ./settings/settings.csv)                          # session account
acc=$(awk -F ";" '$1 == "account" { print $2}' ./settings/settings.csv)                           # HPC account
part=$(awk -F ";" '$1 == "partition" { print $2}' ./settings/settings.csv)                        # HPC partition

# Clean and/or create output directories
rm -r $wp/outputs/$project/* 2>/dev/null
mkdir -p $svp/outputs/$project/ 2>/dev/null
mkdir -p $sop/outputs/$project/ 2>/dev/null
mkdir -p $sop/tmp/$project/ 2>/dev/null
mkdir -p $wp/outputs/$project/settings/tmp/ 2>/dev/null
mkdir -p $wp/outputs/$project/sacct/ 2>/dev/null

# Permissions
chmod -R 777 $wp/data/$project 2>/dev/null
chmod -R 777 $wp/scripts/$project 2>/dev/null
chmod -R 777 $wp/outputs/$project 2>/dev/null

# Clean existing 0_mainPRE log files if any
rm $wp/scripts/$project/main/0_mainPRE/*.err 2>/dev/null
rm $wp/scripts/$project/main/0_mainPRE/*.out 2>/dev/null

# Create N-SDM unique identifier
ssl_id=$(openssl rand -hex 3)
echo $ssl_id > $wp/outputs/$project/settings/tmp/ssl_id.txt

echo Welcome to this new $project run $ssl_id

# Define N-SDM settings and disaggregate species occurence data
pre_A_m=$(awk -F ";" '$1 == "pre_A_m" { print $2}' ./settings/settings.csv)                           # memory
pre_A_t=$(awk -F ";" '$1 == "pre_A_t" { print $2}' ./settings/settings.csv)                           # time
pre_A_c=$(awk -F ";" '$1 == "pre_A_c" { print $2}' ./settings/settings.csv)                           # cores
sbatch --wait --account=$acc --partition=$part --mem=$pre_A_m --time=$pre_A_t --cpus-per-task=$pre_A_c --ntasks=1 ./0_mainPRE/job_pre_A.sh
echo N-SDM settings defined and species occurence data for "$(cat $wp/outputs/$project/settings/tmp/n_spe.txt)" species disaggregated 

### Loop over spe_runs (n_spe/n_spe_max_hpc) to prevent scratch path saturation
spe_runs="$(cat $wp/outputs/$project/settings/tmp/spe_runs.txt)"

for i in $(seq 1 $spe_runs)
do
echo $i > $wp/outputs/$project/settings/tmp/run_id.txt
echo "Starting N-SDM run $i out of $spe_runs runs"

# Retrieve time/date
dt=$(date +"%FT%T")

# Update N-SDM settings
cd $wp/scripts/$project/main/
Rscript ./0_mainPRE/nsdm_update.R 1>/dev/null 2>&1
echo N-SDM settings updated

# Update n_spe (number of species to be modelled in this run)
n_spe="$(cat $wp/outputs/$project/settings/tmp/n_spe.txt)"

# Dimensions for array definitions
mod_algo=$(awk 'BEGIN{FS=";"}/mod_algo/ {print $2}' ./settings/settings.csv)        # modelling algorithms evaluated
n_algo=$(($(grep -o "'" <<<"$mod_algo" | grep -c .) / 2))                           # number of modelling algorithms evaluated
nest_met=$(awk 'BEGIN{FS=";"}/nesting_methods/ {print $2}' ./settings/settings.csv) # modelling algorithms evaluated
n_nesting=$(($(grep -o "'" <<<"$nest_met" | grep -c .) / 2))                        # number of nesting methods evaluated
scenars=$(awk 'BEGIN{FS=";"}/proj_scenarios/ {print $2}' ./settings/settings.csv)   # alternative scenarios evaluated
n_scenarios=$(($(grep -o "'" <<<"$scenars" | grep -c .) / 2))                       # number of alternative scenarios evaluated
periods=$(awk 'BEGIN{FS=";"}/proj_periods/ {print $2}' ./settings/settings.csv)          # alternative periods evaluated
n_periods=$(($(grep -o "'" <<<"$periods" | grep -c .) / 2))                         # number of alternative periods evaluated

# Memory/Time/Cores/Array definitions
## PRE_B
pre_B_m=$(awk -F ";" '$1 == "pre_B_m" { print $2}' ./settings/settings.csv)                           # memory
pre_B_t=$(awk -F ";" '$1 == "pre_B_t" { print $2}' ./settings/settings.csv)                           # time
pre_B_c=$(awk -F ";" '$1 == "pre_B_c" { print $2}' ./settings/settings.csv)                           # cores

## GLO
glo_A_m=$(awk -F ";" '$1 == "glo_A_m" { print $2}' ./settings/settings.csv)                           # A memory
glo_A_t=$(awk -F ";" '$1 == "glo_A_t" { print $2}' ./settings/settings.csv)                           # A time
glo_A_c=$(awk -F ";" '$1 == "glo_A_c" { print $2}' ./settings/settings.csv)                           # A cores
glo_A_a=$n_spe                                                                                        # A array extent

glo_B_m=$(awk -F ";" '$1 == "glo_B_m" { print $2}' ./settings/settings.csv)                           # B memory
glo_B_t=$(awk -F ";" '$1 == "glo_B_t" { print $2}' ./settings/settings.csv)                           # B time
glo_B_c=$(awk -F ";" '$1 == "glo_B_c" { print $2}' ./settings/settings.csv)                           # B cores
glo_B_a=`expr $n_spe \* $n_algo`                                                                      # B array extent

glo_C_m=$(awk -F ";" '$1 == "glo_C_m" { print $2}' ./settings/settings.csv)                           # C memory
glo_C_t=$(awk -F ";" '$1 == "glo_C_t" { print $2}' ./settings/settings.csv)                           # C time
glo_C_c=$(awk -F ";" '$1 == "glo_C_c" { print $2}' ./settings/settings.csv)                           # C cores
glo_C_a=$n_spe                                                                                        # C array extent

## LOC
loc_A_m=$(awk -F ";" '$1 == "loc_A_m" { print $2}' ./settings/settings.csv)                           # A memory
loc_A_t=$(awk -F ";" '$1 == "loc_A_t" { print $2}' ./settings/settings.csv)                           # A time
loc_A_c=$(awk -F ";" '$1 == "loc_A_c" { print $2}' ./settings/settings.csv)                           # A cores
loc_A_a=$n_spe                                                                                        # A array extent

loc_B_m=$(awk -F ";" '$1 == "loc_B_m" { print $2}' ./settings/settings.csv)                           # B memory
loc_B_t=$(awk -F ";" '$1 == "loc_B_t" { print $2}' ./settings/settings.csv)                           # B time
loc_B_c=$(awk -F ";" '$1 == "loc_B_c" { print $2}' ./settings/settings.csv)                           # B cores
loc_B_a=`expr $n_spe \* $n_algo \* $n_nesting`                                                        # B array extent

loc_C_m=$(awk -F ";" '$1 == "loc_C_m" { print $2}' ./settings/settings.csv)                           # C memory
loc_C_t=$(awk -F ";" '$1 == "loc_C_t" { print $2}' ./settings/settings.csv)                           # C time
loc_C_c=$(awk -F ";" '$1 == "loc_C_c" { print $2}' ./settings/settings.csv)                           # C cores
loc_C_a=`expr $n_spe \* $n_nesting`                                                                   # C array extent   

## FUT
fut_A_m=$(awk -F ";" '$1 == "fut_A_m" { print $2}' ./settings/settings.csv)                           # A memory
fut_A_t=$(awk -F ";" '$1 == "fut_A_t" { print $2}' ./settings/settings.csv)                           # A time
fut_A_c=$(awk -F ";" '$1 == "fut_A_c" { print $2}' ./settings/settings.csv)                           # A cores
fut_A_a=`expr $n_spe \* $n_algo \* $n_scenarios`                                                      # A array extent

fut_B_m=$(awk -F ";" '$1 == "fut_B_m" { print $2}' ./settings/settings.csv)                           # B memory
fut_B_t=$(awk -F ";" '$1 == "fut_B_t" { print $2}' ./settings/settings.csv)                           # B time
fut_B_c=$(awk -F ";" '$1 == "fut_B_c" { print $2}' ./settings/settings.csv)                           # B cores
fut_B_a=`expr $n_spe \* $n_scenarios`                                                                 # B array extent

fut_C_m=$(awk -F ";" '$1 == "fut_C_m" { print $2}' ./settings/settings.csv)                           # C memory
fut_C_t=$(awk -F ";" '$1 == "fut_C_t" { print $2}' ./settings/settings.csv)                           # C time
fut_C_c=$(awk -F ";" '$1 == "fut_C_c" { print $2}' ./settings/settings.csv)                           # C cores
fut_C_a=`expr $n_spe \* $n_algo \* $n_nesting \* $n_scenarios`                                        # C array extent

fut_D_m=$(awk -F ";" '$1 == "fut_D_m" { print $2}' ./settings/settings.csv)                           # D memory
fut_D_t=$(awk -F ";" '$1 == "fut_D_t" { print $2}' ./settings/settings.csv)                           # D time
fut_D_c=$(awk -F ";" '$1 == "fut_D_c" { print $2}' ./settings/settings.csv)                           # D cores
fut_D_a=`expr $n_spe \* $n_nesting \* $n_scenarios`                                                   # D array extent

## END_A
end_A_m=$(awk -F ";" '$1 == "end_A_m" { print $2}' ./settings/settings.csv)                           # D memory
end_A_t=$(awk -F ";" '$1 == "end_A_t" { print $2}' ./settings/settings.csv)                           # D time
end_A_c=$(awk -F ";" '$1 == "end_A_c" { print $2}' ./settings/settings.csv)                           # D cores
end_A_a=$n_spe                                                                                        # D array extent

# Clean existing log files
rm $wp/scripts/$project/main/0_mainPRE/pre_B*.err 2>/dev/null
rm $wp/scripts/$project/main/0_mainPRE/pre_B*.out 2>/dev/null
rm $wp/scripts/$project/main/1_mainGLO/*.err 2>/dev/null
rm $wp/scripts/$project/main/1_mainGLO/*.out 2>/dev/null
rm $wp/scripts/$project/main/2_mainLOC/*.err 2>/dev/null
rm $wp/scripts/$project/main/2_mainLOC/*.out 2>/dev/null
rm $wp/scripts/$project/main/3_mainFUT/*.err 2>/dev/null
rm $wp/scripts/$project/main/3_mainFUT/*.out 2>/dev/null
rm $wp/scripts/$project/main/4_mainEND/*.err 2>/dev/null
rm $wp/scripts/$project/main/4_mainEND/*.out 2>/dev/null

# Clean scratch output folder if requested
clear_sop=$(awk -F ";" '$1 == "clear_sop" { print $2}' ./settings/settings.csv)
if [ $clear_sop = "TRUE" ]
then
rm -r $sop/outputs/$project/* 2>/dev/null
rm -r $sop/tmp/$project/* 2>/dev/null
fi

# Start running jobs
## n_levels of analyses (1=GLO; 2=GLO+LOC)?
n_levels=$(awk -F ";" '$1 == "n_levels" { print $2}' ./settings/settings.csv) # Skip future analyses ?

## Do future analyses?
do_proj=$(awk -F ";" '$1 == "do_proj" { print $2}' ./settings/settings.csv) # Skip future analyses ?

## PRE_B
cd $wp/scripts/$project/main/0_mainPRE
sbatch --wait --account=$acc --partition=$part --mem=$pre_B_m --time=$pre_B_t --cpus-per-task=$pre_B_c --ntasks=1 job_pre_B.sh
echo PRE modelling datasets generated

## GLO level
cd $wp/scripts/$project/main/1_mainGLO
sbatch --wait --account=$acc --partition=$part --mem=$glo_A_m --time=$glo_A_t --cpus-per-task=$glo_A_c --ntasks=1 --array [1-$glo_A_a] job_glo_A.sh
echo GLO data preparation and covariate selection done
sbatch --wait --account=$acc --partition=$part --mem=$glo_B_m --time=$glo_B_t --cpus-per-task=$glo_B_c --ntasks=1 --array [1-$glo_B_a] job_glo_B.sh
echo GLO modelling done
sbatch --wait --account=$acc --partition=$part --mem=$glo_C_m --time=$glo_C_t --cpus-per-task=$glo_C_c --ntasks=1 --array [1-$glo_C_a] job_glo_C.sh
echo GLO ensembling done

if [ $n_levels -gt 1 ]
then 
## LOC level
cd $wp/scripts/$project/main/2_mainLOC
sbatch --wait --account=$acc --partition=$part --mem=$loc_A_m --time=$loc_A_t --cpus-per-task=$loc_A_c --ntasks=1 --array [1-$loc_A_a] job_loc_A.sh
echo LOC data preparation and covariate selection done
sbatch --wait --account=$acc --partition=$part --mem=$loc_B_m --time=$loc_B_t --cpus-per-task=$loc_B_c --ntasks=1 --array [1-$loc_B_a] job_loc_B.sh
echo LOC modelling done
sbatch --wait --account=$acc --partition=$part --mem=$loc_C_m --time=$loc_C_t --cpus-per-task=$loc_C_c --ntasks=1 --array [1-$loc_C_a] job_loc_C.sh
echo LOC ensembling and scale nesting done
fi

## FUT predictions
if [ $do_proj = "TRUE" ]
then
cd $wp/scripts/$project/main/3_mainFUT
sbatch --wait --account=$acc --partition=$part --mem=$fut_A_m --time=$fut_A_t --cpus-per-task=$fut_A_c --ntasks=1 --array [1-$fut_A_a] job_fut_A.sh
echo individual FUT GLO predictions done
sbatch --wait --account=$acc --partition=$part --mem=$fut_B_m --time=$fut_B_t --cpus-per-task=$fut_B_c --ntasks=1 --array [1-$fut_B_a] job_fut_B.sh
echo FUT GLO ensembling done
sbatch --wait --account=$acc --partition=$part --mem=$fut_C_m --time=$fut_C_t --cpus-per-task=$fut_C_c --ntasks=1 --array [1-$fut_C_a] job_fut_C.sh
echo individual FUT LOC predictions done
sbatch --wait --account=$acc --partition=$part --mem=$fut_D_m --time=$fut_D_t --cpus-per-task=$fut_D_c --ntasks=1 --array [1-$fut_D_a] job_fut_D.sh
echo FUT LOC ensembling and scale nesting done
fi

## END analyses
cd $wp/scripts/$project/main/4_mainEND
sbatch --wait --account=$acc --partition=$part --mem=$end_A_m --time=$end_A_t --cpus-per-task=$end_A_c --ntasks=1 --array [1-$end_A_a] job_end_A.sh
echo Final evaluation done
sacct --starttime $dt -u $own --format JobID,JobName,Elapsed,NCPUs,TotalCPU,CPUTime,ReqMem,MaxRSS,MaxDiskRead,MaxDiskWrite,State,ExitCode > $wp/outputs/$project/sacct/"${ssl_id}_${i}_sacct.txt"
Rscript end_B.R 1>/dev/null 2>&1
echo Sacct outputs analysis done

# Permissions
chmod -R 777 $wp/scripts/$project/main

# rsync to saving location before cleaning scratch folder
cd $sop/outputs/$project/
find d2_models/ -name '*glm.rds' -o -name '*gam.rds' -o -name '*rf.rds' -o -name '*max.rds' -o -name '*gbm.rds' -o -name '*esm.rds' > $wp/outputs/$project/settings/tmp/modfiles.txt
rsync -a --files-from=$wp/outputs/$project/settings/tmp/modfiles.txt . $svp/outputs/$project
rsync -a --exclude={'d1_covsels','d2_models','d6_preds','d13_preds-fut','d14_maps-fut','d17_nested-ensembles-fut','d18_nested-ensembles-cv-fut'} $sop/outputs/$project/ $svp/outputs/$project
echo Main outputs sync to saving location
done

# Download, (pre-fill) and save ODMAP protocol
mkdir $wp/outputs/$project/ODMAP 2>/dev/null
curl -o $wp/outputs/$project/ODMAP/ODMAP.xlsx https://damariszurell.github.io/files/Zurell_etal_ODMAP.v1.0_TableS1.xlsx 2>/dev/null
rsync -a $wp/outputs/$project/ODMAP $svp/outputs/$project
echo ODMAP protocol generated

echo Finished!