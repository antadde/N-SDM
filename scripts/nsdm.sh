#!/bin/bash

##############################################
##############################################
## NSDM - Core N-SDM script
## Author: Antoine Adde (antoine.adde@eawag.ch)
## Created: 20-05-2022
## Last Update: 20-09-2025
##############################################
##############################################

echo "##############################################"
echo "#                                            #"
echo "#    N   N    SSSSS    DDDDDD    M     M     #"
echo "#    NN  N   S         D     D   MM   MM     #"
echo "#    N N N    SSSSS    D     D   M M M M     #"
echo "#    N  NN         S   D     D   M  M  M     #"
echo "#    N   N    SSSSS    DDDDDD    M     M     #"
echo "#                                            #"
echo "#          N-SDM SESSION STARTING            #"
echo "##############################################"

##############################################
#        SESSION MANAGEMENT                  #
##############################################
# Load nsdm functions
source "./helpers/functions.sh"

trap cleanup EXIT

# Load the software stack if needed
module load "$(get_value "software_stack")"

# Load required modules
module load "$(get_value "module_r")"
module load "$(get_value "module_proj")"
module load "$(get_value "module_perl")"
module load "$(get_value "module_curl")"
module load "$(get_value "module_geos")"
module load "$(get_value "module_gdal")"
module load "$(get_value "module_sqlite")"
module load "$(get_value "module_udunits")"

# Retrieve main paths from settings.psv
wp=$(get_value "w_path")    # Working path
sop=$(get_value "scr_path")  # Scratch output path
svp=$(get_value "svp_path")  # Saving output path
rlibs=$(get_value "lib_path")  # Saving output path

# General definitions from settings.psv
own=$(get_value "sess_own")  # Session account
acc=$(get_value "account")   # HPC account
part=$(get_value "partition") # HPC partition

# ssl_id directory
ssl_dir="$wp/tmp/ssl_id"
ssl_id_file="$wp/tmp/ssl_id/ssl_id.txt"

# Checkpoint directory
checkpoint_dir="$wp/tmp/checkpoints"

# SACCT log
sacct_dir="$wp/tmp/sacct"
sacct_log="$wp/tmp/sacct/sacct_log.txt"

# Check if a new session needs to be started (if "new_sess" is TRUE)
new_sess=$(get_value "new_sess")
if [ "$new_sess" = "TRUE" ]; then
echo "New session requested."
rm -f "$ssl_id_file"  # Delete the existing ssl_id file if new session is requested
fi

# Check if the file exists and load the existing session ID
if [ -f "$ssl_id_file" ]; then
ssl_id=$(cat "$ssl_id_file")
echo "Loaded existing Session ID: $ssl_id"
else
# Clean and/or create necessary directories
# Remove old temporary files, if any
rm -r "$wp/tmp/" 2>/dev/null || true

# Create required directories if they don't already exist
mkdir -p "$svp/outputs/" 2>/dev/null || true
mkdir -p "$sop/outputs/" 2>/dev/null || true
mkdir -p "$sop/tmp/" 2>/dev/null || true
mkdir -p "$wp/tmp/settings/" 2>/dev/null || true
mkdir -p "$wp/tmp/sacct/" 2>/dev/null || true

# Generate a new session ID if none exists
mkdir -p "$ssl_dir" 2>/dev/null || true
ssl_id=$(openssl rand -hex 3)
echo "$ssl_id" > "$ssl_id_file"
echo "Generated new Session ID: $ssl_id"

# Checkpoint directory
rm -r "$checkpoint_dir" 2>/dev/null || true
mkdir -p "$checkpoint_dir" 2>/dev/null || true

# SACCT directory and log
rm -r "$sacct_dir" 2>/dev/null || true
mkdir -p "$sacct_dir" 2>/dev/null || true

# Set directory permissions
chmod -R 777 "$wp/data" 2>/dev/null || true
chmod -R 777 "$wp/scripts" 2>/dev/null || true
chmod -R 777 "$wp/tmp" 2>/dev/null || true

# Clean up old log files from the mainPRE script if they exist
rm "$wp/scripts/0_mainPRE/logs/"*.err 2>/dev/null || true
rm "$wp/scripts/0_mainPRE/logs/"*.out 2>/dev/null || true
fi

##############################################
#        RUN PRE_A JOB                        
##############################################
PRE_A_m=$(get_value "pre_A_m")  # Memory
PRE_A_t=$(get_value "pre_A_t")  # Time
PRE_A_c=$(get_value "pre_A_c")  # Cores

pre_A_job() {
local job_name="pre_A_${ssl_id}"
local log_dir="./0_mainPRE/logs"
local job_command="export OMP_NUM_THREADS=1; Rscript ./0_mainPRE/pre_A.R"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$PRE_A_m" "$PRE_A_t" "$PRE_A_c" 1 "$job_command" "" "$log_dir"
check_exit "$job_name" $?
}
pre_A_job

##############################################
#        SESSION UPDATE                      #
##############################################
# Loop over species runs to prevent scratch path saturation
spe_runs="$(cat $wp/tmp/settings/ref_species_runs.txt)"

# Iterate through each run
for i in $(seq 1 "$spe_runs"); do
# Save the current run ID
echo "$i" > "$wp/tmp/settings/tmp_run_id.txt"
echo "Starting N-SDM run $i out of $spe_runs runs"

# Retrieve current time/date for logging
dt=$(date +"%FT%T")
echo "Run $i started at $dt"

# Update N-SDM settings using the R script
cd "$wp/scripts/"
Rscript ./0_mainPRE/nsdm_update.R >/dev/null 2>&1 || true

# Retrieve the number of species to model in this run
output_list="$wp/tmp/settings/tmp_species_list.txt"
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for this run: $n_spe"
check_species_count "$n_spe"

# Retrieve modelling algorithms, nesting methods, scenarios, and periods using get_value
mod_algo=$(get_value "mod_algo")  # Modelling algorithms evaluated
n_algo=$(count_items "$mod_algo") # Number of modelling algorithms

nest_met=$(get_value "nesting_methods")  # Nesting methods evaluated
n_nesting=$(count_items "$nest_met")     # Number of nesting methods

scenars=$(get_value "proj_scenarios")    # Projection scenarios evaluated
n_scenarios=$(count_items "$scenars")    # Number of projection scenarios

periods=$(get_value "proj_periods")      # Projection periods evaluated
n_periods=$(count_items "$periods")      # Number of projection periods

# Define file patterns to clean
file_patterns=("*.err" "*.out")

# Define log_dirs
log_dirs=(
"$wp/scripts/1_mainGLO/"
"$wp/scripts/2_mainREG/"
"$wp/scripts/3_mainSCE/"
)

# Loop through directories and file patterns to remove log files
for dir in "${log_dirs[@]}"; do
for pattern in "${file_patterns[@]}"; do
find "$dir" -type f -name "$pattern" -exec rm -f {} \; 2>/dev/null || true
done
done

# Clean scratch output folder if requested
clear_sop=$(get_value "clear_sop")
if [ "$clear_sop" = "TRUE" ] && [ "$new_sess" = "TRUE" ]; then
echo "Clearing scratch directories..."
rm -r "$sop/outputs/" 2>/dev/null || true
rm -r "$sop/tmp/" 2>/dev/null || true
fi

# Retrieve n_levels and do_proj values from settings.psv
n_levels=$(get_value "n_levels")  # Number of levels of analyses
do_proj=$(get_value "do_proj")    # Do scenario analyses?

##############################################
#        RUN PRE_B JOB                          
##############################################
# PRE_B Job
cd $wp/scripts
PRE_B_m=$(get_value "pre_B_m")
PRE_B_t=$(get_value "pre_B_t")
PRE_B_c=$(get_value "pre_B_c")

pre_B_job() {
local job_name="pre_B_${ssl_id}_${i}" 
local log_dir="./logs"
local job_command="export OMP_NUM_THREADS=1; Rscript pre_B.R"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$PRE_B_m" "$PRE_B_t" "$PRE_B_c" 1 "$job_command" "" "$log_dir"
check_exit "$job_name" $?
}
cd "$wp/scripts/0_mainPRE"
pre_B_job

# Check number of species processed after pre_B
save_dir="$sop/outputs/d0_datasets/base"
output_list="$wp/tmp/settings/tmp_species_list.txt"
find "$save_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" > "$output_list"
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for GLO: $n_spe"
check_species_count "$n_spe"

##############################################
#        SET GLO_A JOB                          
##############################################
cd "$wp/scripts/"
# GLO_A job
GLO_A_m=$(get_value "glo_A_m")
GLO_A_t=$(get_value "glo_A_t")
GLO_A_c=$(get_value "glo_A_c")
GLO_A_a=$n_spe

glo_A_job() {
local job_name="glo_A_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$GLO_A_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript glo_A.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$GLO_A_m" "$GLO_A_t" "$GLO_A_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET GLO_B JOB                          
##############################################
GLO_B_m=$(get_value "glo_B_m")
GLO_B_t=$(get_value "glo_B_t")
GLO_B_c=$(get_value "glo_B_c")
GLO_B_a=$((n_spe * n_algo))

glo_B_job() {
local job_name="glo_B_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$GLO_B_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript glo_B.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$GLO_B_m" "$GLO_B_t" "$GLO_B_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET GLO_C JOB                           
##############################################
GLO_C_m=$(get_value "glo_C_m")
GLO_C_t=$(get_value "glo_C_t")
GLO_C_c=$(get_value "glo_C_c")
GLO_C_a=$n_spe

glo_C_job() {
local job_name="glo_C_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$GLO_C_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript glo_C.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$GLO_C_m" "$GLO_C_t" "$GLO_C_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        RUN GLO JOBS                           
##############################################
cd "$wp/scripts/1_mainGLO"
glo_A_job
glo_B_job
glo_C_job

##############################################
#        SET REG_A JOB                           
##############################################
if [ $n_levels -gt 1 ]; then 
# Check number of species processed after glo_C
save_dir="$sop/outputs/d8_ensembles/glo"
output_list="$wp/tmp/settings/tmp_species_list.txt"
find "$save_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" > "$output_list"
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for REG: $n_spe"
check_species_count "$n_spe"

## REGIONAL level
cd "$wp/scripts/"
# REG_A job
REG_A_m=$(get_value "reg_A_m")
REG_A_t=$(get_value "reg_A_t")
REG_A_c=$(get_value "reg_A_c")
REG_A_a=$n_spe

reg_A_job() {
local job_name="reg_A_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$REG_A_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript reg_A.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$REG_A_m" "$REG_A_t" "$REG_A_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET REG_B JOB                           
##############################################
REG_B_m=$(get_value "reg_B_m")
REG_B_t=$(get_value "reg_B_t")
REG_B_c=$(get_value "reg_B_c")
REG_B_a=$((n_spe * n_algo * n_nesting))

reg_B_job() {
local job_name="reg_B_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$REG_B_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript reg_B.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$REG_B_m" "$REG_B_t" "$REG_B_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET REG_C JOB                           
##############################################
REG_C_m=$(get_value "reg_C_m")
REG_C_t=$(get_value "reg_C_t")
REG_C_c=$(get_value "reg_C_c")
REG_C_a=$n_spe

reg_C_job() {
local job_name="reg_C_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$REG_C_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript reg_C.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$REG_C_m" "$REG_C_t" "$REG_C_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        RUN REG JOBS                           
##############################################
cd "$wp/scripts/2_mainREG"
reg_A_job
reg_B_job
reg_C_job
fi

##############################################
#        SET SCE_A JOB                           
##############################################
if [ $do_proj = "TRUE" ]; then
# Check number of species processed after reg_C
save_dir="$sop/outputs/d8_ensembles/reg"
output_list="$wp/tmp/settings/tmp_species_list.txt"
find "$save_dir" -mindepth 2 -maxdepth 2 -type d -printf "%f\n" | sort | uniq > "$output_list"
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for SCE: $n_spe"
check_species_count "$n_spe"

## Projections
cd "$wp/scripts/"
# SCE_A job
SCE_A_m=$(get_value "sce_A_m")
SCE_A_t=$(get_value "sce_A_t")
SCE_A_c=$(get_value "sce_A_c")
SCE_A_a=$((n_spe * n_scenarios))

sce_A_job() {
local job_name="sce_A_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$SCE_A_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript sce_A.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$SCE_A_m" "$SCE_A_t" "$SCE_A_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET SCE_B JOB                           
##############################################
SCE_B_m=$(get_value "sce_B_m")
SCE_B_t=$(get_value "sce_B_t")
SCE_B_c=$(get_value "sce_B_c")
SCE_B_a=$n_spe

sce_B_job() {
local job_name="sce_B_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$SCE_B_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript sce_B.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$SCE_B_m" "$SCE_B_t" "$SCE_B_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        RUN SCE_A_B JOBS                           
##############################################
cd "$wp/scripts/3_mainSCE"
sce_A_job
sce_B_job

##############################################
#        SET SCE_C JOB                           
##############################################
if [ $n_levels -gt 1 ]; then
cd "$wp/scripts/"
# SCE_C job
SCE_C_m=$(get_value "sce_C_m")
SCE_C_t=$(get_value "sce_C_t")
SCE_C_c=$(get_value "sce_C_c")
SCE_C_a=$((n_spe * n_nesting * n_scenarios))

sce_C_job() {
local job_name="sce_C_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$SCE_C_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript sce_C.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$SCE_C_m" "$SCE_C_t" "$SCE_C_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        SET SCE_D JOB                           
##############################################
SCE_D_m=$(get_value "sce_D_m")
SCE_D_t=$(get_value "sce_D_t")
SCE_D_c=$(get_value "sce_D_c")
SCE_D_a=$((n_spe * n_nesting))

sce_D_job() {
local job_name="sce_D_${ssl_id}_${i}"
local log_dir="./logs"
local array_flag="[1-$SCE_D_a]"
local job_command="export OMP_NUM_THREADS=1; Rscript sce_D.R \$SLURM_ARRAY_TASK_ID"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$SCE_D_m" "$SCE_D_t" "$SCE_D_c" 1 "$job_command" "$array_flag" "$log_dir"
check_exit "$job_name" $?
}

##############################################
#        RUN SCE_C_D JOBS                           
##############################################
cd "$wp/scripts/3_mainSCE"
sce_C_job
sce_D_job
fi
fi

##############################################
#        SYNC TO SAVEPATH                           
##############################################
job_name="sync_files_${ssl_id}_${i}"
cd "$wp/scripts"

# Check if the job has already been completed
if job_completed "$job_name"; then 
echo "Job $job_name has already been completed successfully. Skipping..."
else
echo "Starting job $job_name..."

# Set permissions
chmod -R 777 "$wp/scripts/"

# Get the rsync exclude list
rsync_exclude=$(get_value "rsync_exclude" | tr ',' ' ')

# Check if "d2_models" is NOT in the exclude list to filter only final model objects
if [[ ! " $rsync_exclude " =~ " d2_models " ]]; then
cd "$sop/outputs/" || exit 1
find d2_models/ -name '*glm.rds' -o -name '*gam.rds' -o -name '*rf.rds' \
-o -name '*max.rds' -o -name '*gbm.rds' -o -name '*esm.rds' > "$wp/tmp/settings/tmp_modfiles.txt"

# Perform rsync only if tmp_modfiles.txt is not empty
[[ -s "$wp/tmp/settings/tmp_modfiles.txt" ]] && rsync -a --files-from="$wp/tmp/settings/tmp_modfiles.txt" . "$svp/outputs"
fi

# Generate exclusion list and sync excluding files
awk -F "|" '$1 == "rsync_exclude" { print $2 }' "$wp/scripts/settings/settings.psv" | tr ',' '\n' > "$wp/tmp/settings/tmp_exclfiles.txt"
rsync -a --exclude-from="$wp/tmp/settings/tmp_exclfiles.txt" "$sop/outputs/" "$svp/outputs"

# Sync sacct file
sacct_name="sacct_${ssl_id}_${i}.txt"
sacct_dir="$svp/outputs/sacct"
mkdir -p "$sacct_dir" && rsync -a "$wp/tmp/sacct/sacct_log.txt" "$sacct_dir/$sacct_name"

echo "Outputs synced to saving location."

# Create a checkpoint after successful completion
create_checkpoint "$job_name"
fi
done

# Inform that nsdm.sh is completed
echo "N-SDM completed."