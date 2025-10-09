#!/bin/bash

##############################################
##############################################
## nsdm.sh - Core N-SDM script
## Author: Antoine Adde (antoine.adde@eawag.ch)
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
# Load helper functions
for f in ./helpers/functions/*.sh; do
  [ -f "$f" ] && source "$f"
done

# Run cleanup function on script exit
trap cleanup EXIT

# Retrieve main paths from settings.psv
wp=$(get_value "w_path")    # Working path
sop=$(get_value "scr_path")  # Scratch output path
svp=$(get_value "svp_path")  # Saving output path
rlibs=$(get_value "lib_path")  # Saving output path

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
rm "$wp/scripts/0_mainPRE/logs/"* 2>/dev/null || true
fi

##############################################
#        RUN PRE_A JOB                        
##############################################
PRE_A_m=$(get_value "pre_A_m")  # Memory
PRE_A_t=$(get_value "pre_A_t")  # Time
PRE_A_c=$(get_value "pre_A_c")  # Cores
PRE_A_p=$(get_value "pre_A_p")  # Partition (optional)

pre_A_job() {
local job_name="pre_A_${ssl_id}"
local log_dir="./0_mainPRE/logs"
local job_command="Rscript ./0_mainPRE/pre_A.R"

# Check if the job has already been completed
if job_completed "$job_name"; then
echo "$job_name has already been completed successfully. Skipping..."
return 0
fi

# Create logs directory if it doesn't exist
mkdir -p "$log_dir"

# Submit the job
submit_and_monitor_job "$job_name" "$PRE_A_m" "$PRE_A_t" "$PRE_A_c" 1 "$PRE_A_p" "$job_command" "" "$log_dir"
check_exit "$job_name" $? "$i" "$ssl_id"
}
pre_A_job

##############################################
#        SESSION UPDATE                      #
##############################################
# Loop over species runs to prevent scratch path saturation
spe_runs="$(cat $wp/tmp/settings/ref_species_runs.txt)"

# Iterate through each run
for i in $(seq 1 "$spe_runs"); do

cd "$wp/scripts/"

# Save the current run ID
echo "$i" > "$wp/tmp/settings/tmp_run_id.txt"
echo "Starting N-SDM run $i out of $spe_runs runs"

# Retrieve current time/date for logging
dt=$(date +"%FT%T")
echo "Run $i started at $dt"

# Update N-SDM settings
run_update

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
file_patterns=("*.err" "*.out" "*.sbatch")

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
clear_sop=$(get_value "clear_scr")
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
PRE_B_m=$(get_value "pre_B_m")  # Memory
PRE_B_t=$(get_value "pre_B_t")  # Time
PRE_B_c=$(get_value "pre_B_c")  # Cores
PRE_B_p=$(get_value "pre_B_p")  # Partition (optional)

pre_B_job() {
    local job_name="pre_B_${ssl_id}_${i}" 
    local log_dir="./0_mainPRE/logs"
    local job_command="Rscript ./0_mainPRE/pre_B.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$PRE_B_m" "$PRE_B_t" "$PRE_B_c" 1 "$PRE_B_p" "$job_command" "" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}
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
GLO_A_m=$(get_value "glo_A_m")  # Memory
GLO_A_t=$(get_value "glo_A_t")  # Time
GLO_A_c=$(get_value "glo_A_c")  # Cores
GLO_A_p=$(get_value "glo_A_p")  # Partition (optional)
GLO_A_a=$n_spe                  # Array size

glo_A_job() {
    local job_name="glo_A_${ssl_id}_${i}"
    local log_dir="./1_mainGLO/logs"
    local array_flag="[1-$GLO_A_a]"
    local job_command="Rscript ./1_mainGLO/glo_A.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$GLO_A_m" "$GLO_A_t" "$GLO_A_c" 1 "$GLO_A_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET GLO_B JOB                          
##############################################
GLO_B_m=$(get_value "glo_B_m")  # Memory
GLO_B_t=$(get_value "glo_B_t")  # Time
GLO_B_c=$(get_value "glo_B_c")  # Cores
GLO_B_p=$(get_value "glo_B_p")  # Partition (optional)
GLO_B_a=$((n_spe * n_algo))     # Array size

glo_B_job() {
    local job_name="glo_B_${ssl_id}_${i}"
    local log_dir="./1_mainGLO/logs"
    local array_flag="[1-$GLO_B_a]"
    local job_command="Rscript ./1_mainGLO/glo_B.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$GLO_B_m" "$GLO_B_t" "$GLO_B_c" 1 "$GLO_B_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET GLO_C JOB                           
##############################################
GLO_C_m=$(get_value "glo_C_m")  # Memory
GLO_C_t=$(get_value "glo_C_t")  # Time
GLO_C_c=$(get_value "glo_C_c")  # Cores
GLO_C_p=$(get_value "glo_C_p")  # Partition (optional)
GLO_C_a=$n_spe                  # Array size

glo_C_job() {
    local job_name="glo_C_${ssl_id}_${i}"
    local log_dir="./1_mainGLO/logs"
    local array_flag="[1-$GLO_C_a]"
    local job_command="Rscript ./1_mainGLO/glo_C.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$GLO_C_m" "$GLO_C_t" "$GLO_C_c" 1 "$GLO_C_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        RUN GLO JOBS                           
##############################################
glo_A_job
glo_B_job
glo_C_job

## REGIONAL level
if [ $n_levels -gt 1 ]; then 
# Check number of species processed after glo_C
save_dir="$sop/outputs/d8_ensembles/glo"
output_list="$wp/tmp/settings/tmp_species_list.txt"
find "$save_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" > "$output_list"
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for REG: $n_spe"
check_species_count "$n_spe"

##############################################
#        SET REG_A JOB                        
##############################################
REG_A_m=$(get_value "reg_A_m")  # Memory
REG_A_t=$(get_value "reg_A_t")  # Time
REG_A_c=$(get_value "reg_A_c")  # Cores
REG_A_p=$(get_value "reg_A_p")  # Partition (optional)
REG_A_a=$n_spe                  # Array size

reg_A_job() {
    local job_name="reg_A_${ssl_id}_${i}"
    local log_dir="./2_mainREG/logs"
    local array_flag="[1-$REG_A_a]"
    local job_command="Rscript ./2_mainREG/reg_A.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$REG_A_m" "$REG_A_t" "$REG_A_c" 1 "$REG_A_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET REG_B JOB                           
##############################################
REG_B_m=$(get_value "reg_B_m")  # Memory
REG_B_t=$(get_value "reg_B_t")  # Time
REG_B_c=$(get_value "reg_B_c")  # Cores
REG_B_p=$(get_value "reg_B_p")  # Partition (optional)
REG_B_a=$((n_spe * n_algo * n_nesting))  # Array size

reg_B_job() {
    local job_name="reg_B_${ssl_id}_${i}"
    local log_dir="./2_mainREG/logs"
    local array_flag="[1-$REG_B_a]"
    local job_command="Rscript ./2_mainREG/reg_B.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$REG_B_m" "$REG_B_t" "$REG_B_c" 1 "$REG_B_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET REG_C JOB                           
##############################################
REG_C_m=$(get_value "reg_C_m")  # Memory
REG_C_t=$(get_value "reg_C_t")  # Time
REG_C_c=$(get_value "reg_C_c")  # Cores
REG_C_p=$(get_value "reg_C_p")  # Partition (optional)
REG_C_a=$n_spe                  # Array size

reg_C_job() {
    local job_name="reg_C_${ssl_id}_${i}"
    local log_dir="./2_mainREG/logs"
    local array_flag="[1-$REG_C_a]"
    local job_command="Rscript ./2_mainREG/reg_C.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$REG_C_m" "$REG_C_t" "$REG_C_c" 1 "$REG_C_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        RUN REG JOBS                           
##############################################
reg_A_job
reg_B_job
reg_C_job
fi

## Scenarios
if [ $do_proj = "TRUE" ]; then
  # Check number of species processed after reg_C
  if [ "$n_levels" -gt 1 ]; then
    save_dir="$sop/outputs/d8_ensembles/reg"
    output_list="$wp/tmp/settings/tmp_species_list.txt"
    find "$save_dir" -mindepth 2 -maxdepth 2 -type d -printf "%f\n" | sort | uniq > "$output_list"
  else
    # glo_C if single level
    save_dir="$sop/outputs/d8_ensembles/glo"
	output_list="$wp/tmp/settings/tmp_species_list.txt"
    find "$save_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort | uniq > "$output_list"
  fi
n_spe=$(wc -l < "$output_list")
echo "Number of species considered for SCE: $n_spe"
check_species_count "$n_spe"

##############################################
#        SET SCE_A JOB                           
##############################################
SCE_A_m=$(get_value "sce_A_m")  # Memory
SCE_A_t=$(get_value "sce_A_t")  # Time
SCE_A_c=$(get_value "sce_A_c")  # Cores
SCE_A_p=$(get_value "sce_A_p")  # Partition (optional)
SCE_A_a=$((n_spe * n_scenarios))  # Array size

sce_A_job() {
    local job_name="sce_A_${ssl_id}_${i}"
    local log_dir="./3_mainSCE/logs"
    local array_flag="[1-$SCE_A_a]"
    local job_command="Rscript ./3_mainSCE/sce_A.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$SCE_A_m" "$SCE_A_t" "$SCE_A_c" 1 "$SCE_A_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET SCE_B JOB                           
##############################################
SCE_B_m=$(get_value "sce_B_m")  # Memory
SCE_B_t=$(get_value "sce_B_t")  # Time
SCE_B_c=$(get_value "sce_B_c")  # Cores
SCE_B_p=$(get_value "sce_B_p")  # Partition (optional)
SCE_B_a=$n_spe                  # Array size

sce_B_job() {
    local job_name="sce_B_${ssl_id}_${i}"
    local log_dir="./3_mainSCE/logs"
    local array_flag="[1-$SCE_B_a]"
    local job_command="Rscript ./3_mainSCE/sce_B.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$SCE_B_m" "$SCE_B_t" "$SCE_B_c" 1 "$SCE_B_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        RUN SCE_A_B JOBS                           
##############################################
sce_A_job
sce_B_job


if [ $n_levels -gt 1 ]; then
##############################################
#        SET SCE_C JOB                           
##############################################
SCE_C_m=$(get_value "sce_C_m")  # Memory
SCE_C_t=$(get_value "sce_C_t")  # Time
SCE_C_c=$(get_value "sce_C_c")  # Cores
SCE_C_p=$(get_value "sce_C_p")  # Partition (optional)
SCE_C_a=$((n_spe * n_nesting * n_scenarios))  # Array size

sce_C_job() {
    local job_name="sce_C_${ssl_id}_${i}"
    local log_dir="./3_mainSCE/logs"
    local array_flag="[1-$SCE_C_a]"
    local job_command="Rscript ./3_mainSCE/sce_C.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$SCE_C_m" "$SCE_C_t" "$SCE_C_c" 1 "$SCE_C_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        SET SCE_D JOB                           
##############################################
SCE_D_m=$(get_value "sce_D_m")  # Memory
SCE_D_t=$(get_value "sce_D_t")  # Time
SCE_D_c=$(get_value "sce_D_c")  # Cores
SCE_D_p=$(get_value "sce_D_p")  # Partition (optional)
SCE_D_a=$((n_spe * n_nesting))  # Array size

sce_D_job() {
    local job_name="sce_D_${ssl_id}_${i}"
    local log_dir="./3_mainSCE/logs"
    local array_flag="[1-$SCE_D_a]"
    local job_command="Rscript ./3_mainSCE/sce_D.R"

    # Check if the job has already been completed
    if job_completed "$job_name"; then
        echo "$job_name has already been completed successfully. Skipping..."
        return 0
    fi

    # Create logs directory if it doesn't exist
    mkdir -p "$log_dir"

    # Submit the job
    submit_and_monitor_job "$job_name" "$SCE_D_m" "$SCE_D_t" "$SCE_D_c" 1 "$SCE_D_p" "$job_command" "$array_flag" "$log_dir"
    check_exit "$job_name" $? "$i" "$ssl_id"
}

##############################################
#        RUN SCE_C_D JOBS                           
##############################################
sce_C_job
sce_D_job
fi
fi

##############################################
#        SYNC TO SAVEPATH                           
##############################################
job_name="sync_files_${ssl_id}_${i}"

# Check if the job has already been completed
if job_completed "$job_name"; then 
echo "Job $job_name has already been completed successfully. Skipping..."
else
echo "Starting job $job_name..."

# Call the sync function
sync_outputs

# Create a checkpoint after successful completion
create_checkpoint "$job_name"
fi
done

# Inform that nsdm.sh is completed
echo "N-SDM completed."