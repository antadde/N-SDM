#!/bin/bash

##############################################
# FUNCTION: Count the number of items in a list
##############################################
count_items() {
    local list=$1
    echo $(( $(grep -o "'" <<< "$list" | wc -l) / 2 ))
}

##############################################
# FUNCTION: Create a checkpoint for a job
##############################################
create_checkpoint() {
    local job_name=$1
    local i=$2
    local ssl_id=$3
    echo "$ssl_id $i" > "$checkpoint_dir/${job_name}_done"
}

##############################################
# FUNCTION: Check if a job has already been completed
##############################################
job_completed() {
    local job_name=$1
    [ -f "$checkpoint_dir/${job_name}_done" ]
}

##############################################
# FUNCTION: Submit a job to SLURM
##############################################
submit_and_monitor_job() {
    local job_name=$1
    local mem=$2
    local time=$3
    local cpus=$4
    local ntasks=$5
    local job_command=$6
    local array_flag=$7
    local log_dir=$8

# Create timestamp
    local timestamp=$(date +%Y%m%d_%H%M%S)

    # Customize output and error filenames
    if [ -z "$array_flag" ]; then
        local output_file="${log_dir}/${job_name}_$timestamp.out"
        local error_file="${log_dir}/${job_name}_$timestamp.err"
    else
        local output_file="${log_dir}/${job_name}_%A_%a_$timestamp.out"
        local error_file="${log_dir}/${job_name}_%A_%a_$timestamp.err"
    fi

    echo "Submitting job: $job_name..."
    local job_id
    if [ -z "$array_flag" ]; then
        job_id=$(sbatch --wait --job-name="$job_name" --output="$output_file" --error="$error_file" \
                      --mem-per-cpu="$mem" --time="$time" --cpus-per-task="$cpus" --ntasks="$ntasks" --wrap="$job_command" | awk '{print $NF}')
    else
        job_id=$(sbatch --wait --job-name="$job_name" --output="$output_file" --error="$error_file" \
                      --mem-per-cpu="$mem" --time="$time" --cpus-per-task="$cpus" --ntasks="$ntasks" --array="$array_flag" --wrap="$job_command" | awk '{print $NF}')
    fi

    # Wait a moment for job to register in sacct
    sleep 5

# Fetch the job status from sacct (remove empty lines and consolidate)
job_status=$(sacct -j "$job_id" --format=State --noheader | awk '{$1=$1};1')

# Store job failure states if present (but prevent grep from exiting)
failed_status=$(echo "$job_status" | grep -E 'FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|OUT_OF_ME|NODE_FAIL' || true)

# Print the full job status
echo "Job $job_name ended."

# If any failure is detected, return non-zero to caller
if [[ -n "$failed_status" ]]; then
    echo "Error: Detected a failure in job $job_name"
    return 1
fi

    # Check if the log file already exists and is not empty
    if [[ ! -s "$sacct_log" ]]; then
        # If the file does not exist or is empty, print the header
        sacct_summary "$job_id" >> "$sacct_log"
    else
        # If the file exists and is not empty, append only the data without the header
        sacct_summary "$job_id" | tail -n +2 >> "$sacct_log"
    fi
}


##############################################
# FUNCTION: SACCT SUMMARY
#############################################
sacct_summary() {
    if [ -z "$1" ]; then
        echo "Usage: sacct_summary <JobID>"
        return 1
    fi

    job_id="$1"

    sacct -j "$job_id" --format=JobID,JobName,State,Elapsed,TotalCPU,ReqCPUS,Timelimit,MaxRSS,NNodes,NTasks,ReqMem --parsable2 | awk -F '|' '
    NR==1 { print "JobID|JobName|State|Elapsed|Timelimit|TotalCPU|ReqCPUS|UsedCPUS|MaxRSS|NNodes|NTasks|ReqMem"; next }
    $2 !~ /batch|extern/ { jobname=$2; reqmem=$11; timelimit=$7; reqcpus=$6 }  # Store main job info
    $2 == "batch" {
        sub(/\.batch/, "", $1)

        mem_kb = ($8+0)
        mem_gb = (mem_kb / (1024^2))
        mem_rec = (mem_gb * 1.2) * reqcpus

        split($4, elapsed, ":")
        if (length(elapsed) == 3) { elapsed_sec = elapsed[1] * 3600 + elapsed[2] * 60 + elapsed[3] }
        else if (length(elapsed) == 2) { elapsed_sec = elapsed[1] * 60 + elapsed[2] }
        else { elapsed_sec = elapsed[1] }

        split($5, totalcpu, ":")
        if (length(totalcpu) == 3) { totalcpu_sec = totalcpu[1] * 3600 + totalcpu[2] * 60 + totalcpu[3] }
        else if (length(totalcpu) == 2) { totalcpu_sec = totalcpu[1] * 60 + totalcpu[2] }
        else { totalcpu_sec = totalcpu[1] }

        if (elapsed_sec > 0) {
            used_cpus = totalcpu_sec / elapsed_sec
        } else {
            used_cpus = 0
        }

        split(timelimit, reqtime_arr, ":")
        if (length(reqtime_arr) == 3) { req_sec = reqtime_arr[1] * 3600 + reqtime_arr[2] * 60 + reqtime_arr[3] }
        else if (length(reqtime_arr) == 2) { req_sec = reqtime_arr[1] * 60 + reqtime_arr[2] }
        else { req_sec = reqtime_arr[1] }

        rec_sec = int(elapsed_sec * 1.2)
        rec_hh = int(rec_sec / 3600)
        rec_mm = int((rec_sec % 3600) / 60)
        rec_ss = rec_sec % 60
        recommended_time = sprintf("%02d:%02d:%02d", rec_hh, rec_mm, rec_ss)

        printf "%s|%s|%s|%s|%s|%s|%s|%.1f|%s|%s|%s|%s\n", $1, jobname, $3, $4, timelimit, $5, reqcpus, used_cpus, $8, $9, $10, reqmem
    }'
}

##############################################
# FUNCTION: OTHERS
#############################################
# Function: retrieve values from settings.psv
get_value() {
  local key=$1
  awk -F "|" -v search_key="$key" '$1 == search_key { print $2 }' ./settings/settings.psv
}

# Exit if species count is zero
check_species_count() {
  if [ "$1" -eq 0 ]; then
    echo "No species found."
    exit 1
  fi
}

# Exit on any error (based on settings)
exit_on_error=$(get_value "exit_any")
if [ "$exit_on_error" = "TRUE" ]; then
  set -e
  echo "[INFO] Exiting on any error is ENABLED (exit_any=TRUE)"
else
  echo "[INFO] Exiting on any error is DISABLED (exit_any=$exit_on_error)"
fi

# Check exiting status
check_exit() {
    local job_name=$1
    local exit_code=$2
    if [ $exit_code -ne 0 ]; then
        echo "Job $job_name failed."
        if [ "$exit_on_error" = "TRUE" ]; then
            echo "Exiting due to exit_any=TRUE"
            exit $exit_code
        fi
    fi
}

# Clean up on script exit
cleanup() {
    local exit_code=$?
    # Always do cleanup actions here (e.g. remove temp files)

    if [ "$exit_on_error" = "TRUE" ] && [ $exit_code -ne 0 ]; then
        echo "An error occurred. Cleaning up and exiting."
    fi
}

