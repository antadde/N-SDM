#!/bin/bash

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

    # Resolve the R module BEFORE writing the script
    local module_r
    module_r=$(get_value "module_r")

    # Create timestamp for logs
    local timestamp
    timestamp=$(date +%Y%m%d_%H%M%S)

    # Log files
    local output_file error_file
    if [ -z "$array_flag" ]; then
        output_file="${log_dir}/${job_name}_${timestamp}.out"
        error_file="${log_dir}/${job_name}_${timestamp}.err"
    else
        output_file="${log_dir}/${job_name}_%A_%a_${timestamp}.out"
        error_file="${log_dir}/${job_name}_%A_%a_${timestamp}.err"
    fi

    mkdir -p "$log_dir"

    # Temporary job script
    local job_script="${log_dir}/${job_name}_${timestamp}.sbatch"
	
	# Decide what command will actually be written to the script
    local final_command="$job_command"
    if [ -n "$array_flag" ]; then
        final_command="$job_command \$SLURM_ARRAY_TASK_ID"
    fi
	
    cat > "$job_script" <<EOF
#!/bin/bash -l
#SBATCH --job-name=$job_name
#SBATCH --output=$output_file
#SBATCH --error=$error_file
#SBATCH --mem-per-cpu=$mem
#SBATCH --time=$time
#SBATCH --cpus-per-task=$cpus
#SBATCH --ntasks=$ntasks
${array_flag:+#SBATCH --array=$array_flag}

module load $module_r
export OMP_NUM_THREADS=1

$final_command
EOF

    echo "Submitting job: $job_name..."
    local job_id
    job_id=$(sbatch --wait "$job_script" | awk '{print $NF}')

    sleep 5

    job_status=$(sacct -j "$job_id" --format=State --noheader | awk '{$1=$1};1')
    failed_status=$(echo "$job_status" | grep -E 'FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|OUT_OF_ME|NODE_FAIL' || true)

    echo "Job $job_name finished."

    if [[ -n "$failed_status" ]]; then
        echo "Error: Detected a failure in job $job_name"
        return 1
    fi

    if [[ -n "$sacct_log" ]]; then
        if [[ ! -s "$sacct_log" ]]; then
            sacct_summary "$job_id" >> "$sacct_log"
        else
            sacct_summary "$job_id" | tail -n +2 >> "$sacct_log"
        fi
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

        printf "%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n", $1, jobname, $3, $4, timelimit, $5, reqcpus, $8, $9, $10, reqmem
    }'
}

#############################################
# OTHERS
#############################################

# FUNCTION: retrieve values from settings.psv
get_value() {
  local key=$1
  awk -F "|" -v search_key="$key" '$1 == search_key { print $2 }' ./settings/settings.psv
}

# FUNCTION: Count the number of items in a list
count_items() {
    local list=$1
    echo $(( $(grep -o "'" <<< "$list" | wc -l) / 2 ))
}

# FUNCTION: Create a checkpoint for a job
create_checkpoint() {
    local job_name=$1
    local i=$2
    local ssl_id=$3
    echo "$ssl_id $i" > "$checkpoint_dir/${job_name}_done"
}

# FUNCTION: Check if a job has already been completed
job_completed() {
    local job_name=$1
    [ -f "$checkpoint_dir/${job_name}_done" ]
}

# FUNCTION: Exit if species count is zero
check_species_count() {
  if [ "$1" -eq 0 ]; then
    echo "No species found."
    exit 1
  fi
}

# FUNCTION: Exit on any error (based on settings)
exit_on_error=$(get_value "exit_any")
if [ "$exit_on_error" = "TRUE" ]; then
  set -e
  echo "[INFO] Exiting on any error is ENABLED (exit_any=TRUE)"
else
  echo "[INFO] Exiting on any error is DISABLED (exit_any=$exit_on_error)"
fi

# FUNCTION: Check exiting status
check_exit() {
    local job_name=$1
    local exit_code=$2
    local i=$3
    local ssl_id=$4

    if [ $exit_code -eq 0 ]; then
        create_checkpoint "$job_name" "$i" "$ssl_id"
    else
        echo "Job $job_name failed."
        if [ "$exit_on_error" = "TRUE" ]; then
            echo "Exiting due to exit_any=TRUE"
            exit $exit_code
        fi
    fi
}

# FUNCTION: Clean up on script exit
cleanup() {
    local exit_code=$?
    # Always do cleanup actions here (e.g. remove temp files)

    if [ "$exit_on_error" = "TRUE" ] && [ $exit_code -ne 0 ]; then
        echo "An error occurred. Cleaning up and exiting."
    fi
}