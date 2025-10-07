##############################################
# FUNCTION: submit_and_monitor_job
# DESCRIPTION:
#   Submit and monitor a job on a SLURM-managed HPC cluster.
#   Automatically prepares an sbatch script, configures resources,
#   supports array jobs, manages logs, and checks for job failures.
#
# NOTES:
#   - Loads required modules including R.
#   - Creates timestamped .out/.err log files.
#   - Waits for job completion and verifies status.
#   - Appends job summary to $sacct_log if available.
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

    # Resolve the module list
    local module_r module_others
    module_r=$(get_value "module_r")
    module_others=$(get_value "module_others" | tr ',' ' ')

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

$(if [ -n "$module_others" ]; then
    for mod in $module_others; do
        echo "module load $mod"
    done
fi)

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