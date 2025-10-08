##############################################
# FUNCTION: submit_and_monitor_job
# DESCRIPTION:
#   Submit and monitor a job on a SLURM-managed HPC cluster.
#   Automatically prepares an sbatch script, configures resources,
#   supports array jobs, manages logs, and checks for job failures.
##############################################
submit_and_monitor_job() {
    local job_name=$1 mem=$2 time=$3 cpus=$4 ntasks=$5 partition=$6 job_command=$7 array_flag=$8 log_dir=$9
    local template_path="./helpers/job_template.sbatch.in"
    local timestamp job_script output_file error_file module_r module_others

    module_r=$(get_value "module_r")
    module_others=$(get_value "module_others" | tr ',' ' ')
    timestamp=$(date +%Y%m%d_%H%M%S)
    mkdir -p "$log_dir"

    if [ -z "$array_flag" ]; then
        output_file="${log_dir}/${job_name}_${timestamp}.out"
        error_file="${log_dir}/${job_name}_${timestamp}.err"
    else
        output_file="${log_dir}/${job_name}_%A_%a_${timestamp}.out"
        error_file="${log_dir}/${job_name}_%A_%a_${timestamp}.err"
    fi
    job_script="${log_dir}/${job_name}_${timestamp}.sbatch"

    ##########################################################
    # UNIVERSAL MEMORY MODE DETECTION (cached)
    ##########################################################
    local mem_mode_cache="/tmp/slurm_mem_mode_detected"
    if [[ ! -f "$mem_mode_cache" ]]; then
        local test_script
        test_script=$(mktemp)
        echo -e "#!/bin/bash\n#SBATCH --mem=100M\necho OK" > "$test_script"
        sbatch --test-only "$test_script" 2>&1 | grep -q "not supported" \
            && echo "percpu" > "$mem_mode_cache" || echo "node" > "$mem_mode_cache"
        rm -f "$test_script"
    fi

    local SLURM_MEM_MODE mem_directive
    SLURM_MEM_MODE=$(<"$mem_mode_cache")
    if [[ "$SLURM_MEM_MODE" == "percpu" ]]; then
        local mem_value mem_unit per_cpu_mem
        mem_value=$(echo "$mem" | grep -oE '[0-9]+')
        mem_unit=$(echo "$mem" | grep -oE '[A-Za-z]+')
        per_cpu_mem=$(( (mem_value + cpus - 1) / cpus ))
        mem_directive="#SBATCH --mem-per-cpu=${per_cpu_mem}${mem_unit}"
    else
        mem_directive="#SBATCH --mem=$mem"
    fi

    ##########################################################
    # SLURM DIRECTIVES
    ##########################################################
    local partition_line="" array_line=""
    [[ -n "$partition" ]] && partition_line="#SBATCH --partition=$partition"
    [[ -n "$array_flag" ]] && array_line="#SBATCH --array=$array_flag"

    ##########################################################
    # COMMAND HANDLING (pass array ID if needed)
    ##########################################################
    local final_command="$job_command"
    [[ -n "$array_flag" ]] && final_command="$job_command \$SLURM_ARRAY_TASK_ID"

    ##########################################################
    # LOAD TEMPLATE (fallback if missing)
    ##########################################################
    local template_content
    if [[ -f "$template_path" ]]; then
        template_content=$(<"$template_path")
    else
        read -r -d '' template_content <<'EOF'
#!/bin/bash -l
#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{OUTPUT_FILE}}
#SBATCH --error={{ERROR_FILE}}
#SBATCH --time={{TIME}}
#SBATCH --cpus-per-task={{CPUS}}
#SBATCH --ntasks={{NTASKS}}
{{PARTITION_LINE}}
{{ARRAY_LINE}}
{{MEM_LINE}}

{{MODULE_OTHERS}}
module load {{MODULE_R}}
export OMP_NUM_THREADS=1

{{COMMAND}}
EOF
    fi

    ##########################################################
    # SUBSTITUTE PLACEHOLDERS
    ##########################################################
    local module_lines=""
    if [ -n "$module_others" ]; then
        for mod in $module_others; do
            module_lines+="module load $mod"$'\n'
        done
    fi

    awk -v job_name="$job_name" \
        -v output_file="$output_file" \
        -v error_file="$error_file" \
        -v time="$time" \
        -v cpus="$cpus" \
        -v ntasks="$ntasks" \
        -v partition_line="$partition_line" \
        -v array_line="$array_line" \
        -v mem_line="$mem_directive" \
        -v module_others="$module_lines" \
        -v module_r="$module_r" \
        -v command="$final_command" '
        { gsub(/{{JOB_NAME}}/,job_name);
          gsub(/{{OUTPUT_FILE}}/,output_file);
          gsub(/{{ERROR_FILE}}/,error_file);
          gsub(/{{TIME}}/,time);
          gsub(/{{CPUS}}/,cpus);
          gsub(/{{NTASKS}}/,ntasks);
          gsub(/{{PARTITION_LINE}}/,partition_line);
          gsub(/{{ARRAY_LINE}}/,array_line);
          gsub(/{{MEM_LINE}}/,mem_line);
          gsub(/{{MODULE_OTHERS}}/,module_others);
          gsub(/{{MODULE_R}}/,module_r);
          gsub(/{{COMMAND}}/,command);
          print }' <<< "$template_content" > "$job_script"

    ##########################################################
    # SUBMIT AND MONITOR
    ##########################################################
    echo "Submitting job: $job_name..."
    local job_id
    job_id=$(sbatch --wait "$job_script" | awk '{print $NF}')
    sleep 5

    local job_status failed_status
    job_status=$(sacct -j "$job_id" --format=State --noheader | awk '{$1=$1};1')
    failed_status=$(grep -E 'FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|OUT_OF_ME|NODE_FAIL' <<< "$job_status" || true)

    echo "Job $job_name finished."
    if [[ -n "$failed_status" ]]; then
        echo "Error: Detected a failure in job $job_name"
        return 1
    fi

    if [[ -n "$sacct_log" ]]; then
        sacct_summary "$job_id" >> "$sacct_log"
    fi
}