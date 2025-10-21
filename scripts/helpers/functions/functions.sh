# FUNCTION: retrieve values from settings.psv
get_value() {
  local key=$1
  local val
  val=$(awk -F"|" -v k="$key" '$1==k{print $2}' ./settings/settings.psv | tr -d '\r\t' | xargs)

  if [[ $key == "nesting_methods" ]]; then
    echo "$val" | tr ',' '\n' | grep -E 'covariate|multiply|average' \
      | sed 's/^\(multiply.*\|average.*\)$/posthoc/' | sort -u | paste -sd',' -
  else
    echo "$val"
  fi
}

# FUNCTION: Count the number of items in a list
count_items() {
    local list="$1"
    list=$(echo "$list" | tr -d '\r\t' | sed 's/^ *//;s/ *$//;s/ *, */,/g;s/,,*/,/g')
    if [[ -z "$list" ]]; then
        echo 0
        return
    fi
    echo "$list" | tr ',' '\n' | grep -v '^[[:space:]]*$' | wc -l
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
    echo "No species found, exiting."
    exit 1
  fi
}

# FUNCTION: Exit on any error (based on settings)
exit_on_error=$(get_value "exit_any")
if [ "$exit_on_error" = "TRUE" ]; then
  set -e
  echo "[INFO] Exiting N-SDM on any error is ENABLED (setting exit_any=TRUE)"
else
  echo "[INFO] Exiting N-SDM on any error is DISABLED (setting exit_any=FALSE)"
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
    if [ "$exit_on_error" = "TRUE" ] && [ $exit_code -ne 0 ]; then
        echo "An error occurred. Exiting."
    fi
}