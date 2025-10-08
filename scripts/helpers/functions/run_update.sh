# ============================================================
# FUNCTION: run_update
# Description:
#   Split the reference species list into chunks according to
#   n_mx_spe and current run ID, then save the subset to tmp_species_list.txt
# Arguments:
#   None (uses global variable, and helper function get_value)
# ============================================================
run_update() {
  # --- Inputs ---
  local species_file="${wp}/tmp/settings/ref_species_list.txt"
  local run_id_file="${wp}/tmp/settings/tmp_run_id.txt"
  local out_file="${wp}/tmp/settings/tmp_species_list.txt"

  # --- Safety checks ---
  [[ -f "$species_file" ]] || { echo "Error: species file not found: $species_file"; return 1; }
  [[ -f "$run_id_file" ]] || { echo "Error: run ID file not found: $run_id_file"; return 1; }

  # --- Read run ID ---
  local run_id
  run_id=$(cat "$run_id_file")

  # --- Read species into array ---
  local species
  mapfile -t species < "$species_file"

  # --- Compute split boundaries ---
  local n_mx_spe
  n_mx_spe=$(get_value "n_mx_spe")
  local total=${#species[@]}
  local chunks=$(( (total + n_mx_spe - 1) / n_mx_spe ))
  local start=$(( (run_id - 1) * n_mx_spe ))
  local end=$(( run_id * n_mx_spe - 1 ))
  (( end >= total )) && end=$(( total - 1 ))
  
  # --- Extract subset and write to file ---
	> "$out_file"
	local i
	for ((i=start; i<=end; i++)); do
		echo "${species[$i]}" >> "$out_file"
	done
}