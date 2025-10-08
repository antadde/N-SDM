# ============================================================
# FUNCTION: sync_outputs
# Description:
#   Sync model outputs and logs to the saving directory ($svp)
#   using rsync and dynamically generated include/exclude lists.
# Arguments:
#   None (uses global vars: wp, sop, svp, ssl_id, i)
# ============================================================
sync_outputs() {
  # Ensure scripts are executable
  chmod -R 777 "$wp/scripts/"

  # Retrieve rsync exclude list from settings.psv
  rsync_exclude=$(get_value "rsync_exclude" | tr ',' ' ')

  # If d2_models is NOT excluded, filter final model objects
  if [[ ! " $rsync_exclude " =~ " d2_models " ]]; then
    cd "$sop/outputs/" || { echo "Error: cannot cd to $sop/outputs/"; return 1; }

    find d2_models/ \( -name '*glm.rds' -o -name '*gam.rds' -o -name '*rf.rds' \
      -o -name '*max.rds' -o -name '*gbm.rds' -o -name '*esm.rds' \) \
      > "$wp/tmp/settings/tmp_modfiles.txt"

    # Perform rsync only if tmp_modfiles.txt is not empty
    if [[ -s "$wp/tmp/settings/tmp_modfiles.txt" ]]; then
      rsync -a --files-from="$wp/tmp/settings/tmp_modfiles.txt" . "$svp/outputs"
    fi
  fi

  # Generate exclusion list and sync excluding files
  awk -F "|" '$1 == "rsync_exclude" { print $2 }' "$wp/scripts/settings/settings.psv" \
    | tr ',' '\n' > "$wp/tmp/settings/tmp_exclfiles.txt"

  rsync -a --exclude-from="$wp/tmp/settings/tmp_exclfiles.txt" "$sop/outputs/" "$svp/outputs"

  # Sync sacct log
  sacct_name="sacct_${ssl_id}_${i}.psv"
  sacct_dir="$svp/outputs/sacct"
  mkdir -p "$sacct_dir"
  rsync -a "$wp/tmp/sacct/sacct_log.txt" "$sacct_dir/$sacct_name"

  echo "Outputs synced to saving location: $svp/outputs"
}
