#############################################################################
## Script: nsdm_update
## Author: Antoine Adde
#############################################################################

### =========================================================================
### A. Initialization
### =========================================================================

# Load N-SDM settings
load(file.path(gsub("scripts", "tmp", getwd()), "settings", "ref_nsdm_settings.RData"))

# global reproducibility seed
set.seed(seed)

### =========================================================================
### B. Species Management
### =========================================================================

# Load species list
species <- readLines(file.path(w_path, "tmp", "settings", "ref_species_list.txt"))

# Split species into cluster run chunks to match n_mx_spe
splits <- split(seq_len(length(species)), ceiling(seq_len(length(species)) / n_mx_spe))

# Retrieve run ID
run_id <- readLines(file.path(w_path, "tmp", "settings", "tmp_run_id.txt"))

# Subset species for the current run
species <- species[splits[[run_id]]]

# Save the species list for the current run
writeLines(species, file.path(w_path, "tmp", "settings", "tmp_species_list.txt"))

### =========================================================================
### C. Update Settings
### =========================================================================

# Update N-SDM settings
save.image(file.path(w_path, "tmp", "settings", "tmp_nsdm_settings.RData"))