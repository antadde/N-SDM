#############################################################################
## Script: 2_mainREG_C
## Purpose: Ensembling and mapping
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
# Load N-SDM settings
load(file.path(gsub("scripts", "tmp", gsub("/2_mainREG", "", getwd())), "settings", "tmp_nsdm_settings.RData"))

# Set permissions for new files
Sys.umask(mode = "000")

# Set working directory
setwd(w_path)

# Set library path
.libPaths(lib_path)

# Load N-SDM package
require(nsdm2)

### =========================================================================
### B- Definitions
### =========================================================================
# SBATCH param
args <- commandArgs(trailingOnly = TRUE)
arrayID <- as.numeric(args[1])

# Target species
species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

# Scale-nesting methods for combining GLO and REG predictions
nesting_methods<-nesting_methods

# SBATCH
ispi_name <- species[arrayID]

cat(paste0("Starting mapping and ensembling of REG predictions for ", ispi_name, "...\n"))

# Loop on nesting methods
for(nesting_method in nesting_methods){
	
cat(paste0("Starting mapping and ensembling of REG predictions for ", ispi_name, " with the ", nesting_method, " nesting method...\n"))
		   
### =========================================================================
### C- Save prediction raster
### =========================================================================

for(i in seq_along(mod_algo)){
  # Load raw prediction data
  model_name <- mod_algo[i]
  pred_path <- file.path(scr_path, "outputs", "d6_preds/reg", nesting_method)
  full_pred_path <- file.path(pred_path, ispi_name, model_name)
  pred_file <- list.files(full_pred_path, pattern="\\.rds$", full.names=TRUE)
  
  if (length(pred_file) == 0) {
    warning(paste("No .rds file found for", model_name, "in", full_pred_path))
    next  # Skip to the next iteration if no file is found
  }
  
  pred <- readRDS(pred_file)
  pred$template<-terra::unwrap(pred$template)

  # Predict
  map_i <- nsdm.map(template=pred$template,
                     nona_ix=pred$nona_ix,
                     species_name=ispi_name, model_name=model_name, nesting_name=nesting_method, level="reg",
                     pred=pred$ndata_bck) 
  
  # Save
  save_path <- file.path(scr_path, "outputs", "d7_maps/reg", nesting_method)
  nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="tif", save_path=save_path)
  
  cat(paste0(model_name, " predictions saved\n"))
}

### =========================================================================
### D- Ensemble predictions
### =========================================================================

ensemble_reg <- nsdm.ensemble(model_names = mod_algo,
                              species_name = ispi_name,
                              level = "reg", nesting_name=nesting_method,
                              map_path = file.path(scr_path, "outputs", "d7_maps/reg", nesting_method), 
                              score_path = file.path(scr_path, "outputs", "d3_evals/reg", nesting_method), 
                              weighting = as.logical(do_weighting), 
                              weight_metric = weight_metric, 
                              discthre = disc_thre)

nsdm.savemap(maps = ensemble_reg$ensemble, species_name = ispi_name, model_name = NULL, format="tif", 
               save_path = file.path(scr_path, "outputs", "d8_ensembles/reg", nesting_method))

nsdm.savemap(maps = ensemble_reg$ensemble_cv, species_name = ispi_name, model_name = NULL, format="tif", 
               save_path = file.path(scr_path, "outputs", "d9_ensembles-cv/reg", nesting_method))

cat("Ensemble predictions saved \n")
}

### =========================================================================
### D.2- Ensemble predictions evaluation
### =========================================================================
# Load test sets
d0_test_train <- readRDS(file.path(scr_path, "outputs", "d0_datasets", "base", ispi_name, paste0(ispi_name, ".rds")))$all_sets

scores_array<-nsdm.ensembleeval(sets = d0_test_train, level = "reg", model_names = mod_algo, species_name = ispi_name, scratch_path=scr_path, nesting_name=nesting_methods)

print(lapply(scores_array$scores_array, round, 2))

nsdm.savethis(object = scores_array,
              species_name = ispi_name,
              compression = TRUE,
              save_path = file.path(scr_path, "outputs", "d11_evals-final", "reg"))

cat("Ensemble predictions evaluated \n")

### =========================================================================
### E- Combine REG and GLO predictions
### =========================================================================
for(nesting_method in nesting_methods){
# E.1.1.1 "Multiply" (geometric mean) nesting
if (nesting_method == "multiply") {
  ## Ensembling
  # response
  ensemble_glo <- rast(list.files(
    file.path(scr_path, "outputs", "d8_ensembles", "glo", ispi_name),
    pattern = ".tif", full.names = TRUE
  ))
  
  ensemble_reg <- rast(list.files(
    file.path(scr_path, "outputs", "d8_ensembles", "reg", nesting_method, ispi_name),
    pattern = ".tif", full.names = TRUE
  ))
  
  ensemble_nested <- sqrt(ensemble_glo * ensemble_reg)
  names(ensemble_nested) <- names(ensemble_reg)
   
  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", nesting_method))


# E.1.1.2 "Multiply" (weighted geometric mean) nesting
if (multiply_weighted == TRUE) {
  # Define weights
  w_glo <- scores_array$w_glo
  w_reg <- scores_array$w_reg
  
  # Weighted geometric mean: (glo^w1 * reg^w2)^(1 / (w1 + w2))
  weighted_product <- (ensemble_glo ^ w_glo) * (ensemble_reg ^ w_reg)
  ensemble_nested <- weighted_product ^ (1 / (w_glo + w_reg))

  # Match names
  names(ensemble_nested) <- gsub("_multiply_", "_multiplyw_", names(ensemble_reg))

  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", paste0(nesting_method,"w")))
}}

# E.1.2 "Covariate" nesting
if (nesting_method == "covariate") {
  ## Ensembling
  ensemble_nested <- rast(list.files(
    file.path(scr_path, "outputs", "d8_ensembles", "reg", nesting_method, ispi_name),
    pattern = ".tif", full.names = TRUE
  ))
  
  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", "covariate"))
}
}

cat("GLO and REG predictions nested \n")

cat("Finished!\n")
