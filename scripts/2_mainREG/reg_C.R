#############################################################################
## Script: 2_mainREG_C
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================

# Load N-SDM settings
load(file.path(gsub("scripts", "tmp", getwd()), "settings", "ref_nsdm_settings.RData"))

# global reproducibility seed
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

# Set permissions for new files
Sys.umask(mode = "000")

# Set working directory
setwd(w_path)

# Set library path
.libPaths(Rlib_path)

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
nesting_methods <- nesting_methods

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
                              discthre = disc_thre_reg)

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

scores_array<-nsdm.ensembleeval(sets = d0_test_train, level = "reg", model_names = mod_algo, species_name = ispi_name, scratch_path=scr_path, nesting_name=nesting_methods, posthoc_nesting_name=posthoc_nesting_methods)

print(lapply(scores_array$scores_array, round, 2))

scores_array_df = do.call(rbind, lapply(names(scores_array$scores_array), function(level)
  data.frame(
    Level = level,
    Metric = names(scores_array$scores_array[[level]]),
    Value = as.numeric(scores_array$scores_array[[level]]),
    row.names = NULL
  )))

nsdm.savethis(object = scores_array_df,
              species_name = ispi_name,
              compression = TRUE,
              save_path = file.path(scr_path, "outputs", "d11_evals-ensembles", "reg"),
			  format = "psv")

cat("Ensemble predictions evaluated \n")

### =========================================================================
### E- Combine REG and GLO predictions
### =========================================================================
for(nesting_method in nesting_methods){
# E.1 "Posthoc" nesting
if (nesting_method == "posthoc") {
	
  # Load ensembles
  ## Global
  ensemble_glo <- rast(list.files(
    file.path(scr_path, "outputs", "d8_ensembles", "glo", ispi_name),
    pattern = ".tif", full.names = TRUE
  ))
  
  ## Regional
  ensemble_reg <- rast(list.files(
    file.path(scr_path, "outputs", "d8_ensembles", "reg", nesting_method, ispi_name),
    pattern = ".tif", full.names = TRUE
  ))
  
# E.1.1 "Multiply" (geometric mean) nesting
if (any(c("multiply") %in% posthoc_nesting_methods)){
	
  ensemble_nested <- sqrt(ensemble_glo * ensemble_reg)
  
  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_multiply_", names(ensemble_reg))
   
  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", "multiply"))
}

# E.1.2 "Multiply weighted" (weighted geometric mean) nesting
if (any(c("multiplyw") %in% posthoc_nesting_methods)){
	
  # Define weights
  w_glo <- scores_array_df[scores_array_df$Level == "GLO" & scores_array_df$Metric == weight_metric, ]$Value
  w_reg <- scores_array_df[scores_array_df$Level == "REG" & scores_array_df$Metric == weight_metric, ]$Value
  
  # Weighted geometric mean
  weighted_product <- (ensemble_glo ^ w_glo) * (ensemble_reg ^ w_reg)
  ensemble_nested <- weighted_product ^ (1 / (w_glo + w_reg))

  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_multiplyw_", names(ensemble_reg))

  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", "multiplyw"))
}

# E.1.3 "Average" (arithmetic mean) nesting
if (any(c("average") %in% posthoc_nesting_methods)){
	
  ensemble_nested <- mean(ensemble_glo, ensemble_reg)
  
  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_average_", names(ensemble_reg))
   
  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", "average"))
}

# E.1.4 "Average weighted" (weighted arithmetic mean) nesting
if (any(c("averagew") %in% posthoc_nesting_methods)){
  # Define weights
  w_glo <- scores_array_df[scores_array_df$Level == "GLO" & scores_array_df$Metric == weight_metric, ]$Value
  w_reg <- scores_array_df[scores_array_df$Level == "REG" & scores_array_df$Metric == weight_metric, ]$Value
  
  # Weighted arithmetic mean
  ensemble_nested <- (w_glo * ensemble_glo + w_reg * ensemble_reg) / (w_glo + w_reg)


  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_averagew_", names(ensemble_reg))

  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d10_nested-ensembles", "averagew"))
}
}

# E.2 "Covariate" nesting
if (nesting_method == "covariate") {
  ## Load ensemble
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
