#############################################################################
## Script: 1_mainGLO_C
## Author: Antoine Adde
#############################################################################

### =========================================================================
### A. Preparation
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

ispi_name<-species[arrayID]

cat(paste0("Starting mapping and ensembling of GLO predictions for ", ispi_name, "...\n"))

### =========================================================================
### C- Save prediction raster
### =========================================================================

for(i in seq_along(mod_algo)){
  # Load raw prediction data
  model_name <- mod_algo[i]
  pred_path <- file.path(scr_path, "outputs", "d6_preds/glo")
  full_pred_path <- file.path(pred_path, ispi_name, model_name)
  pred_file <- list.files(full_pred_path, pattern="\\.rds$", full.names=TRUE)
  
  if (length(pred_file) == 0) {
    warning(paste("No .rds file found for", model_name, "in", full_pred_path))
    next  
  }
  
  pred <- readRDS(pred_file)
  pred$template<-terra::unwrap(pred$template)

  # Predict
  map_i <- nsdm.map(template=pred$template,
                     nona_ix=pred$nona_ix,
                     species_name=ispi_name, model_name=model_name, level="glo",
                     pred=pred$ndata_bck) 
  
  # Save
  save_path <- file.path(scr_path, "outputs", "d7_maps/glo")
  nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="tif", save_path=save_path)
  
  cat(paste0(model_name, " predictions saved\n"))
}

### =========================================================================
### D- Ensemble predictions
### =========================================================================

ensemble_glo <- nsdm.ensemble(model_names = mod_algo,
                              species_name = ispi_name,
                              level = "glo",
                              map_path = file.path(scr_path, "outputs", "d7_maps/glo"), 
                              score_path = file.path(scr_path, "outputs", "d3_evals/glo"), 
                              weighting = as.logical(do_weighting), 
                              weight_metric = weight_metric, 
                              discthre = disc_thre_glo)

nsdm.savemap(maps = ensemble_glo$ensemble, species_name = ispi_name, model_name = NULL, format="tif",
               save_path = file.path(scr_path, "outputs", "d8_ensembles/glo"))

nsdm.savemap(maps = ensemble_glo$ensemble_cv, species_name = ispi_name, model_name = NULL, format="tif",
               save_path = file.path(scr_path, "outputs", "d9_ensembles-cv/glo"))

cat("Ensemble predictions saved \n")

### =========================================================================
### D- Ensemble predictions evaluation
### =========================================================================
# Load selected covariates
d1_covsels <- readRDS(file.path(scr_path, "outputs", "d1_covsels", "glo", ispi_name, paste0(ispi_name, ".rds")))

# Load test sets
d0_test_train <- readRDS(file.path(scr_path, "outputs", "d0_datasets", "base", ispi_name, paste0(ispi_name, ".rds")))$all_sets

# Evaluate
scores_array<-nsdm.ensembleeval(sets = d0_test_train, level = "glo", model_names = mod_algo, species_name = ispi_name, scratch_path = scr_path)

scores_array_df = do.call(rbind, lapply(names(scores_array$scores_array), function(level)
  data.frame(
    Level = level,
    Metric = names(scores_array$scores_array[[level]]),
    Value = as.numeric(scores_array$scores_array[[level]]),
    row.names = NULL
  )))

print(lapply(scores_array$scores_array, round, 2))

### =========================================================================
### Save
### =========================================================================

nsdm.savethis(object = scores_array_df,
              species_name = ispi_name,
              compression = TRUE,
              save_path = file.path(scr_path, "outputs", "d11_evals-ensembles", "glo"),
			  format = "psv")

cat("Ensemble predictions evaluated \n")

cat("Finished!\n")
