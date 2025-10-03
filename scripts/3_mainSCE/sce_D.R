#############################################################################
## Script: 3_mainSCE_D
## Purpose: Ensembling and mapping
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A. Preparation
### =========================================================================

# global reproducibility seed
set.seed(123)

# Load N-SDM settings
load(file.path(gsub("scripts", "tmp", gsub("/3_mainSCE", "", getwd())), "settings", "tmp_nsdm_settings.RData"))

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

# SBATCH array
array<-expand.grid(nesting=nesting_methods, species=species)
ispi_name <- as.character(array[arrayID,"species"])
nesting_method <- array[arrayID,"nesting"]

for (scenar in proj_scenarios) {
  for (per in proj_periods) {
    
cat(paste0("Starting mapping and ensembling of ", scenar, " ", per, 
           " scenario REG projections for ", ispi_name, 
           " under the ", nesting_method, " nesting method...\n"))
    
    ### =========================================================================
    ### C- Save prediction raster
    ### =========================================================================
    
    for (i in seq_along(mod_algo)) {
      # Load raw prediction data
      model_name <- mod_algo[i]
      pred_path <- file.path(scr_path, "outputs", "d12_preds-sce", "reg", nesting_method, scenar, per)
      full_pred_path <- file.path(pred_path, ispi_name, model_name)
      
      pred_file <- list.files(full_pred_path, pattern = ".rds", full.names = TRUE)
      pred <- readRDS(pred_file)
      
      # Predict 
      map_i <- nsdm.map(
        template = unwrap(pred$template),
        nona_ix = pred$nona_ix, 
        species_name = ispi_name,
        model_name = model_name,
        level = "reg",
        scenar_name = scenar,
        period_name = per,
        nesting_name = nesting_method,
        pred = pred$ndata_bck
      ) 
      
      # Save
      nsdm.savemap(
        maps = map_i, 
        species_name = ispi_name, 
        model_name = model_name, 
        format = "tif", 
        save_path = file.path(scr_path, "outputs", "d13_maps-sce", "reg", nesting_method, scenar, per)
      )
      
      cat(paste0(model_name, " projections saved \n"))
    }
    
    ### =========================================================================
    ### D- Ensemble projections
    ### =========================================================================
    
    ensemble_reg <- nsdm.ensemble(
      model_names = mod_algo,
      species_name = ispi_name,
      level = "reg",
      scenar_name = scenar,
      period_name = per,
      nesting_name = nesting_method,
      map_path = file.path(scr_path, "outputs", "d13_maps-sce", "reg", nesting_method, scenar, per),
      score_path = file.path(scr_path, "outputs", "d3_evals", "reg", nesting_method),
      weighting = as.logical(do_weighting),
      weight_metric = weight_metric,
      discthre = disc_thre_reg
    )
    
    nsdm.savemap(
      maps = ensemble_reg$ensemble, 
      species_name = ispi_name, 
      model_name = NULL, format="tif",
      save_path = file.path(scr_path, "outputs", "d14_ensembles-sce", "reg", nesting_method, scenar, per)
    )
    
    nsdm.savemap(
      maps = ensemble_reg$ensemble_cv, 
      species_name = ispi_name, 
      model_name = NULL, 
      format = "tif", 
      save_path = file.path(scr_path, "outputs", "d15_ensembles-cv-sce", "reg", nesting_method, scenar, per)
    )
    
    ### =========================================================================
    ### E- Combine REG and GLO projections
    ### =========================================================================
    
    # E.1 "Multiply" nesting
    if (nesting_method == "multiply") {
      # Load response data
      ensemble_glo <- rast(list.files(
        file.path(scr_path, "outputs", "d14_ensembles-sce", "glo", scenar, per, ispi_name), 
        pattern = ".tif", full.names = TRUE
      ))
      
      ensemble_nested <- sqrt(ensemble_glo * ensemble_reg$ensemble)
      names(ensemble_nested) <- names(ensemble_reg$ensemble)
      
     # Save
      nsdm.savemap(
        map = ensemble_nested, 
        species_name = ispi_name, format="tif",
        save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", nesting_method, scenar, per)
      )
    
	
	# E.1.2 "Multiply" (weighted geometric mean) nesting
	if (multiply_weighted == TRUE) {
	
	scores_array <- nsdm.loadthis(
        species_name = ispi_name,
        read_path = file.path(scr_path, "outputs", "d11_evals-final", "reg")
      )
	  
	  # Define weights
	  w_glo <- scores_array$w_glo
	  w_reg <- scores_array$w_reg
	  
	  # Load response data
	  ensemble_glo <- rast(list.files(
		file.path(scr_path, "outputs", "d14_ensembles-sce", "glo", scenar, per, ispi_name), 
		pattern = ".tif", full.names = TRUE
	  ))

	  # Compute weighted geometric mean
	  weighted_product <- (ensemble_glo ^ w_glo) * (ensemble_reg$ensemble ^ w_reg)
	  ensemble_nested <- weighted_product ^ (1 / (w_glo + w_reg))
	  
	  # Match names
	  names(ensemble_nested) <- gsub("_multiply_", "_multiplyw_", names(ensemble_reg$ensemble))
	  
	  # Save
	  nsdm.savemap(
		map = ensemble_nested, 
		species_name = ispi_name, format="tif",
		save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", paste0(nesting_method,"w"), scenar, per)
	  )
	}}
    
    # E.2 "Covariate" nesting
    if (nesting_method == "covariate") {
      ensemble_nested <- ensemble_reg$ensemble
      
      # Save
      nsdm.savemap(
        map = ensemble_nested, 
        species_name = ispi_name, format="tif",
        save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", nesting_method, scenar, per)
      )
	}
  }
}

cat("GLO and REG projections nested and saved \n")
cat("Finished!\n")
