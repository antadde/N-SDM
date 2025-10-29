#############################################################################
## Script: 3_mainSCE_D
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

# SBATCH array
array<-expand.grid(nesting=nesting_methods, species=species)
ispi_name <- as.character(array[arrayID,"species"])
nesting_method <- array[arrayID,"nesting"]

for (scenar in proj_scenarios) {
  for (per in proj_periods) {
    
cat(paste0("Starting mapping and ensembling of ", scenar, " ", per, 
           " scenario REG projections for ", ispi_name, 
           " under the ", nesting_method, " nesting method(s)...\n"))
    
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
    
	if (!is.null(ensemble_reg$ensemble_cv)) {
    nsdm.savemap(
      maps = ensemble_reg$ensemble_cv, 
      species_name = ispi_name, 
      model_name = NULL, 
      format = "tif", 
      save_path = file.path(scr_path, "outputs", "d15_ensembles-cv-sce", "reg", nesting_method, scenar, per)
    )}
    
    ### =========================================================================
    ### E- Combine REG and GLO projections
    ### =========================================================================
    
    # E.1 "Multiply" nesting
	if (nesting_method == "posthoc") {
	
  # Load ensembles
  ## Global
      ensemble_glo <- rast(list.files(
        file.path(scr_path, "outputs", "d14_ensembles-sce", "glo", scenar, per, ispi_name), 
        pattern = ".tif", full.names = TRUE
      ))
	  
  ## Regional
	ensemble_reg <- ensemble_reg$ensemble  
	
	if (any(c("multiply") %in% posthoc_nesting_methods)){

     ensemble_nested <- sqrt(ensemble_glo * ensemble_reg)
     
	 # Rename
     names(ensemble_nested) <- gsub("_posthoc_", "_multiply_", names(ensemble_reg))
      
     # Save
      nsdm.savemap(
        map = ensemble_nested, 
        species_name = ispi_name, format="tif",
        save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", "multiply", scenar, per)
      )
    
	
	# E.1.2 "Multiply" (weighted geometric mean) nesting
	if (any(c("multiplyw") %in% posthoc_nesting_methods)){

	
	scores_array <- nsdm.loadthis(
        species_name = ispi_name,
        read_path = file.path(scr_path, "outputs", "d11_evals-ensembles", "reg"),
		format = "psv"
      )
	  
	  # Define weights
  w_glo <- mean(scores_array[scores_array$Level == "GLO" & scores_array$Metric == weight_metric, ]$Value)
  w_reg <- mean(scores_array[scores_array$Level == "REG" & scores_array$Metric == weight_metric, ]$Value)
	  
  # Weighted geometric mean
  weighted_product <- (ensemble_glo ^ w_glo) * (ensemble_reg ^ w_reg)
  ensemble_nested <- weighted_product ^ (1 / (w_glo + w_reg))

  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_multiplyw_", names(ensemble_reg))
	  
	  # Save
	  nsdm.savemap(
		map = ensemble_nested, 
		species_name = ispi_name, format="tif",
		save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", "multiplyw", scenar, per)
	  )
	}
	
# E.1.3 "Average" (arithmetic mean) nesting
if (any(c("average") %in% posthoc_nesting_methods)){
	
  ensemble_nested <- mean(ensemble_glo, ensemble_reg)
  
  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_average_", names(ensemble_reg))
   
  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", "average", scenar, per))
}	

# E.1.4 "Average weighted" (weighted arithmetic mean) nesting
if (any(c("averagew") %in% posthoc_nesting_methods)){
  # Define weights
  w_glo <- mean(scores_array[scores_array$Level == "GLO" & scores_array$Metric == weight_metric, ]$Value)
  w_reg <- mean(scores_array[scores_array$Level == "REG" & scores_array$Metric == weight_metric, ]$Value)
  
  # Weighted arithmetic mean
  ensemble_nested <- (w_glo * ensemble_glo + w_reg * ensemble_reg) / (w_glo + w_reg)


  # Rename
  names(ensemble_nested) <- gsub("_posthoc_", "_averagew_", names(ensemble_reg))

  # Save
  nsdm.savemap(map = ensemble_nested, species_name = ispi_name, format="tif",
               save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", "averagew", scenar, per))
}
}}
	    
    # E.2 "Covariate" nesting
    if (nesting_method == "covariate") {
      
	  ensemble_nested <- ensemble_reg$ensemble
      
      # Save
      nsdm.savemap(
        map = ensemble_nested, 
        species_name = ispi_name, format="tif",
        save_path = file.path(scr_path, "outputs", "d16_nested-ensembles-sce", "covariate", scenar, per)
      )
	}
  }
}

cat("GLO and REG projections nested and saved \n")
cat("Finished!\n")
