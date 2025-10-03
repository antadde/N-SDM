#############################################################################
## Script: 3_mainSCE_B
## Purpose: Ensembling and mapping
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
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
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# SBATCH param
args <- commandArgs(trailingOnly = TRUE)
arrayID <- as.numeric(args[1])

# Target species
species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

ispi_name <- species[arrayID]

####

for (scenar in proj_scenarios) {
  for (per in proj_periods) {
    
cat(paste0("Starting mapping and ensembling of ", scenar, " ", per, 
           " scenario GLO projections for ", ispi_name, "...\n"))
    
    ### =========================================================================
    ### C- Save prediction raster
    ### =========================================================================
    
    for (i in seq_along(mod_algo)) {
      # C.1 Load raw prediction data
      model_name <- mod_algo[i]
      pred_path <- file.path(scr_path, "outputs", "d12_preds-sce", "glo", scenar, per)
      full_pred_path <- file.path(pred_path, ispi_name, model_name)
      
      pred_file <- list.files(full_pred_path, pattern = ".rds", full.names = TRUE)
      pred <- readRDS(pred_file)
      
      # C.2 Predict
      map_i <- nsdm.map(
        template = unwrap(pred$template),
        nona_ix = pred$nona_ix,
        species_name = ispi_name,
        model_name = model_name,
        level = "glo",
        scenar_name = scenar,
        period_name = per,
        pred = pred$ndata_bck
      ) 
      
      # C.3 Save
      nsdm.savemap(
        maps = map_i, 
        species_name = ispi_name, 
        model_name = model_name, 
        format = "tif", 
        save_path = file.path(scr_path, "outputs", "d13_maps-sce", "glo", scenar, per)
      )
      
      cat(paste0(model_name, " projections saved \n"))
    }
    
    ### =========================================================================
    ### D- Ensemble projections
    ### =========================================================================
    
    ensemble_glo <- nsdm.ensemble(
      model_names = mod_algo,
      species_name = ispi_name,
      level = "glo",
      scenar_name = scenar,
      period_name = per,
      map_path = file.path(scr_path, "outputs", "d13_maps-sce", "glo", scenar, per), 
      score_path = file.path(scr_path, "outputs", "d3_evals", "glo"),
      weighting = as.logical(do_weighting),
      weight_metric = weight_metric, 
      discthre = disc_thre_glo
    )
    
    nsdm.savemap(
      maps = ensemble_glo$ensemble, 
      species_name = ispi_name, 
      model_name = NULL, format = "tif",
      save_path = file.path(scr_path, "outputs", "d14_ensembles-sce", "glo", scenar, per)
    )
    
    nsdm.savemap(
      maps = ensemble_glo$ensemble_cv, 
      species_name = ispi_name, 
      model_name = NULL, format = "tif",
      save_path = file.path(scr_path, "outputs", "d15_ensembles-cv-sce", "glo", scenar, per)
    )
  }
}

cat("Ensemble projections saved \n")
cat("Finished!\n")
