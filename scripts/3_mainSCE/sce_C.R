#############################################################################
## Script: 3_mainSCE_1
## Purpose: REG projections
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
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# SBATCH param
args <- commandArgs(trailingOnly = TRUE)
arrayID <- as.numeric(args[1])

# Target species
species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

# SBATCH array
array <- expand.grid(nesting = nesting_methods, species = species, scenarios = proj_scenarios)
ispi_name <- as.character(array[arrayID, "species"])
nesting_method <- array[arrayID, "nesting"]
scenar <- array[arrayID, "scenarios"]

  for (model_name in mod_algo) {
  for (per in proj_periods) {
    
cat(paste0("Starting computation of ", scenar, " ", per, " ", toupper(model_name),  
           " REG projections for ", ispi_name, 
           " using the ", nesting_method, " method...\n"))
		   
    ### =========================================================================
    ### C- Load REG model
    ### =========================================================================
    
    # Load REG data
    d0_datasets <- nsdm.loadthis(
      species_name = ispi_name,
      read_path = file.path(scr_path, "outputs", "d0_datasets", "reg")
    )
     
    # Specific loading strategy for lgb.booster
    if (model_name == "gbm") {
      prmod2 <- lgb.load(
        file.path(scr_path, "outputs", "d2_models", "reg", 
                  nesting_method, ispi_name, "gbm", paste0(ispi_name, "_gbm.rds"))
      )
      prmod <- nsdm.loadthis(
        model_name = "glm",
        species_name = ispi_name,
        read_path = file.path(scr_path, "outputs", "d2_models", "reg", nesting_method)
      )
      prmod@fits[[1]] <- prmod2
    } else {
	    # Load REG model
    prmod <- nsdm.loadthis(
      model_name = model_name,
      species_name = ispi_name,
      read_path = file.path(scr_path, "outputs", "d2_models", "reg", nesting_method)
    )}
	
    # List covariates
    cov <- unlist(strsplit(prmod@meta$env_vars, ", "))
	if ("mainGLO" %in% cov) cov <- c(setdiff(cov, "mainGLO"), "mainGLO")
    
### =========================================================================
### D - Load scenario layers
### =========================================================================

# Retrieve list of candidate covariates and covinfo table
lr_file <- file.path(w_path, "tmp", "settings", "ref_covariates.rds")
lr <- readRDS(lr_file)
cov_info <- lr$cov_info
cov_info$ID <- paste(cov_info$cada, cov_info$variable, cov_info$attribute, 
                     cov_info$focal, sep = "_")
    
# Align IDs with cov order
cov_match <- gsub(".tif", "", basename(cov_info$file))
cov_ID <- cov_info$ID[match(cov, cov_match)]

# Presence covariates
if (n_levels == 1) {
  cov_info_pres <- cov_info[cov_info$ID %in% cov_ID & 
                            (is.na(cov_info$scenario) | trimws(cov_info$scenario) == ""), ]
}

if (n_levels == 2) {
  cov_info_pres <- cov_info[cov_info$ID %in% cov_ID & 
                            cov_info$level == "reg" & 
                            (is.na(cov_info$scenario) | trimws(cov_info$scenario) == ""), ]
}

# Deduplicate (keep last occurrence) and reorder to cov_ID
cov_info_pres <- cov_info_pres[!duplicated(cov_info_pres$ID, fromLast = TRUE), ]
cov_info_pres <- cov_info_pres[match(cov_ID, cov_info_pres$ID), ]
lr_pres_ID <- cov_info_pres$ID

# --- D.1 List available scenario layers for cov_ID ---
cov_info_sce <- cov_info[
  cov_info$ID %in% cov_ID & 
  cov_info$scenario == scenar & 
  cov_info$year == per, ]
    
# Deduplicate and reorder to cov_ID
cov_info_sce <- cov_info_sce[!duplicated(cov_info_sce$ID, fromLast = TRUE), ]
cov_info_sce <- cov_info_sce[match(cov_ID, cov_info_sce$ID), ]

lr_sce <- cov_info_sce$file
lr_sce_ID <- cov_info_sce$ID

# --- Skip alignment and loading if lr_sce is all NA ---
if (length(lr_sce) > 0 && !all(is.na(lr_sce))) {
  
  # Align names with pres files
  lr_sce_names <- gsub(".tif", "", basename(cov_info_pres$file))[match(lr_sce_ID, lr_pres_ID)]
  
  # D.2 Load target scenario layers
  stk_sce <- lapply(lr_sce, function(f) toMemory(rast(f)))
  stk_sce <- rast(stk_sce)
  
} else {
  lr_sce_names <- character(0)  # or NA, depending on what you want downstream
  stk_sce <- list()
}

 # --- D.3 Complete with non-scenario layers if needed ---
remainders <- setdiff(cov_ID, lr_sce_ID)

if (length(remainders) > 0) {
  cov_info_remain <- cov_info_pres[!duplicated(cov_info_pres$ID, fromLast = TRUE), ]
  cov_info_remain <- cov_info_remain[match(remainders, cov_info_remain$ID), ]
  lr_remain <- cov_info_remain$file
  
  stk_remain <- lapply(lr_remain, function(f) toMemory(rast(f)))
  stk_remain <- rast(stk_remain)

  # Combine scenario and remainder layers, preserving cov_ID order
  stk_sce <- rast(c(stk_sce, stk_remain))
}

# --- D.4 Load scenario GLO output if needed ---
if ("mainGLO" %in% cov) {
  glo_out <- list.files(
    file.path(scr_path, "outputs", "d14_ensembles-sce", "glo", scenar, per, ispi_name), 
    pattern = ".tif", full.names = TRUE
  )
  glo <- toMemory(rast(glo_out))
  names(glo) <- "mainGLO"
  stk_sce <- c(stk_sce, glo)
  }

# --- Final rename ---
names(stk_sce) <- cov
    
    ### =========================================================================
    ### E- Spatial projections
    ### =========================================================================
    
    ## E.1 Prepare covariate data for projections
    if (length(cov_observ) > 0) {
      cov_obs <- grep(paste0(cov_observ, collapse = "|"), names(stk_sce), value = TRUE)
    } else {
      cov_obs <- NULL
    }
    
    hab_df_reg <- nsdm.retrieve4pred(
      covstk = stk_sce,
      observational = cov_obs,
      obsval = cov_observ_val,
      mask = mask_pred,
      scaleparam = attributes(d0_datasets$env_vars)[c("scaled:center", "scaled:scale")]
    )
    
    ## E.2 Clean workspace to free some memory before predicting
	template <- stk_sce[[1]]
	suppressWarnings(suppressMessages({
	  rm(d0_datasets, stk_sce)
	  invisible(gc())  # Suppress the output of garbage collection
	}))
		
    ## E.3 Predict
    ndata_bck <- nsdm.predict(
      models = prmod,
      nwdata = hab_df_reg$covdf,
      nsplits = ncores
    )
    
    ## E.4 Save
	template <- wrap(template)

    nsdm.savethis(
      object = list(
        ndata_bck = ndata_bck,
        template = template,
        nona_ix = hab_df_reg$covdf_ix
      ),
      model_name = model_name,
      species_name = ispi_name,
      save_path = file.path(scr_path, "outputs", "d12_preds-sce", "reg", nesting_method, scenar, per)
    )
  }
}

cat("Predictions calculated and saved \n")
cat("Finished!\n")
