#############################################################################
## Script: 3_mainSCE_A
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

# Target scenario for projection
scenars<-proj_scenarios

# SBATCH array
array<-expand.grid(species=species, scenarios=scenars)
ispi_name <- as.character(array[arrayID,"species"])
scenar<-as.character(array[arrayID,"scenarios"])

  for (model_name in mod_algo) {
  for (per in proj_periods) {
cat(paste0("Starting computation of ", scenar, " ", per, " ", toupper(model_name),  
           " GLO projections for ", ispi_name, "...\n"))
    
    ### =========================================================================
    ### C- Load GLO model
    ### =========================================================================
    
    # Load GLO data
    d0_datasets <- nsdm.loadthis(
      species_name = ispi_name,
      read_path = file.path(scr_path, "outputs", "d0_datasets", "glo")
    )
    
    # Specific loading strategy for lgb.booster
    if (model_name == "gbm") {
	
	substitut <- c("glm", "gam", "max", "rf")[c("glm", "gam", "max", "rf") %in% mod_algo][1]

      prmod2 <- lgb.load(
        file.path(scr_path, "outputs", "d2_models", "glo", 
                  ispi_name, "gbm", paste0(ispi_name, "_gbm.rds"))
      )
      prmod <- nsdm.loadthis(
        model_name = substitut,
        species_name = ispi_name,
		tag = substitut,
        read_path = file.path(scr_path, "outputs", "d2_models", "glo")
      )
      prmod@fits[[1]] <- prmod2
    } else {
	
	    # Load GLO model
    prmod <- nsdm.loadthis(
      model_name = model_name,
      species_name = ispi_name,
	  tag = model_name,
      read_path = file.path(scr_path, "outputs", "d2_models", "glo"))
	}
    
    # List covariates
    cov <- unlist(strsplit(prmod@meta$env_vars, ", "))
    
    ### =========================================================================
    ### B. Load scenario layers
    ### =========================================================================
    # Retrieve list of candidate covariates and covinfo table
    lr_file <- file.path(w_path, "tmp", "settings", "ref_covariates.rds")
    lr <- readRDS(lr_file)
    cov_info <- lr$cov_info
    cov_info$ID <- paste(cov_info$cada, cov_info$variable, cov_info$attribute, 
                         cov_info$focal, sep = "_")
						 
# Build the matching key
cov_match <- gsub(".tif", "", basename(cov_info$file))

# IDs in the same order as cov
cov_ID <- cov_info$ID[match(cov, cov_match)]

# Filter cov_info_pres
if (n_levels == 1) {
  cov_info_pres <- cov_info[cov_info$ID %in% cov_ID & 
                            (is.na(cov_info$scenario) | trimws(cov_info$scenario) == ""), ]
}

if (n_levels == 2) {
  cov_info_pres <- cov_info[cov_info$ID %in% cov_ID & 
                            cov_info$level == "reg" & 
                            (is.na(cov_info$scenario) | trimws(cov_info$scenario) == ""), ]
}

# Reorder according to the cov order
cov_info_pres <- cov_info_pres[match(cov_ID, cov_info_pres$ID), ]
 
# --- B.1 List available scenario layers for cov_ID ---
cov_info_sce <- cov_info[
  cov_info$ID %in% cov_ID &
  cov_info$scenario == scenar &
  cov_info$year == per, ]

# Reorder according to cov_ID
cov_info_sce <- cov_info_sce[match(cov_ID, cov_info_sce$ID), ]

lr_sce <- cov_info_sce$file
lr_sce_ID <- cov_info_sce$ID

# --- B.2 Load target scenario layers if available ---
if (length(lr_sce) > 0) {	
  stk_sce <- lapply(lr_sce, function(f) toMemory(rast(f)))
  stk_sce <- rast(stk_sce)
} else {
  stk_sce <- list()
}

# --- B.3 Complete with non-scenario layers if needed ---
remainders <- setdiff(cov_ID, lr_sce_ID)

if (length(remainders) > 0) {
  cov_info_remain <- cov_info_pres[match(remainders, cov_info_pres$ID), ]
  lr_remain <- cov_info_remain$file
  stk_remain <- lapply(lr_remain, function(f) toMemory(rast(f)))
  stk_remain <- rast(stk_remain)

  # Combine scenario and remainder layers, preserving cov_ID order
  stk_sce <- rast(c(stk_sce, stk_remain))
}

# --- Final rename according to model covariates ---
names(stk_sce) <- cov
    
    ### =========================================================================
    ### C- Spatial projections
    ### =========================================================================
    
    ## C.1 Prepare covariate data for projections
    clim_df_reg <- nsdm.retrieve4pred(
      covstk = stk_sce,
      scaleparam = attributes(d0_datasets$env_vars)[c("scaled:center", "scaled:scale")]
    )
    
    ## C.2 Clean workspace to free some memory before predicting
	template <- stk_sce[[1]]
	suppressWarnings(suppressMessages({
	  rm(d0_datasets, stk_sce)
	  invisible(gc()) 
	}))
    
    ## C.3 Predict
    ndata_bck <- nsdm.predict(
      models = prmod,
      nwdata = clim_df_reg$covdf,
      nsplits = ncores
    )
    
    ## C.4 Save
	template <- wrap(template)
	
    nsdm.savethis(
      object = list(
        ndata_bck = ndata_bck,
        template = template,
        nona_ix = clim_df_reg$covdf_ix
      ),
      model_name = model_name,
      species_name = ispi_name,
      save_path = file.path(scr_path, "outputs", "d12_preds-sce", "glo", scenar, per)
    )
  }
}

cat("Projections calculated and saved \n")
cat("Finished!\n")