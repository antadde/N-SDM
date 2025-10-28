#############################################################################
## Script: 1_mainGLO_B
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
### B. Definitions
### =========================================================================

# Number of cores to be used during parallel operations
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# SBATCH param
args <- commandArgs(trailingOnly = TRUE)
arrayID <- as.numeric(args[1])

# Target species
species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

# Target model algorithms
models<-mod_algo

# SBATCH array
array<-expand.grid(model=models, species=species)
ispi_name <- array[arrayID,"species"]
model_name <- array[arrayID,"model"]

cat(paste0("Starting global-level ", toupper(model_name), " modelling for ", ispi_name, "...\n"))

### =========================================================================
### C. Load glo_A Outputs
### =========================================================================

# Load datasets
d0_datasets <- readRDS(file.path(scr_path, "outputs", "d0_datasets", "glo", ispi_name, paste0(ispi_name, ".rds")))

# Load training and testing sets
d0_test_train <- readRDS(file.path(scr_path, "outputs", "d0_datasets", "base", ispi_name, paste0(ispi_name, ".rds")))$all_sets

# Load selected covariates
d1_covsels <- readRDS(file.path(scr_path, "outputs", "d1_covsels", "glo", ispi_name, paste0(ispi_name, ".rds")))

cat("glo_A outputs loaded\n")

### =========================================================================
### D. Define model parameters 
### =========================================================================
modinp <- nsdm.setparam(
  model_name     = model_name,
  covariate_names = names(d1_covsels$pseu.abs_i@env_vars),
  param_grid     = param_grid,
  tmp_path       = file.path(scr_path, "tmp"),
  weights        = d0_datasets$weights,
  ncov.esm       = ncov_esm,
  comb.esm       = comb_esm,
  nthreads       = ncores       
)

cat(paste0('Modelling parameters defined \n'))

### =========================================================================
### E. Fit, Evaluate, and Save Models
### =========================================================================
## E.1 Fit Model
mod <- try(
  nsdm.fitting(
    x = d1_covsels$pseu.abs_i,
	sets = d0_test_train,
    mod_args = modinp,
    ncores = ncores,
    level = "glo",
    tmp_path = file.path(scr_path, "tmp")
  ), 
  silent = TRUE
)

## E.2 Evaluate Model
if (!inherits(mod, "try-error")) {
  
  ### Assessment metrics
  evals <- try(
    nsdm.eval(
	  x = d1_covsels$pseu.abs_i,
      models = mod,
	  sets = d0_test_train,
      ncores = ncores, 
      level = "glo",
      tmp_path = file.path(scr_path, "tmp")
    ), 
    silent = TRUE
  )
  
  if (!inherits(evals, "try-error")) {  
    ### Bind results
    smev <- nsdm.summary(evals)    
    print(smev)
  }
}

## E.3 Save Model
if (!inherits(mod, "try-error")) {
  suppressWarnings(
  nsdm.savethis(
      object = list(model = mod, parameters = modinp),
      species_name = ispi_name, 
      model_name = model_name,
      tag = paste(model_name, "tune", sep = "_"),
      save_path = file.path(scr_path, "outputs", "d2_models", "glo")
    )
  )
}

## E.4 Save Evaluation Table
if (exists("smev")) {
  suppressWarnings(
    nsdm.savethis(
      object =  data.frame(Metric = rownames(smev), smev, row.names = NULL, check.names = FALSE),
      model_name = model_name, 
      species_name = ispi_name,
      compression = TRUE,
      save_path = file.path(scr_path, "outputs", "d3_evals", "glo"),
	  format = "psv"
    )
  )
}

cat("\n\nModels fitted and evaluated\n")


### =========================================================================
### F - Refit Top Model for Prediction
### =========================================================================

## F.1.1 Identify Best Model
if (model_name != "esm") {
  if (ncol(smev) > 1) {
    ord <- sort(smev[best_met, ], decreasing = TRUE)
    modinp_top <- modinp[names(ord[1])]
  } else {
    modinp_top <- modinp
  }
}

## F.1.2 ... or Discard Models with best_thre_esm < (for esm)
if (model_name == "esm") {
  ord <- sort(smev[best_met, ], decreasing = TRUE)
  ord <- ord[ord > best_thre_esm]
  modinp_top <- modinp[names(ord)]
}

## F.2a Refit Model(s) Using Full Dataset
suppressWarnings(
  prmod <- nsdm.fitfull(
    x = d1_covsels$pseu.abs_i,
    mod_args = modinp_top)
)

## Save Model
suppressWarnings(
  nsdm.savethis(
    object = prmod,
    species_name = ispi_name,
    model_name = model_name,
    tag = model_name,
    compression = TRUE,
    save_path = file.path(scr_path, "outputs", "d2_models", "glo")
  )
)

## Save GBM Model Separately
if (model_name == "gbm") {
  gbm_path <- file.path(scr_path, "outputs", "d2_models", "glo", ispi_name, "gbm")
  dir.create(gbm_path, recursive = TRUE, showWarnings = FALSE)
  lgb.save(
    prmod@fits[[1]], 
    file.path(gbm_path, paste0(ispi_name, "_", model_name, ".rds"))
  )
}

cat("Top model ", names(modinp_top), " refitted on full dataset for predictions \n")


### =========================================================================
### G - Compute Covariate Importance and Response Curves
### =========================================================================

## G.1 Covariate Importance
imp <- nsdm.varimp(prmod)
print(imp)

if (model_name == "esm") {
imp <- do.call(rbind, lapply(names(imp), function(esm_name) {
     data.frame(
         esm = esm_name,
         imp[[esm_name]],
         stringsAsFactors = FALSE
     )
 })) }

nsdm.savethis(
  object = imp,
  model_name = model_name, 
  species_name = ispi_name,
  compression = TRUE,
  save_path = file.path(scr_path, "outputs", "d4_covimps", "glo"),
  format = "psv"
)

## G.2 Response Curves
Data <- d1_covsels$pseu.abs_i@env_vars

respcurves <- nsdm.respcurve(
  prmod,
  Data = Data,
  scaleparam = attributes(d0_datasets$env_vars)[c("scaled:center", "scaled:scale")],
  model_name = model_name, 
  species_name = ispi_name,
  plotting = TRUE, 
  ncores = ncores, 
  save_path = file.path(scr_path, "outputs", "plots", "respcurves", "glo")
)

nsdm.savethis(
  object = respcurves,
  model_name = model_name, 
  species_name = ispi_name,
  compression = TRUE,
  save_path = file.path(scr_path, "outputs", "d5_respcurves", "glo")
)

cat("\n\nVariable importance scores and response curves computed\n")

### =========================================================================
### H. Spatial Predictions
### =========================================================================
## H.0 Prepare Covariate Data for Predictions (optional glo extent projections)
if(n_levels == 2 && glo_full_extent == TRUE){
covstk_glo<-unwrap(d1_covsels$covstk_glo)

if (length(cov_observ) > 0) {
  cov_obs_glo <- grep(paste(cov_observ, collapse = "|"), names(covstk_glo), value = TRUE)
} else {
  cov_obs_glo <- NULL
}

stk_df_glo <- nsdm.retrieve4pred(
  covstk = covstk_glo, 
  observational = cov_obs_glo,
  obsval = cov_observ_val,
  mask = mask_pred_glo,
  scaleparam = attributes(d0_datasets$env_vars)[c("scaled:center", "scaled:scale")]
)

template_glo <- covstk_glo[[1]]
template_glo <- wrap(template_glo)
}


## H.1 Prepare Covariate Data for Predictions
covstk_reg<-unwrap(d1_covsels$covstk)

if (length(cov_observ) > 0) {
  cov_obs_reg <- grep(paste(cov_observ, collapse = "|"), names(covstk_reg), value = TRUE)
} else {
  cov_obs_reg <- NULL
}

stk_df_reg <- nsdm.retrieve4pred(
  covstk = covstk_reg, 
  observational = cov_obs_reg,
  obsval = cov_observ_val,
  mask = mask_pred_reg,
  scaleparam = attributes(d0_datasets$env_vars)[c("scaled:center", "scaled:scale")]
)

template_reg <- covstk_reg[[1]]
template_reg <- wrap(template_reg)

## H.2 Clean Workspace to Free Memory Before Predicting
suppressWarnings(suppressMessages({
  rm(d0_datasets, d1_covsels, covstk_glo, covstk_reg, respcurves, imp, eval_list)
  invisible(gc())
}))

## H.3 Predict
if(n_levels == 2 && glo_full_extent == TRUE){
ndata_bck_glo <- nsdm.predict(
  models = prmod,
  nwdata = stk_df_glo$covdf,
  nsplits = ncores
)

nsdm.savethis(
  object = list(ndata_bck = ndata_bck_glo, template = template_glo, nona_ix = stk_df_glo$covdf_ix),
  model_name = model_name, 
  species_name = ispi_name,
  save_path = file.path(scr_path, "outputs", "d6_preds", "glo_full_extent")
)
}

ndata_bck_reg <- nsdm.predict(
  models = prmod,
  nwdata = stk_df_reg$covdf,
  nsplits = ncores
)

nsdm.savethis(
  object = list(ndata_bck = ndata_bck_reg, template = template_reg, nona_ix = stk_df_reg$covdf_ix),
  model_name = model_name, 
  species_name = ispi_name,
  save_path = file.path(scr_path, "outputs", "d6_preds", "glo")
)

cat("Predictions calculated and saved\n")
cat("Finished!\n")
