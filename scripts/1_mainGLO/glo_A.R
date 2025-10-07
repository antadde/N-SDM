#############################################################################
## Script: 1_mainGLO_A
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
### B. Definitions
### =========================================================================

# Number of cores to be used during parallel operations
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# SBATCH param
args <- commandArgs(trailingOnly = TRUE)
arrayID <- as.numeric(args[1])

### =========================================================================
### C- Species Data
### =========================================================================

# Target species arrayID
species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

ispi_name <- species[arrayID]

# Load species data
sp_dat_file <- file.path(scr_path, "outputs", "d0_datasets", "base", ispi_name, paste0(ispi_name, ".rds"))

sp_dat <- readRDS(sp_dat_file)
group_name <- sp_dat$group

cat(paste0("Starting global-level covariate extraction and selection for ", ispi_name, "...\n"))

### =========================================================================
### D- Covariate Data
### =========================================================================

# Retrieve lists of candidate covariates and cov_info table
lr_file <- file.path(w_path, "tmp", "settings", "ref_covariates.rds")
lr <- readRDS(lr_file)
cov_info <- lr$cov_info

# Refine global (glo) set
cov_info_glo <- cov_info[cov_info$level == "glo"]

if (n_levels > 1) {
cov_info_reg <- cov_info[cov_info$level == "reg"]

# Thin glo for matchin glo/reg files only
reg_files <- gsub("reg", "", cov_info_reg$file, fixed = T)
glo_files <- gsub("glo", "", cov_info_glo$file, fixed = T)
cov_info_glo <- cov_info_glo[na.omit(match(reg_files, glo_files)),]
}

# Subset with expert-filtered candidate covariates, if available
if (file.exists(expert_table)) {
  
  expert_tab <- try(fread(expert_table, sep="|"), silent = TRUE)

  if (!inherits(expert_tab, "try-error")) {
  
    # Ensure `group_name` exists before filtering
    if (group_name %in% colnames(expert_tab)) {
	expert_tab <- expert_tab[rowSums(expert_tab[, ..group_name] == 1) > 0]

    # Refine global covariates
      cov_info_glo <- merge(cov_info_glo, expert_tab, by = intersect(names(cov_info_glo), names(expert_tab)))
      lr_glo <- cov_info_glo$file

    } else {
      cat("Warning: Column '", group_name, "' not found in expert table. No filtering applied.\n", sep = "")
    }
  } else {
    cat("Warning: Expert table could not be read.\n")
  }
} else {
  lr_glo <- cov_info_glo$file
  cat("Warning: Expert table file not found. Using all candidate covariates.\n")
}

# Retrieve glo reference rasters
rsts_ref_file <- file.path(w_path, "tmp", "settings", "ref_rasters.rds")
rsts_ref <- readRDS(rsts_ref_file)
rsts_ref$rst_glo <- unwrap(rsts_ref$rst_glo)

### =========================================================================
### E - Covariate extraction
### =========================================================================

# Extract global data
pseu.abs_i_glo <- nsdm.bigextract(
  cov = gsub("\\.tif$", ".fst", lr_glo),
  data = sp_dat$pseu.abs_i_glo,
  rst_ref = rsts_ref$rst_glo,
  cov_info = cov_info_glo,
  t_match = as.logical(tmatch_glo),
  nsplits = ncores
)

# Filter unmatched covariates for GLO
filtered_out_layers <- pseu.abs_i_glo$filtered_out_layers
pseu.abs_i_glo <- unlist(pseu.abs_i_glo$data)
if (!is.null(filtered_out_layers$unmatched_covariates) && length(filtered_out_layers$unmatched_covariates) > 0) {
  cov_info_glo_filtered <- cov_info_glo[!gsub("\\.tif$", "", basename(cov_info_glo$file)) %in% filtered_out_layers$unmatched_covariates, ]
} else {
  cov_info_glo_filtered <- cov_info_glo
}

# Save a copy for later use
pseu.abs_i_glo_copy<-pseu.abs_i_glo

# Combine and scale environmental variables
env_vars <- scale(pseu.abs_i_glo@env_vars)
pseu.abs_i_glo@env_vars <- data.frame(env_vars)

# Define weights
wi <- which(pseu.abs_i_glo@pa == 1)
wt <- rep(1, length(pseu.abs_i_glo@pa))
wt[wi] <- round((length(pseu.abs_i_glo@pa) - length(wi)) / length(wi))
wt[wi][wt[wi] == 0] <- 1  # Ensure no zero weights for presences

# E.6 Save modeling set
l <- list(
  group = group_name,
  pseu.abs_i_glo_copy = pseu.abs_i_glo_copy,
  pseu.abs_i_glo = pseu.abs_i_glo,
  weights = wt,
  env_vars = env_vars
)

# Save dataset
nsdm.savethis(
  l,
  species_name = ispi_name,
  compression = TRUE,
  save_path = file.path(scr_path, "outputs/d0_datasets/glo")
)

cat("Modeling dataset prepared\n")


### =========================================================================
### F- Covariate Selection
### =========================================================================

counter <- 0
while (TRUE) {
  counter <- counter + 1
  
  # F.1 Step 1: Filtering for collinearity
  cat("Covariate selection S1: Filtering...\n")
  cov.filter_i <- try(
    covsel.filter(
      pseu.abs_i_glo@env_vars,
      pseu.abs_i_glo@pa,
      variables = cov_info_glo_filtered$variable,
      categories = cov_info_glo_filtered$cada,
      weights = wt,
      corcut = cor_cut
    ), 
    silent = TRUE
  )
  
  # F.2 Step 2: Model-specific embedding
  cat("Covariate selection S2: Embedding...\n")
  cov.embed_i <- try(
    covsel.embed(
      cov.filter_i,
      pseu.abs_i_glo@pa,
      weights = wt,
      ncov = if (!"esm" %in% mod_algo) {
        ceiling(log2(table(pseu.abs_i_glo@pa)["1"])) - 1
      } else {
        ncov_esm
      },
      maxncov = max_thre,
      nthreads = ncores
    ), 
    silent = TRUE
  )
  
  # Break loop if successful or after 5 attempts
  if (counter > 5 || !inherits(cov.embed_i, "try-error")) break
  cat(paste0("An error occurred; attempt number ", counter, " for covariate selection...\n"))
}

# F.3 Finalize: Create raster stack for selected covariates
files <- if (n_levels > 1) {
  na.omit(cov_info_reg$file[match(
    gsub("glo", "reg", cov.embed_i$ranks_2$covariate), 
    gsub("\\.tif$", "", basename(cov_info_reg$file))
  )])
} else {
  na.omit(cov_info_glo$file[match(
    cov.embed_i$ranks_2$covariate, 
    gsub("\\.tif$", "", basename(cov_info_glo$file))
  )])
} 
		
stk <- lapply(files, function(f) { 
  toMemory(rast(f))
})

# Extract SpatRasters from nested lists
raster_list <- unlist(stk, recursive = TRUE)

# Stack them into a single SpatRaster
raster_stack<-rast(raster_list)
names(raster_stack) <- gsub("reg_", "glo_", sub("\\.[^.]+$", "", basename(files)))
raster_stack <- wrap(raster_stack)

# F.4 Save Results
pseu.abs_i_glo@env_vars <- cov.embed_i$covdata

nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i_glo, filter=names(cov.filter_i), embed=cov.embed_i, covstk=raster_stack),
              species_name=ispi_name,
              save_path=file.path(scr_path,"outputs/d1_covsels/glo"))
			  
cat("Covariate preparation and selection completed successfully.\n")
cat("Finished.\n")