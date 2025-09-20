#############################################################################
## Script: 2_mainREG_A
## Purpose: Covariate extraction and selection
## Date: 20-05-2022
## Author: Antoine Adde
#############################################################################

### =========================================================================
### A. Preparation
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

cat(paste0("Starting regional-level covariate extraction and selection for ", ispi_name, "...\n"))

### =========================================================================
### D- Covariate data
### =========================================================================
# Retrieve lists of candidate covariates and cov_info table
lr_file <- file.path(w_path, "tmp", "settings", "ref_covariates.rds")
lr <- readRDS(lr_file)
cov_info <- lr$cov_info

# Add catvar column
cov_info[, catvar := paste(category, variable, sep = "_")]

# Identify catvars common to both 'reg' and 'glo' levels
common_catvar <- intersect(
  cov_info[level == "reg", unique(catvar)],
  cov_info[level == "glo", unique(catvar)]
)

# Get file paths for those non common catvars
common_catvar_files <- cov_info[!catvar %in% common_catvar, file]

# Filter lr_reg: exclude scenario paths and retain only those in common_catvar_files
lr_reg <- lr$lr_reg[!grepl("/scenario/", lr$lr_reg) & lr$lr_reg %in% common_catvar_files]
cov_info_reg <- cov_info[match(lr_reg, cov_info$file), , drop = FALSE]

# Subset with expert-filtered candidate covariates, if available
if (file.exists(expert_table)) {
  expert_tab <- try(fread(expert_table, sep="|"), silent = TRUE)

  if (!inherits(expert_tab, "try-error")) {
  
    # Ensure `group_name` exists before filtering
    if (group_name %in% colnames(expert_tab)) {
	expert_tab <- expert_tab[rowSums(expert_tab[, ..group_name] == 1) > 0]

    # Refine regional covariates if applicable
      if (n_levels > 1) {
        cov_info_reg <- merge(cov_info_reg, expert_tab, by = intersect(names(cov_info_reg), names(expert_tab)), all.x = TRUE)
        lr_reg <- cov_info_reg$file
      }
    } else {
      cat("Warning: Column '", group_name, "' not found in expert table. No filtering applied.\n", sep = "")
    }
  } else {
    cat("Warning: Expert table could not be read.\n")
  }
} else {
  cat("Warning: Expert table file not found. Using all candidate covariates.\n")
}

# Retrieve glo and reg reference rasters
rsts_ref_file <- file.path(w_path, "tmp", "settings", "ref_rasters.rds")
rsts_ref <- readRDS(rsts_ref_file)
rsts_ref$rst_reg <- unwrap(rsts_ref$rst_reg)

cat(paste0('Covariate data listed \n'))

### =========================================================================
### E- Prepare modelling dataset
### =========================================================================

# E.1 Retrieve predictions from mainGLO
glo_out_f <- list.files(
  file.path(scr_path, "outputs", "d8_ensembles/glo", ispi_name), 
  pattern = "\\.tif$", 
  full.names = TRUE
)

if (length(glo_out_f) == 0) stop("Error: No file found in GLO output folder!")

# E.2 Extract new habitat covariates
pseu.abs_i <- nsdm.bigextract(
  cov = c(gsub("\\.tif", ".fst", lr_reg), glo_out_f),
  data = sp_dat$pseu.abs_i_reg,
  rst_ref = rsts_ref$rst_reg,
  cov_info = cov_info_reg,
  t_match = as.logical(tmatch_reg),
  tmatch_scheme = tmatch_scheme_reg,
  nzvt = 16,
  nsplits = ncores
)

# Update cov_info_reg
filtered_out_layers <- pseu.abs_i$filtered_out_layers
pseu.abs_i <- unlist(pseu.abs_i$data)
if (!is.null(filtered_out_layers$unmatched_covariates) && length(filtered_out_layers$unmatched_covariates) > 0) {
  
  # Extract filenames without the `.tif` extension to match `cov_info$file`
  cov_info_reg_filtered <- cov_info_reg[!gsub("\\.tif$", "", basename(cov_info_reg$file)) %in% filtered_out_layers$unmatched_covariates, ]

} else {
  cov_info_reg_filtered <- cov_info_reg  # No changes if no unmatched covariates
}

# Scale env_vars
colnames(pseu.abs_i@env_vars)[ncol(pseu.abs_i@env_vars)] <- "mainGLO"
env_vars <- scale(pseu.abs_i@env_vars)
pseu.abs_i@env_vars <- as.data.frame(env_vars)

# E.3 Define weights
wi <- which(pseu.abs_i@pa == 1)
wt <- rep(1, length(pseu.abs_i@pa))
wt[wi] <- round((length(pseu.abs_i@pa) - length(wi)) / length(wi))
if (any(wt[wi] == 0)) wt[wi] <- 1

# E.4 Save modelling set
save_path <- file.path(scr_path, "outputs", "d0_datasets/reg")

nsdm.savethis(
  object = list(
    pseu.abs_i = pseu.abs_i,
    weights = wt,
    cov_info = cov_info_reg_filtered,
    env_vars = env_vars,
    group = sp_dat$group
  ),
  species_name = ispi_name,
  compression = TRUE,
  save_path = save_path
)

cat("Modelling dataset prepared \n")

### =========================================================================
### F- Covariate selection with mainGLO forced
### =========================================================================
covstk_res <- list()
covdata_res <- list()

if ("covariate" %in% nesting_methods) {
  counter <- 0
  while (TRUE) {
    counter <- counter + 1
    
    # Step 1: Filtering for collinearity
    cat('Covariate selection with mainGLO forced S1: Filtering...\n')
    cov.filter_i <- try(
      covsel.filter(
        pseu.abs_i@env_vars,
        pseu.abs_i@pa,
        variables = c(cov_info_reg_filtered$variable, "mainGLO"),
        categories = c(cov_info_reg_filtered$cada, "mainGLO"),
        weights = wt,
        force = c("mainGLO"),
        corcut = cor_cut
      ), silent = TRUE
    )
   
    # Step 2: Model-specific embedding
    cat('Covariate selection with mainGLO forced S2: Embedding...\n')
    cov.embed_i <- try(
      covsel.embed(
        cov.filter_i,
        pseu.abs_i@pa,
        weights = wt,
        algorithms = c("glm", "gam", "rf"),
        force = c("mainGLO"),
        ncov = if (!"esm" %in% mod_algo) {
          ceiling(log2(table(pseu.abs_i@pa)['1'])) - 1
        } else {
          ncov_esm
        },
        maxncov = max_thre,
        nthreads = ncores
      ), silent = TRUE
    )
    
    # If embedding succeeds or counter reaches max attempts, break loop
    if (!inherits(cov.embed_i, "try-error") || counter > 5) break
    cat(paste0("An error occurred; attempt number ", counter, " for covariate selection...\n"))
  }
  
  # Step 3: Finalizing the selection process
  if (!inherits(cov.embed_i, "try-error")) {
    files <- na.omit(cov_info_reg_filtered$file[match(
      cov.embed_i$ranks_2$covariate, 
      gsub("\\.tif$", "", basename(cov_info_reg_filtered$file))
    )])
    
    stk <- lapply(files, function(f) { 
      toMemory(rast(f))
    })
    
    glo_out <- rast(glo_out_f)
    names(glo_out) <- "mainGLO"
    
	raster_stack<-rast(stk)
    names(raster_stack)<-basename(file_path_sans_ext(files))
	covstk_res[["cov"]] <- wrap(c(raster_stack, glo_out))
    covdata_res[["cov"]] <- cov.embed_i$covdata
  } else {
    cat("Covariate embedding failed after 5 attempts.\n")
  }
}

### =========================================================================
### G- Covariate selection without main GLO
### =========================================================================
if (any(c("multiply") %in% nesting_methods)) {
  counter <- 0
  while (TRUE) {
    counter <- counter + 1
    
    # Remove mainGLO
    pseu.abs_i@env_vars <- subset(pseu.abs_i@env_vars, select = -c(mainGLO))
    
    # Step 1: Filtering for collinearity
    cat('Covariate selection S1: Filtering...\n')
    cov.filter_i <- try(
      covsel.filter(
        pseu.abs_i@env_vars,
        pseu.abs_i@pa,
        variables = cov_info_reg_filtered$variable,
        categories = cov_info_reg_filtered$cada,
        weights = wt,
        corcut = cor_cut
      ), silent = TRUE
    )
    
    # Step 2: Model-specific embedding
    cat('Covariate selection S2: Embedding...\n')
    cov.embed_i <- try(
      covsel.embed(
        cov.filter_i,
        pseu.abs_i@pa,
        weights = wt,
        algorithms = c("glm", "gam", "rf"),
        ncov = if (!"esm" %in% mod_algo) {
          ceiling(log2(table(pseu.abs_i@pa)['1'])) - 1
        } else {
          ncov_esm
        },
        maxncov = max_thre,
        force = NULL,
        nthreads = ncores
      ), silent = TRUE
    )
    
    # If embedding succeeds or counter reaches max attempts, break loop
    if (!inherits(cov.embed_i, "try-error") || counter > 5) break
    cat(paste0("An error occurred; attempt number ", counter, " for covariate selection...\n"))
  }
  
  # Step 3: Finalizing the selection process
  if (!inherits(cov.embed_i, "try-error")) {
    files <- na.omit(cov_info_reg_filtered$file[match(
      cov.embed_i$ranks_2$covariate, 
      gsub("\\.tif$", "", basename(cov_info_reg_filtered$file))
    )])
    
    stk <- lapply(files, function(f) { 
      toMemory(rast(f))
    })

    raster_stack<-rast(stk)
    names(raster_stack)<-basename(file_path_sans_ext(files))
    covstk_res[["mul"]] <- wrap(raster_stack)
    covdata_res[["mul"]] <- cov.embed_i$covdata
	
  } else {
    cat("Covariate embedding failed after 5 attempts.\n")
  }
}

# Save
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i,
              covstk=covstk_res,
			  covdata=covdata_res,
			  env_vars=env_vars),
              species_name=ispi_name,
              save_path=file.path(scr_path,"outputs/d1_covsels/reg"))

cat(paste0('Covariate selection for ',ispi_name,' done \n'))
cat(paste0('Finished \n'))