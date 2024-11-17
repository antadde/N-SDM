#############################################################################
## 1_mainGLO
## A: covariate extraction and covariate selection
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/1_mainGLO","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","tmp",gsub("/main/1_mainGLO","",getwd())),"/settings/nsdm-settings.RData"))

# Set permissions for new files
Sys.umask(mode="000")

# Set your working directory
setwd(w_path)

# Set lib path
.libPaths(lib_path)

# Load nsdm package
require(nsdm)

### =========================================================================
### B- Definitions
### =========================================================================
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# SBATCH param
args<-eval(parse(text=args))
arrayID<-eval(parse(text=arrayID))

### =========================================================================
### C- Species data
### =========================================================================
# Target species arrayID
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))
ispi_name<-species[arrayID]

# Load species data
sp_dat<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/base/",ispi_name,"/",ispi_name,".rds"))
group_name<-sp_dat$group

cat(paste0('Ready for modelling dataset preparation and covariate selection for ', ispi_name, '...\n'))

### =========================================================================
### D- Covariate Data
### =========================================================================
# Retrieve lists of candidate covariates and cov_info table
lr <- readRDS(paste0(w_path, "tmp/", project, "/settings/covariates-list.rds"))
cov_info <- lr$cov_info

# Refine global (glo) set to exclude future scenarios
lr_glo <- lr$lr_glo[grep("/future/", lr$lr_glo, invert = TRUE)]
cov_info_glo <- cov_info[match(lr_glo, cov_info$file), ]

# Refine regional (reg) set to include only CHELSA present scenarios
if (n_levels > 1) {
  lr_reg <- lr$lr_reg[grep("/chelsa/present", lr$lr_reg)]
  cov_info_reg <- cov_info[match(lr_reg, cov_info$file), ]
}

# Subset with expert-filtered candidate covariates, if available
expert_tab <- try(read_excel(expert_table, .name_repair = "minimal"), silent = TRUE)

if (!inherits(expert_tab, "try-error")) {
  # Clean column names by removing text in parentheses
  colnames(expert_tab) <- gsub("\\s*\\([^\\)]+\\)", "", colnames(expert_tab))

  # Remove duplicated columns, if any
  if (anyDuplicated(colnames(expert_tab)) > 0) {
    expert_tab <- expert_tab[, !duplicated(colnames(expert_tab))]
  }

  # Filter rows based on the group name condition
  expert_tab <- expert_tab[expert_tab[, group_name] == "1", ]

  # Refine global covariates
  cov_info_glo <- merge(cov_info_glo, expert_tab, by = intersect(names(cov_info_glo), names(expert_tab)))
  lr_glo <- cov_info_glo$file

  # Refine regional covariates if applicable
  if (n_levels > 1) {
    cov_info_reg <- merge(cov_info_reg, expert_tab, by = intersect(names(cov_info_reg), names(expert_tab)))
    lr_reg <- cov_info_reg$file
  }
}

# Retrieve glo and reg reference rasters
rsts_ref<-readRDS(paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))

### =========================================================================
### E- Covariate extraction
### =========================================================================
# E.1 GLO: Extract global data
pseu.abs_i_glo <- nsdm.bigextract(
  cov = gsub(".rds", ".fst", lr_glo),
  data = sp_dat$pseu.abs_i_glo,
  rst_ref = rsts_ref$rst_glo,
  cov_info = cov_info_glo,
  t_match = tmatch_glo,
  nsplits = ncores
)

# Create a copy of pseu.abs_i_glo for later use
pseu.abs_i_glo_copy <- pseu.abs_i_glo

# E.2 REG: Extract regional data (if applicable)
if (n_levels > 1) {
  pseu.abs_i_reg <- nsdm.bigextract(
    cov = gsub(".rds", ".fst", lr_reg),
    data = sp_dat$pseu.abs_i_reg,
    rst_ref = rsts_ref$rst_reg,
    cov_info = cov_info_reg,
    t_match = tmatch_reg,
    nsplits = ncores
  )
  
  # E.3 Combine GLO and REG data
  pseu.abs_i_glo@pa <- c(pseu.abs_i_glo@pa, pseu.abs_i_reg@pa)
  pseu.abs_i_glo@years <- c(pseu.abs_i_glo@years, pseu.abs_i_reg@years)
  pseu.abs_i_glo@xy <- rbind(pseu.abs_i_glo@xy, pseu.abs_i_reg@xy)

  # Extract variable names from pseu.abs_i_glo@env_vars
  glo_vars <- sapply(names(pseu.abs_i_glo@env_vars), function(x) sub(".*_", "", x))
  
  # Extract variable names from pseu.abs_i_reg@env_vars
  reg_vars <- sapply(names(pseu.abs_i_reg@env_vars), function(x) sub(".*_", "", x))
  
  # Find common variables between global and regional datasets
  common_vars <- intersect(glo_vars, reg_vars)
  
  # Match indices of common variables in both datasets
  glo_indices <- match(common_vars, glo_vars)
  reg_indices <- match(common_vars, reg_vars)
  
  # Update the names of pseu.abs_i_glo@env_vars using corresponding regional names
  names(pseu.abs_i_glo@env_vars)[glo_indices] <- names(pseu.abs_i_reg@env_vars)[reg_indices]
  
  # Scale environmental variables
  env_vars <- scale(rbind(pseu.abs_i_glo@env_vars, pseu.abs_i_reg@env_vars))
  pseu.abs_i_glo@env_vars <- data.frame(env_vars)
} else {
  # Scale environmental variables for GLO only
  env_vars <- scale(pseu.abs_i_glo@env_vars)
  pseu.abs_i_glo@env_vars <- data.frame(env_vars)
}

# E.4 Update cov_info table
cov_info_glo <- na.omit(
  cov_info_glo[match(
    colnames(pseu.abs_i_glo@env_vars),
    gsub(".rds", "", basename(cov_info_glo$file))
  ), ]
)

# E.5 Define weights
wi <- which(pseu.abs_i_glo@pa == 1)
wt <- rep(1, length(pseu.abs_i_glo@pa))
wt[wi] <- round((length(pseu.abs_i_glo@pa) - length(wi)) / length(wi))
if (any(wt[wi] == 0)) wt[wi] <- 1

# E.6 Save modeling set
l <- if (n_levels > 1) {
  list(
    group = sp_dat$group,
    pseu.abs_i_glo_copy = pseu.abs_i_glo_copy,
    pseu.abs_i_glo = pseu.abs_i_glo,
    pseu.abs_i_reg = pseu.abs_i_reg,
    weights = wt,
    env_vars = env_vars
  )
} else {
  list(
    group = sp_dat$group,
    pseu.abs_i_glo_copy = pseu.abs_i_glo_copy,
    pseu.abs_i_glo = pseu.abs_i_glo,
    weights = wt,
    env_vars = env_vars
  )
}

# Save the dataset
nsdm.savethis(
  l,
  species_name = ispi_name,
  compression = TRUE,
  save_path = paste0(scr_path, "/outputs/", project, "/d0_datasets/glo")
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
      variables = cov_info_glo$variable,
      categories = cov_info_glo$cada,
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
if (n_levels > 1) {
  stk <- try(
    nsdm.fastraster(
      files = na.omit(cov_info_reg$file[match(
        cov.embed_i$ranks_2$covariate, 
        gsub(".rds", "", basename(cov_info_reg$file))
      )]), 
      nsplits = ncores
    ), 
    silent = TRUE
  )
} else {
  stk <- try(
    nsdm.fastraster(
      files = na.omit(cov_info_glo$file[match(
        cov.embed_i$ranks_2$covariate, 
        gsub(".rds", "", basename(cov_info_glo$file))
      )]), 
      nsplits = ncores
    ), 
    silent = TRUE
  )
}

# F.4 Save Results
pseu.abs_i_glo@env_vars <- cov.embed_i$covdata

nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i_glo, filter=names(cov.filter_i), embed=cov.embed_i, covstk=stk),
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/glo"))

cat("Covariate selection completed successfully.\n")
cat("Finished.\n")