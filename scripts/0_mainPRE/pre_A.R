#############################################################################
## Script: 0_mainPRE_A
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Main N-SDM Settings
### =========================================================================

# Set permissions for new files
Sys.umask(mode = "000")

# Load and retrieve main settings
settings <- read.csv("./settings/settings.psv", sep = "|", stringsAsFactors = FALSE)
parameters <- settings$Parameter
values <- settings$Value

# Iterate through parameters and assign values
for (i in seq_along(parameters)) {
  param <- parameters[i]
  value <- trimws(values[i])

  # Case 1: NULL or empty/whitespace-only
  if (toupper(value) == "NULL" || grepl("^\\s*$", value)) {
    do.call("<-", list(param, NULL))

  # Case 2: comma-separated list (e.g., glm,gam,gbm)
  } else if (grepl(",", value)) {
    items <- trimws(unlist(strsplit(value, ",")))
    do.call("<-", list(param, items))

  # Case 3: numeric value
  } else if (!is.na(suppressWarnings(as.numeric(value)))) {
    do.call("<-", list(param, as.numeric(value)))

  # Case 4: plain string
  } else {
    do.call("<-", list(param, value))
  }
}

# Define data path
if (is.null(data_path)) data_path <- file.path(w_path, "data")

# Define species data
spe_glo <- list.files(file.path(data_path, "species", "glo"), 
                      full.names = TRUE, pattern = "\\.psv$")

if (n_levels > 1) {
  spe_reg <- list.files(file.path(data_path, "species", "reg"), 
                        full.names = TRUE, pattern = "\\.psv$")
}

# Define parameter grid, expert tables and forced species list
param_grid <- file.path(w_path, "scripts", "settings", param_grid)

if (!is.null(expert_table)) expert_table <- file.path(w_path, "scripts", "settings", expert_table)

if (!is.null(forced_species)) forced_species <- file.path(w_path, "scripts", "settings", forced_species)

# Check and refine masks
if (!is.null(mask_pred_reg))  mask_pred_reg <- file.path(data_path, "masks", mask_pred_reg)
if (!is.null(mask_pred_glo))  mask_pred_glo <- file.path(data_path, "masks", mask_pred_glo)

# Identify posthoc nesting methods
posthoc_nesting_methods <- intersect(nesting_methods, c("multiply", "multiplyw", "average", "averagew"))
if (length(posthoc_nesting_methods) > 0) {
  nesting_methods <- unique(c(
    intersect(nesting_methods, "covariate"),
    "posthoc"
  ))
}

# Save settings and clean up
rm(settings, parameters, values, i)
save.image(file.path(w_path, "tmp", "settings", "ref_nsdm_settings.RData"))
cat("N-SDM settings defined...\n")

### =========================================================================
### B- Prepare Covariate Data
### =========================================================================

# global reproducibility seed
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

# Set library path
.libPaths(Rlib_path)

# Load nsdm package
require(nsdm2)

# Retrieve number of cores
ncores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# Clean pre-prepared fst files
if (clear_fst) {
fst_files <- list.files(data_path, pattern = "\\.fst$", recursive = TRUE, full.names = TRUE)
invisible(file.remove(fst_files))
}

# Initial data check
nsdm.datacheck(
  data_dir = data_path,
  n_levels = n_levels
)

cat("Initial data check procedure completed..\n")

# Generate and open covariate info table
nsdm.covinfo(
  cov_path = file.path(data_path, "covariates"), 
  save_path = file.path(w_path, "tmp", "settings"))
  cov_info <- fread(
  file.path(w_path, "tmp", "settings", "ref_covariates_available.psv"),
  na.strings = c("NA", "", "NULL")
)

# List available covariates
lr_glo <- cov_info[level == "glo", file]

lr <- if (n_levels > 1) {
  lr_reg <- cov_info[level == "reg", file]
  c(lr_glo, lr_reg)
} else {
  lr_glo
}

# Set FST threads to avoid overloading the system
threads_fst(nr_of_threads = 1)

# Identify files that already have FST versions
fst_exists <- sapply(lr, function(f) file.exists(gsub(".tif", ".fst", f, fixed = TRUE)))
lr2fst <- lr[!fst_exists]  # Filter files that need FST conversion

# Create missing FST files
if (length(lr2fst) > 0) {  
  invisible(mclapply(lr2fst, function(c) {
    r <- rast(c)                  # Load raster
    r_df <- as.data.frame(r, na.rm=FALSE)      # Convert raster to data frame
    write.fst(r_df,               # Write data frame to FST format
              gsub(".tif", ".fst", c, fixed = TRUE), 
              compress = 75)      # Use compression level 75
  }, mc.cores = ncores))
  
} else {
invisible()
}

# Check .fst file sizes compared to .tif files
size_check <- lapply(lr, function(c) {
  tif_size <- file.info(c)$size
  fst_file <- gsub(".tif", ".fst", c, fixed = TRUE)
  fst_size <- if (file.exists(fst_file)) file.info(fst_file)$size else NA
  list(tif_size = tif_size, fst_size = fst_size, redo = !is.na(fst_size) && fst_size < (tif_size / 2))
})

# Identify .fst files that need to be recreated
to_redo <- lr[sapply(size_check, function(x) x$redo)]
if (length(to_redo) > 0) {
  
  invisible(mclapply(to_redo, function(c) {
    r <- rast(c)                  # Load raster
    r_df <- as.data.frame(r)      # Convert raster to data frame
    write.fst(r_df,               # Rewrite data frame to FST format
              gsub(".tif", ".fst", c, fixed = TRUE), 
              compress = 75)      # Use compression level 75
  }, mc.cores = ncores))
  
} else {
invisible()
}

cat("FST conversions and verification done.\n")

# Create regional and global reference rasters
## Global: Use bioclim layer if available, else use the first valid one
rst_glo <- if (length(grep("bio1", lr_glo, value = TRUE)) > 0) {
  rast_file <- grep("bio1", lr_glo, value = TRUE)[1]
  rast(rast_file)
} else {
  rast_file <- lr_glo[1]
  rast(rast_file)
}

## Regional: Use bioclim layer if available, else use the first valid one
if (n_levels > 1) {
  rst_reg <- if (length(grep("bio1", lr_reg, value = TRUE)) > 0) {
    rast_file <- grep("bio1", lr_reg, value = TRUE)[1]
    rast(rast_file)
  } else {
    rast_file <- lr_reg[1]
    rast(rast_file)
  }
  rst_reg_gloproj <- project(rst_reg, crs(rst_glo))
}

# Save covariate data settings
## Reference rasters
ref_rasters <- if (n_levels > 1) {
  list(rst_reg = wrap(rst_reg), rst_glo = wrap(rst_glo), rst_reg_gloproj = wrap(rst_reg_gloproj))
} else {
  list(rst_glo = wrap(rst_glo))
}
saveRDS(ref_rasters, file.path(w_path, "tmp", "settings", "ref_rasters.rds"))

## Covariates list and info
covariates_list <- if (n_levels > 1) {
  list(lr_reg = lr_reg, lr_glo = lr_glo, cov_info = cov_info)
} else {
  list(lr_glo = lr_glo, cov_info = cov_info)
}
saveRDS(covariates_list, file.path(w_path, "tmp", "settings", "ref_covariates.rds"))

cat("Covariate settings defined...\n")

### =========================================================================
### C- Prepare and Filter Species List for Modeling
### =========================================================================
# Check if .fst versions of regional and global species data exist or create them
threads_fst(nr_of_threads = ncores)

if (n_levels > 1) {
spe_reg_fst <- gsub("\\.psv$", ".fst", spe_reg)
   # Check if the .fst version already exists
  if (!file.exists(spe_reg_fst)) {
    if (grepl("\\.psv$", spe_reg, ignore.case = TRUE)) {
      # If spe_reg is a .psv file
      spe_reg_pts_dat <- fread(spe_reg)
    } else {
      stop("Unsupported file format for spe_reg. Only .psv is supported.")
    }
    
    # Write the extracted data to .fst
    write_fst(spe_reg_pts_dat, spe_reg_fst, compress = 0)
  }
}

spe_glo_fst <- gsub("\\.psv$", ".fst", spe_glo)
if (!file.exists(spe_glo_fst)) {
  if (grepl("\\.psv$", spe_glo, ignore.case = TRUE)) {
    # If spe_glo is a .psv file
    spe_glo_pts_dat <- fread(spe_glo)
  } else {
    stop("Unsupported file format for spe_glo. Only .psv is supported.")
  }
  
  # Write the extracted data to .fst
  write_fst(spe_glo_pts_dat, spe_glo_fst, compress = 0)
}

# Load regional species fst data
if (n_levels > 1) {
  spe_reg_fst_path <- gsub("\\.psv$", ".fst", spe_reg, ignore.case = TRUE)
  spe_reg_fst <- read_fst(spe_reg_fst_path)
   
  # Get unique species names for the regional dataset
  spe_reg_names <- unique(spe_reg_fst$species)
  cat(paste0("Initial number of species in REG species dataset is: ", length(spe_reg_names), "...\n"))
}

# Load global species fst data
spe_glo_fst_path <- gsub("\\.psv$", ".fst", spe_glo, ignore.case = TRUE)
spe_glo_fst <- read_fst(spe_glo_fst_path)

# Get unique species names for the global dataset
spe_glo_names <- unique(spe_glo_fst$species)
cat(paste0("Initial number of species in GLO species dataset is: ", length(spe_glo_names), "...\n"))

# Intersect global and regional species lists if both levels exist
if (n_levels > 1) {
  species <- sort(intersect(spe_glo_names, spe_reg_names))
} else {
  species <- spe_glo_names
}

cat(paste0("Number of remaining species after intersecting GLO and REG species dataset is: ", length(species), "...\n"))

# Handle forced species if provided
if (exists("forced_species") && length(forced_species) > 0) {
  if (file.exists(forced_species)) {
    forced_species <- readLines(forced_species)
  }
  if (is.character(forced_species)) {
    species <- forced_species
  }
}

cat(paste0("Total number of species considered for this N-SDM run is: ", length(species), "...\n"))

### =========================================================================
### D- Spatiotemporal Disaggregation and Preparation of Species Data
### =========================================================================

# Subset regional species data
if (n_levels > 1) {
  # Filter regional data for selected species
  spe_reg_pts <- spe_reg_fst[spe_reg_fst$species %in% species, ]
  spe_reg_pts$sid <- paste("reg", 1, 1:nrow(spe_reg_pts), sep="_")
  
  cat(paste0("Ready for spatiotemporal disaggregation of regional data for ", length(species), " species...\n"))

  # Disaggregate regional species data
  if (!as.logical(tmatch_reg)) {
    spe_reg_pts$year <- NULL  # Remove temporal matching if not required
  }
    spe_reg_pts
  if (disag_reg) {
    spe_reg_dis <- nsdm.disaggregate(
      pres = spe_reg_pts,
      rst = rst_reg,
      thindist = thin_dist_reg,
      thinyear = thin_time_reg,
      min_occ = min_occ_reg,
      ncores = ncores
    )
  } else {
    # Use original data if no disaggregation is requested
    spe_reg_dis <- spe_reg_pts
  }

  # Update species list based on regional disaggregation
  species <- unique(spe_reg_dis$species)
}

# Subset global species data
spe_glo_pts <- spe_glo_fst[spe_glo_fst$species %in% species, ]
spe_glo_pts$sid <- paste("glo", 1, 1:nrow(spe_glo_pts), sep="_")
	
cat(paste0("Ready for spatiotemporal disaggregation of global data for ", length(species), " species...\n"))

# Disaggregate global species data
  if (!as.logical(tmatch_glo)) {
  spe_glo_pts$year <- NULL  # Remove temporal matching if not required
}
if (disag_glo) {
  spe_glo_dis <- nsdm.disaggregate(
    pres = spe_glo_pts,
    rst = rst_glo,
    thindist = thin_dist_glo,
    thinyear = thin_time_glo,
    min_occ = min_occ_glo,
    ncores = ncores
  )
} else {
  # Use original data if no disaggregation is requested
  spe_glo_dis <- spe_glo_pts
}

# Use both reg and glo or glo-only species data in the global model if requested
if (n_levels > 1) {
	if (glo_use_reg) {
cat(paste0("Combined use of global and regional species data in the global model requested...\n"))
  # Reproject regional dataset to match global CRS
  spe_reg_dis_reproj <- st_transform(spe_reg_dis, crs = st_crs(spe_glo_dis))

  # Harmonize columns between both datasets
  all_cols <- union(names(spe_glo_dis), names(spe_reg_dis_reproj))
  for (col in setdiff(all_cols, names(spe_glo_dis))) spe_glo_dis[[col]] <- NA
  for (col in setdiff(all_cols, names(spe_reg_dis_reproj))) spe_reg_dis_reproj[[col]] <- NA

  # Bind them
  sp_combined <- rbind(spe_glo_dis[, all_cols], spe_reg_dis_reproj[, all_cols])
  
  # Re-disaggregate after combination on global grid
  if (disag_glo) {
  cat(paste0("Ready for spatiotemporal disaggregation of combined global and regional data for ", length(species), " species...\n"))
  sp_combined <- nsdm.disaggregate(
    pres = sp_combined,
    rst = rst_glo,
    thindist = thin_dist_glo,
    thinyear = thin_time_glo,
    min_occ = min_occ_glo,
    ncores = ncores
  )
} else {
  # Use original data if no disaggregation is requested
  sp_combined <- sp_combined
}
}
}

# Save species data settings
l3 <- if (n_levels > 1 && !is.null(glo_use_reg) && glo_use_reg) {
  list(spe_reg_dis = spe_reg_dis, spe_glo_dis = sp_combined)  # Include both regional and global data
} else if (n_levels > 1) {
	list(spe_reg_dis = spe_reg_dis, spe_glo_dis = spe_glo_dis)  # Use separate global and regional data
} else {
  list(spe_glo_dis = spe_glo_dis)  # Only global data
}
saveRDS(l3, file = file.path(w_path, "tmp", "settings", "ref_species_occurences.rds"))
writeLines(species, file.path(w_path, "tmp", "settings", "ref_species_list.txt"))

# Calculate and save the number of cluster runs required
spe_runs <- length(splitIndices(length(species), ceiling(length(species) / n_mx_spe)))
writeLines(as.character(spe_runs), con = file.path(w_path, "tmp", "settings", "ref_species_runs.txt"))

# Completion message
cat("Finished!...\n")

