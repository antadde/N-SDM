#############################################################################
## Script: 0_mainPRE_B
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

# Retrieve reference rasters
rsts_ref <- readRDS(file.path(w_path, "tmp", "settings", "ref_rasters.rds"))
rsts_ref$rst_glo <- unwrap(rsts_ref$rst_glo)
if (n_levels > 1) rsts_ref$rst_reg <- unwrap(rsts_ref$rst_reg)
if (n_levels > 1) rsts_ref$rst_reg_gloproj <- unwrap(rsts_ref$rst_reg_gloproj)

### =========================================================================
### C. Species Data
### =========================================================================

# Load species data
sp_dat <- readRDS(file.path(w_path, "tmp", "settings", "ref_species_occurences.rds"))

species_file <- file.path(w_path, "tmp", "settings", "tmp_species_list.txt")

species <- readLines(species_file)

# Loop over species
for (ispi_name in species) {

cat(paste0("Starting dataset preparation for: ", ispi_name, "...\n"))

# =============================================================================
# D. Sample Background points
# =============================================================================

# ----------------------
# D.1 GLO Background points
# ----------------------
sp_dat_glo <- sp_dat$spe_glo_dis[sp_dat$spe_glo_dis$species == ispi_name, ]

# Select bias raster if weighted background sampling is applied
if ("group" %in% colnames(sp_dat_glo)) {
  # Extract unique group names
  selected_group <- unique(sp_dat_glo$group)
} else {
  selected_group <- NULL  # No group information available
}

sampf_glo <- if (back_strat == "random_w" & !is.null(selected_group)) {
  selected_files <- list.files(
    file.path(data_path, "background", "random_w"),
    pattern = paste0(selected_group, "_glo\\.tif$"),  
    full.names = TRUE
  )

  if (length(selected_files) > 0) {
    cat("Selected glo background bias file: ", paste(basename(selected_files), collapse = ", "), "\n", sep = "")
	selected_files
  } else {
    cat("No glo background bias file found for group: ", selected_group, "\n", sep = "")
	NULL
  }
}

# Select filtering shapefile if biogeographical background sampling is applied
biog_glo <- if (back_strat == "random_biog") {
  selected_files <- list.files(
    file.path(data_path, "background", "random_biog"),
    pattern = "_glo\\.shp$",  
    full.names = TRUE
  )

  if (length(selected_files) > 0) {
    cat("Selected glo biogeographic background filtering file: ", paste(basename(selected_files), collapse = ", "), "\n", sep = "")
	selected_files
  } else {
    cat("No glo biogeographic background filtering file found \n")
	NULL
  }
}

# Generate background points
pseu.abs_i_glo <- nsdm.absences(
  n = n_back,
  rst_ref = rsts_ref$rst_glo,
  rst_background_weight = sampf_glo,
  shp_background_filter = biog_glo,
  type = pa_po_glo,
  pres = sp_dat_glo,
  rst_reg_gloproj = if (inherits(rsts_ref$rst_reg_gloproj, "SpatRaster")) rsts_ref$rst_reg_gloproj else NULL,
  level = "glo"
)

cat(paste0("GLO dataset prepared (n_occ=", sum(pseu.abs_i_glo@pa == 1), ")...\n"))

# Skip further processing if no occurrences are found
if (sum(pseu.abs_i_glo@pa == 1) == 0) next

# ----------------------
# D.2 REG Background points
# ----------------------
if (n_levels > 1) {
  sp_dat_reg <- sp_dat$spe_reg_dis[sp_dat$spe_reg_dis$species == ispi_name, ]

# Select bias raster if weighted background sampling is applied
if ("group" %in% colnames(sp_dat_reg)) {
  # Extract unique group names and normalize them
  selected_group <- unique(sp_dat_reg$group)
} else {
  selected_group <- NULL  # No group information available
}
sampf_reg <- if (back_strat == "random_w" & !is.null(selected_group)) {
  selected_files <- list.files(
    file.path(data_path, "background", "random_w"),
    pattern = paste0(selected_group, "_reg\\.tif$"),  
    full.names = TRUE
  )
  
  if (length(selected_files) > 0) {
    cat("Selected reg background bias file: ", paste(basename(selected_files), collapse = ", "), "\n", sep = "")
	selected_files
  } else {
    cat("No reg background bias file found for group: ", selected_group, "\n", sep = "")
	NULL
  }
}

# Select filtering shapefile if biogeographical background sampling is applied
biog_reg <- if (back_strat == "random_biog") {
  selected_files <- list.files(
    file.path(data_path, "background", "random_biog"),
    pattern = "_reg\\.shp$",  
    full.names = TRUE
  )

  if (length(selected_files) > 0) {
    cat("Selected biogeographic background filtering file: ", paste(basename(selected_files), collapse = ", "), "\n", sep = "")
	selected_files
  } else {
    cat("No reg biogeographic background filtering file found \n")
	NULL
  }
}

  # Generate background points
  pseu.abs_i_reg <- nsdm.absences(
    n = n_back,
    rst_ref = rsts_ref$rst_reg,
    rst_background_weight = sampf_reg,
	shp_background_filter = biog_reg,
    type = pa_po_reg,
    pres = sp_dat_reg,
	level = "reg"
  )

  cat(paste0("REG dataset prepared (n_occ=", sum(pseu.abs_i_reg@pa == 1), ")...\n"))
  
  if (sum(pseu.abs_i_reg@pa == 1) == 0) next
}

# Skip further processing if no occurrences are found

# ----------------------
# D.3 Save as Shapefiles
# ----------------------
# Create output directory
dir <- file.path(scr_path, "outputs", "d0_datasets", "shp", ispi_name)
dir.create(dir, recursive = TRUE)

# Convert GLO points to sf object and save
glo_pts_pres <- st_as_sf(
  data.frame(pseu.abs_i_glo@xy[pseu.abs_i_glo@pa == 1, , drop = FALSE]),
  coords = c("X", "Y"),  
  crs = st_crs(rsts_ref$rst_glo)
)
st_write(glo_pts_pres, dsn=dir, layer=paste0(ispi_name, "_glo_pres.shp"), driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)

glo_pts_abs <- st_as_sf(
  data.frame(pseu.abs_i_glo@xy[pseu.abs_i_glo@pa == 0, , drop = FALSE]),
  coords = c("X", "Y"),  
  crs = st_crs(rsts_ref$rst_glo)
)
st_write(glo_pts_abs, dsn=dir, layer=paste0(ispi_name, "_glo_abs.shp"), driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)

# Convert REG points to sf object and save (if applicable)
if (n_levels > 1) {
  reg_pts_pres <- st_as_sf(
    data.frame(pseu.abs_i_reg@xy[pseu.abs_i_reg@pa == 1, , drop = FALSE]),
    coords = c("X", "Y"),  
    crs = st_crs(rsts_ref$rst_reg)
  )
  st_write(reg_pts_pres, dsn=dir, layer=paste0(ispi_name, "_reg_pres.shp"), driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)
  
  reg_pts_abs <- st_as_sf(
    data.frame(pseu.abs_i_reg@xy[pseu.abs_i_reg@pa == 0, , drop = FALSE]),
    coords = c("X", "Y"),  
    crs = st_crs(rsts_ref$rst_reg)
  )
  st_write(reg_pts_abs, dsn=dir, layer=paste0(ispi_name, "_reg_abs.shp"), driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE) 

}

# ----------------------
# D.4 Generate training and testing sets
# ----------------------

if (use_spatial_stratification) {
  if (use_random_global) {
    cat("Training and testing sets built using spatial clustering for REG and random splits for GLO (nsdm.preps3 with glo_random = TRUE)", "\n")
    
    if (n_levels > 1) {
      all_sets <- nsdm.preps3(
        pseu_reg   = pseu.abs_i_reg,
        pseu_glo   = pseu.abs_i_glo,
        n_reps     = reps,
        glo_random = TRUE
      )
    } else {
      all_sets <- nsdm.preps3(
        pseu_glo   = pseu.abs_i_glo,
        n_reps     = reps,
        glo_random = TRUE
      )
    }
    
  } else {
    cat("Training and testing sets built using spatially clustered sampling (nsdm.preps3)", "\n")
    
    if (n_levels > 1) {
      all_sets <- nsdm.preps3(
        pseu_reg = pseu.abs_i_reg,
        pseu_glo = pseu.abs_i_glo,
        n_reps   = reps
      )
    } else {
      all_sets <- nsdm.preps3(
        pseu_glo = pseu.abs_i_glo,
        n_reps   = reps
      )
    }
  }

} else {
  cat("Training and testing sets built using random sampling (nsdm.preps2)", "\n")
  
  if (n_levels > 1) {
    all_sets <- nsdm.preps2(
      pseu_reg = pseu.abs_i_reg, 
      pseu_glo = pseu.abs_i_glo,
      n_reps   = reps
    )
  } else {
    all_sets <- nsdm.preps2(
      pseu_glo = pseu.abs_i_glo,
      n_reps   = reps
    )
  }
}


  ### =========================================================================
  ### E. Save
  ### =========================================================================
  
  if (n_levels > 1) {
    l <- list(
      group = selected_group,
      pseu.abs_i_glo = pseu.abs_i_glo,
      pseu.abs_i_reg = pseu.abs_i_reg,
	  all_sets = all_sets
    )
  } else {
    l <- list(
      group = selected_group,
      pseu.abs_i_glo = pseu.abs_i_glo,
	  all_sets = all_sets
    )
  }

  nsdm.savethis(
  object = l,
  species_name = ispi_name,
  save_path = file.path(scr_path, "outputs", "d0_datasets", "base")
  )
}

cat("Finished!", "\n")  