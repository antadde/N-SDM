#' nsdm.datacheck
#'
#' Validates the data for an NSDM-compliant run.  
#' Specifically, this function checks that:
#' \itemize{
#'   \item All required subfolders (e.g., masks, species, covariates) are present based on the number of spatial levels (1 or 2).
#'   \item Covariate data follows the expected folder hierarchy: \code{level/category/dataset[/period/subperiod/scenario]/variable}.
#'   \item For each level ('reg', 'glo'), a subset of raster files is sampled and checked for consistent extent, resolution, CRS, and grid alignment.
#'   \item The CRS and grid of the rasters are compared to a sample of species coordinates to ensure compatibility.
#' }
#'
#' @param data_dir Character. Path to the main project data directory (must contain subfolders like 'covariates', 'species', 'masks').
#' @param cov_dir Character. Path to the covariate directory (typically \code{file.path(data_dir, "covariates")}).
#' @param sp_dir Character. Path to the species data directory (typically \code{file.path(data_dir, "species")}).
#' @param n_levels Integer. Number of spatial levels included in the project (1 for 'reg' only, 2 for both 'reg' and 'glo').
#'
#' @return Invisibly returns a data.frame with consistency results for sampled rasters (extent, resolution, CRS, grid).
#' If inconsistencies are found, the function stops and reports the problematic files or species mismatches.
#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.check_base_structure_and_covariates <- function(data_dir, cov_dir, sp_dir, n_levels = 2) {
  # Step 1: Check required subfolders
  required_dirs <- if (n_levels > 1) {
    c("masks", "species/glo", "species/reg", "covariates/glo", "covariates/reg")
  } else {
    c("masks", "species/reg", "covariates/reg")
  }

  missing_dirs <- required_dirs[!dir.exists(file.path(data_dir, required_dirs))]
  if (length(missing_dirs) > 0) {
    stop("❌ Missing required folders:\n", paste(file.path(data_dir, missing_dirs), collapse = "\n"))
  } else {
    message("✅ All required base data folders are present.")
  }

  # Step 2: Check covariate folder structure
  valid_levels <- if (n_levels > 1) c("reg", "glo") else c("reg")

  is_valid_cov_path <- function(path) {
    rel_path <- gsub(paste0("^", cov_dir, "/?"), "", path)
    parts <- strsplit(rel_path, "/")[[1]]
    depth <- length(parts)
    valid_level <- parts[1] %in% valid_levels
    valid_depth <- depth >= 4 && depth <= 7
    return(valid_level && valid_depth)
  }

  all_dirs <- list.dirs(cov_dir, recursive = TRUE, full.names = TRUE)
  dirs_with_data <- unique(dirname(list.files(cov_dir, pattern = "\\.(tif|fst)$", recursive = TRUE, full.names = TRUE)))
  invalid_paths <- dirs_with_data[!vapply(dirs_with_data, is_valid_cov_path, logical(1))]

  if (length(invalid_paths) > 0) {
    message("🚨 Invalid covariate folder structure found in the following paths:\n",
            paste(invalid_paths, collapse = "\n"))
  } else {
    message("✅ All covariate folders seem to respect the expected structure.")
  }

  # Step 3: Run consistency checks for each level
  for (level in valid_levels) {
    check_covariate_consistency_global(level, cov_dir, sp_dir)
  }
}

# Raster consistency check function
check_covariate_consistency_global <- function(level = "reg", cov_dir, sp_dir) {
  base_path <- file.path(cov_dir, level)

  dataset_paths <- list.dirs(base_path, recursive = TRUE, full.names = TRUE)
  dataset_paths <- dataset_paths[
    sapply(strsplit(dataset_paths, "/"), length) == (length(strsplit(base_path, "/")[[1]]) + 2)
  ]

  sampled_paths <- unlist(lapply(dataset_paths, function(ds_path) {
    tif_files <- list.files(ds_path, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
    if (length(tif_files) >= 2) sample(tif_files, 2) else tif_files
  }))

  if (length(sampled_paths) < 2) stop("❌ Not enough rasters to perform consistency check.")

  raster_list <- lapply(sampled_paths, function(x) tryCatch(rast(x), error = function(e) NULL))
  valid_rasters <- Filter(Negate(is.null), raster_list)
  if (length(valid_rasters) < 2) stop("❌ Failed to load enough valid rasters.")

  ref <- valid_rasters[[1]]
  compare_df <- data.frame(
    path = sampled_paths,
    same_extent = NA,
    same_res    = NA,
    same_crs    = NA,
    same_grid   = NA
  )

  for (i in seq_along(valid_rasters)) {
    r <- valid_rasters[[i]]
    compare_df$path[i]        <- sampled_paths[i]
    compare_df$same_extent[i] <- ext(r) == ext(ref)
    compare_df$same_res[i]    <- all(res(r) == res(ref))
    compare_df$same_crs[i]    <- crs(r) == crs(ref)
    compare_df$same_grid[i]   <- compareGeom(r, ref, stopOnError = FALSE)
  }

  failed <- compare_df[!(compare_df$same_extent & compare_df$same_res & compare_df$same_crs & compare_df$same_grid), ]

  # Check species data alignment
  sp_file <- list.files(file.path(sp_dir, level), pattern = "\\.psv$", recursive = TRUE, full.names = TRUE)[1]
  species_r <- tryCatch(fread(sp_file), error = function(e) NULL)

  if (is.null(species_r) || !"X" %in% names(species_r) || !"Y" %in% names(species_r)) {
    stop("❌ Could not read species data or missing X/Y columns.")
  }

  species_sample <- species_r[sample(nrow(species_r), min(100, nrow(species_r))), ]
  species_coords <- species_sample[, c("X", "Y")]
  species_coords <- species_coords[complete.cases(species_coords), ]
  species_coords$X <- as.numeric(species_coords$X)
  species_coords$Y <- as.numeric(species_coords$Y)

  points_spat <- terra::vect(species_coords, geom = c("X", "Y"), crs = crs(ref))
  vals <- tryCatch(terra::extract(ref, points_spat), error = function(e) NULL)

  # Final report
  if (nrow(failed) > 0) {
    stop("❌ Inconsistent rasters (compared to ", basename(sampled_paths[1]), "):\n",
         paste(failed$path, collapse = "\n"))
  }

  if (is.null(vals)) {
    stop("❌ Species coordinates are not compatible with the raster grid (e.g., CRS or extent mismatch).")
  }

  message("✅ All sampled rasters in level '", level, "' are spatially consistent.")
  message("✅'", level, "'Species coordinates are aligned with the reference raster grid.")

  invisible(compare_df)
}
