#' nsdm.datacheck
#'
#' Validates the data for an NSDM-compliant run.  
#' Specifically, this function checks that:
#' \itemize{
#'   \item All required subfolders (e.g., masks, species, covariates) are present.
#'   \item Covariate data follows the expected folder hierarchy.
#'   \item Covariate data follows the expected naming convention.
#'   \item For each level ('reg', 'glo'), a subset of raster files is sampled and checked for consistent extent, resolution, CRS, and grid alignment.
#'   \item The CRS and grid of the rasters are compared to a sample of species coordinates to ensure compatibility.
#' }
#'
#' @param data_dir Character. Path to the main data directory (must contain subfolders like 'covariates', 'species', 'masks').
#' @param n_levels Integer. Number of spatial levels included (1 for 'reg' only, 2 for both 'reg' and 'glo').
#'
#' @return Invisibly returns a data.frame with consistency results for sampled rasters (extent, resolution, CRS, grid).
#' If inconsistencies are found, the function stops and reports the problematic files or species mismatches.
#'
#' @author Antoine Adde (\email{antoine.adde@eawag.ch})
#' @export

nsdm.datacheck <- function(data_dir, n_levels) {

  cov_dir <- file.path(data_dir, "covariates")
  sp_dir  <- file.path(data_dir, "species")
  tif_paths <- list.files(cov_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
  dirs_with_data <- unique(dirname(tif_paths))
  valid_levels <- if (n_levels > 1) c("reg", "glo") else "reg"

  # Step 1: Check required subfolders
  required_dirs <- if (n_levels > 1) {
    c("masks", "species/glo", "species/reg", "covariates/glo", "covariates/reg")
  } else {
    c("masks", "species/reg", "covariates/reg")
  }

  missing_dirs <- required_dirs[!dir.exists(file.path(data_dir, required_dirs))]
  if (length(missing_dirs) > 0) {
    message("❌ Missing required folders:\n", paste(file.path(data_dir, missing_dirs), collapse = "\n"))
  }

  # Step 2: Check covariate folder structure
  is_valid_cov_path <- function(path) {
    rel_path <- gsub(paste0("^", cov_dir, "/?"), "", path)
    parts <- strsplit(rel_path, "/")[[1]]
    depth <- length(parts)
    parts[1] %in% valid_levels && depth >= 4 && depth <= 7
  }

  invalid_paths <- dirs_with_data[!vapply(dirs_with_data, is_valid_cov_path, logical(1))]
  if (length(invalid_paths) > 0) {
    message("🚨 Invalid covariate folder structure found in the following paths:\n", paste(invalid_paths, collapse = "\n"))
  }

  # Step 2a: Check filenames match folder structure
  mismatched_filenames <- sapply(tif_paths, function(path) {
    path_parts <- strsplit(gsub(paste0("^", cov_dir, "/"), "", path), "/")[[1]]
    expected_base <- paste(path_parts[-length(path_parts)], collapse = "_")
    !startsWith(basename(path), expected_base)
  })

  if (any(mismatched_filenames)) {
    message("🚨 Filenames not mirrored by folder structure:\n", paste(tif_paths[mismatched_filenames], collapse = "\n"))
  }

  # Step 2b: Check for special characters in folder names
  spec_char <- dirs_with_data[grepl("(?<!\\d)_|_(?!\\d)|[^a-zA-Z0-9/._]", dirs_with_data, perl = TRUE)]
  if (length(spec_char) > 0) {
    message("🚨 Suspicious use of special characters in folder names:\n", paste(spec_char, collapse = "\n"))
  }

  # Step 3: Run consistency checks per level
  for (level in valid_levels) {
    check_covariate_consistency(level, cov_dir, sp_dir)
  }

  message("✅ Data check procedure completed successfully.")
}

check_covariate_consistency <- function(level = "reg", cov_dir, sp_dir) {
  base_path <- file.path(cov_dir, level)
  dataset_paths <- list.dirs(base_path, recursive = TRUE, full.names = TRUE)
  dataset_paths <- dataset_paths[
    sapply(strsplit(dataset_paths, "/"), length) == (length(strsplit(base_path, "/")[[1]]) + 2)
  ]

  sampled_paths <- unlist(lapply(dataset_paths, function(ds_path) {
    tifs <- list.files(ds_path, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
    if (length(tifs) >= 2) sample(tifs, 2) else tifs
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
    compare_df[i, 2:5] <- list(
      ext(r) == ext(ref),
      all(res(r) == res(ref)),
      crs(r) == crs(ref),
      compareGeom(r, ref, stopOnError = FALSE)
    )
  }

  failed <- compare_df[!(compare_df$same_extent & compare_df$same_res & compare_df$same_crs & compare_df$same_grid), ]

  # Species data check
  sp_files <- list.files(file.path(sp_dir, level), pattern = "\\.psv$", recursive = TRUE, full.names = TRUE)
  if (length(sp_files) == 0) stop("❌ No species file found for level: ", level)
  species_data <- tryCatch(fread(sp_files[1]), error = function(e) NULL)

  if (is.null(species_data) || !all(c("X", "Y") %in% names(species_data))) {
    stop("❌ Could not read species file or missing required X/Y columns.")
  }

  sample_coords <- species_data[sample(nrow(species_data), min(100, nrow(species_data))), .(X, Y)]
  sample_coords <- sample_coords[complete.cases(sample_coords), ]
  sample_coords[, c("X", "Y") := .(as.numeric(X), as.numeric(Y))]

  points <- terra::vect(sample_coords, geom = c("X", "Y"), crs = crs(ref))
  vals <- tryCatch(terra::extract(ref, points), error = function(e) NULL)

  if (nrow(failed) > 0) {
    stop("❌ Inconsistent rasters (compared to ", basename(sampled_paths[1]), "):\n",
         paste(failed$path, collapse = "\n"))
  }

  if (is.null(vals)) {
    stop("❌ Species coordinates are not compatible with raster grid (CRS or extent mismatch).")
  }

  invisible(compare_df)
}
