#' nsdm.ensemble
#'
#' Ensemble prediction surfaces obtained by individual algorithms.
#'
#' @param model_names A character vector indicating the modelling algorithms to ensemble.
#' @param species_name A character string indicating the name of the taxon.
#' @param map_path A character string indicating the path where algorithm-specific prediction surfaces are stored.
#' @param score_path A character string indicating the path where algorithm-specific evaluation metrics are stored.
#' @param discthre A numeric value specifying the threshold for `weight_metric`, below which an algorithm is discarded. 
#'                 If `NULL`, no models are discarded.
#' @param weighting Logical. If `TRUE`, individual algorithms are weighted by `weight_metric`.
#' @param weight_metric A character string indicating the name of the metric used for weighting or discarding models.
#' @param level A character string indicating the name of the level considered (e.g., "glo" or "reg").
#' @param nesting_name A character string indicating the name of the nesting strategy for which predicted values are provided.
#' @param scenar_name A character string indicating the name of the scenario for which predicted values are provided.
#' @param period_name A character string indicating the name of the period for which predicted values are provided.
#'
#' @return A named list of two `SpatRaster` objects:
#' \describe{
#'   \item{"ensemble"}{A `SpatRaster` representing the average ensemble prediction.}
#'   \item{"ensemble_cv"}{A `SpatRaster` representing the coefficient of variation of the ensemble.}
#' }
#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.ensemble <- function(model_names, species_name, level=NA, nesting_name=NA, scenar_name=NA, period_name=NA, map_path, score_path=NULL, discthre=NULL, weighting=FALSE, weight_metric="Score"){
  
# Retrieve prediction maps
stack_map <- list()

for (i in seq_along(model_names)) {
  model_name <- model_names[i]
  if(model_name=="esm") {
  full_map_paths <- list.files(file.path(map_path, species_name), pattern = "^esm", full.names = TRUE)
  map_files <- list.files(full_map_paths, pattern="\\.tif$", full.names = TRUE, recursive = TRUE)
  } else {
  full_map_path <- file.path(map_path, species_name, model_name)
  map_files <- list.files(full_map_path, pattern="\\.tif$", full.names = TRUE, recursive = TRUE) }

  if (length(map_files) > 0) {
    map2 <- lapply(map_files, rast)
    map <- rast(map2)
    stack_map <- c(stack_map, list(map))  
  }
}

stack_map<-rast(stack_map)

# Initialize results table
res <- data.frame(matrix(nrow = nlyr(stack_map), ncol = 3))
colnames(res) <- c("model_name", "score", "discard")

for (i in seq_along(model_names)) {
  model_name <- model_names[i]
  full_score_path <- file.path(score_path, species_name, model_name, paste0(species_name, "_", model_name, ".psv"))
   
    Score <- fread(full_score_path, sep="|") 

    # Identify selected model(s) and retrieve scores
    if (model_name == "esm") {
      esm_ix <- grep("_esm", names(stack_map))  # Find ESM-related layers
      esm_names <- paste("esm", stri_extract_first_regex(names(stack_map)[esm_ix], "[0-9]+"), sep="-")
      score_val <- Score[Metric == weight_metric, ..esm_names]
      res[esm_ix, "score"] <- as.numeric(score_val)
      res[esm_ix, "model_name"] <- esm_names
     
	 } else {
        vals <- unlist(Score[Metric == weight_metric, -1, with = FALSE])
        score_val  <- max(vals)
        res[i, "score"] <- as.numeric(score_val)
        res[i, "model_name"] <- model_name
    }
}

# Check if models fulfill discard threshold, and remove if needed
if (!"esm" %in% model_names) {
  if (!is.null(discthre)) {
    res[,"discard"] <- res[,"score"] < discthre
  } else {
    res[,"discard"] <- FALSE 
  }

  # Discard models below the threshold
  if (!is.null(discthre)) {
    discard_indices <- which(as.logical(res[,"discard"]))
    if (length(discard_indices) > 0) {
keep_indices <- setdiff(seq_len(nlyr(stack_map)), discard_indices)
stack_map <- terra::subset(stack_map, keep_indices)
      res <- res[-discard_indices, , drop = FALSE]
    }
  }
}

# Perform weighted or unweighted ensemble
if (weighting) {
  discard_indices <- which(as.logical(res[,"discard"]))
  if (length(discard_indices) > 0) {
keep_indices <- setdiff(seq_len(nlyr(stack_map)), discard_indices)
stack_map <- terra::subset(stack_map, keep_indices)
  }
  ensemble <- terra::app(stack_map, mean, w = as.numeric(res[,"score"]), na.rm = TRUE)
} else {
  ensemble <- terra::app(stack_map, mean, na.rm=TRUE)
}

ensemble_mn <- ensemble  # Store mean ensemble

# Compute coefficient of variation using terra
rasterstack_sd_fast <- function(x) {
  s0 <- nlyr(x)
  s1 <- terra::subset(x, 1)
  s2 <- s1^2
  for (ri in 2:s0) {
    r <- terra::subset(x, ri)
    s1 <- s1 + r
    s2 <- s2 + r^2
  }
  sqrt((s0 * s2 - s1 * s1) / (s0 * (s0 - 1)))
}

ensemble_sd <- rasterstack_sd_fast(stack_map)
ensemble_cv <- (ensemble_sd / ensemble_mn) * 100

# Final rounding and type conversion
ensemble_cv <- terra::round(ensemble_cv)
ensemble_mn <- terra::round(ensemble_mn)

storage.mode(values(ensemble_cv)) <- "integer"
storage.mode(values(ensemble_mn)) <- "integer"

# Rename layers
ensemble_name <- paste(ispi_name, level, nesting_name, scenar_name, period_name, "ensemble", sep = "_")
ensemble_cv_name <- paste(ispi_name, level, nesting_name, scenar_name, period_name, "ensemble_cv", sep = "_")

names(ensemble_mn) <- gsub("_NA", "", ensemble_name)
names(ensemble_cv) <- gsub("_NA", "", ensemble_cv_name)

return(list(ensemble=toMemory(ensemble_mn), ensemble_cv=toMemory(ensemble_cv)))
}
