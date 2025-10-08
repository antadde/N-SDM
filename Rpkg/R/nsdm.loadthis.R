#' nsdm.loadthis
#'
#' Easy loading of nsdm objects (.rds or .psv)
#'
#' @param model_name A character string indicating the name of the modelling algorithm related to the target object to be loaded
#' @param species_name A character string indicating the name of the species related to the target object to be loaded
#' @param tag A character string indicating the tag for the hyperparameter X modelling algorithm combination of the target object to be loaded
#' @param read_path A character string indicating the upper path where the target object is stored (above model and species)
#' @param format A character string, either "rds" (default) or "psv", specifying the file format to load
#'
#' @return The desired nsdm object
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.loadthis <- function(model_name = NULL,
                          species_name = NULL,
                          tag = NULL,
                          read_path,
                          format = "rds") {

  # Construct path to object
  read_this_path <- paste(read_path, species_name, model_name, sep = "/")

  # Build base filename
  if (is.null(tag)) {
    base_name <- paste(species_name, model_name, sep = "_")
  } else {
    base_name <- paste(species_name, tag, sep = "_")
  }

  # Normalize potential underscores before dot
  file_rds <- gsub("_\\.", ".", paste0(read_this_path, "/", base_name, ".rds"))
  file_psv <- gsub("_\\.", ".", paste0(read_this_path, "/", base_name, ".psv"))

  # Determine which file to load
  if (tolower(format) == "psv" && file.exists(file_psv)) {
    object <- read.table(file_psv, sep = "|", header = TRUE, stringsAsFactors = FALSE)

  } else if (file.exists(file_rds)) {
    object <- readRDS(file_rds)

  } else if (file.exists(file_psv)) {
    # Fallback if .rds not found but .psv exists
    object <- read.table(file_psv, sep = "|", header = TRUE, stringsAsFactors = FALSE)
    warning("RDS file not found; loaded PSV file instead.")

  } else {
    stop("No valid file found for: ", base_name)
  }

  return(object)
}
