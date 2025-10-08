#' nsdm.savethis
#'
#' Easy saving of .rds or .psv objects
#'
#' @param object Object to be saved
#' @param model_name A character string indicating the name of the target modelling algorithm
#' @param species_name A character string indicating the name of the target taxon
#' @param tag A character string indicating the tag for the target hyperparameter X modelling algorithm combination
#' @param compression Logical (TRUE/FALSE) for saveRDS "compress" argument
#' @param save_path A character string indicating the upper path where to save the object
#' @param format A character string, either "rds" (default) or "psv", specifying the output format
#'
#' @return Saved object in .rds or .psv format
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.savethis <- function(object,
                          model_name = NULL,
                          species_name = NULL,
                          tag = NULL,
                          compression = FALSE,
                          save_path,
                          format = "rds") {

  # Create save directory
  save_this_path <- paste(save_path, species_name, model_name, sep = "/")
  suppressWarnings(dir.create(save_this_path, recursive = TRUE))

  # Define output filename (without extension)
  if (is.null(tag)) {
    base_name <- paste(species_name, model_name, sep = "_")
  } else {
    base_name <- paste(species_name, tag, sep = "_")
  }

  # Handle GBM separately
  if ("lgb.Booster" %in% class(object)) {
    f <- paste0(save_this_path, "/", base_name, ".rds")
    lgb.save(object, file = f)

  } else if (tolower(format) == "psv") {
    # Save as PSV (pipe-separated values)
    f <- paste0(save_this_path, "/", base_name, ".psv")

    if (is.data.frame(object)) {
      write.table(object, file = f, sep = "|", quote = FALSE, row.names = FALSE)
    } else {
      warning("PSV format is intended for data frames; object saved as RDS instead.")
      f <- paste0(save_this_path, "/", base_name, ".rds")
      saveRDS(object, file = f, compress = compression)
    }

  } else {
    # Default RDS saving
    f <- paste0(save_this_path, "/", base_name, ".rds")
    saveRDS(object, file = f, compress = compression)
  }

  return(invisible(f))
}