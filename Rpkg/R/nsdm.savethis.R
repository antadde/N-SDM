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

  # Build save directory (independent of tag)
  dir_parts <- c(save_path, species_name)
  if (!is.null(model_name) && nzchar(as.character(model_name))) {
    dir_parts <- c(dir_parts, as.character(model_name))
  }
  save_this_path <- paste(dir_parts, collapse = "/")
  suppressWarnings(dir.create(save_this_path, recursive = TRUE))

  # --- Clean and safe base name ---
  parts <- c(species_name, tag)
  parts <- parts[!sapply(parts, is.null)]          # drop NULLs
  parts <- as.character(parts)                     # ensure character type
  parts <- parts[nzchar(parts)]                    # drop empty strings safely
  base_name <- paste(parts, collapse = "_")

  # Handle GBM separately
  if ("lgb.Booster" %in% class(object)) {
    f <- file.path(save_this_path, paste0(base_name, ".rds"))
    lgb.save(object, file = f)

  } else if (tolower(format) == "psv") {
    f <- file.path(save_this_path, paste0(base_name, ".psv"))

    if (is.data.frame(object)) {
      write.table(object, file = f, sep = "|", quote = FALSE, row.names = FALSE)
    } else {
      warning("PSV format is intended for data frames; object saved as RDS instead.")
      f <- sub("\\.psv$", ".rds", f)
      saveRDS(object, file = f, compress = compression)
    }

  } else {
    f <- file.path(save_this_path, paste0(base_name, ".rds"))
    saveRDS(object, file = f, compress = compression)
  }

  return(invisible(f))
}

