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

  # --- Safe coercion helper ---
  safe_as_char <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) return(x)
    if (is.factor(x)) return(as.character(x))
    if (is.numeric(x)) return(as.character(x))
  }

  # --- Coerce inputs safely to character ---
  model_name   <- safe_as_char(model_name)
  species_name <- safe_as_char(species_name)
  tag          <- safe_as_char(tag)

  # --- Build save directory (exclude NULL/empty elements) ---
  dir_parts <- c(save_path, species_name, model_name)
  dir_parts <- dir_parts[nzchar(dir_parts)]                # drop empty strings
  save_this_path <- do.call(file.path, as.list(dir_parts)) # join parts safely
  suppressWarnings(dir.create(save_this_path, recursive = TRUE))

  # --- Build clean base file name ---
  parts <- c(species_name, tag)
  parts <- parts[nzchar(parts)]
  base_name <- paste(parts, collapse = "_")

  # --- Save depending on format and object type ---
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


