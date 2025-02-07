#' nsdm.map
#'
#' Fill an empty raster template with predicted values.
#'
#' @param template An empty raster layer to fill with predicted values.
#' @param nona_ix Numeric indices of cells that should be filled with predicted values.
#' @param pred A list where each element corresponds to a modeling algorithm and contains vectors of fitted values.
#' @param species_name A character string indicating the name of the taxon for which models are fitted.
#' @param model_name A character string indicating the name of the modeling algorithm used.
#' @param level A character string indicating the spatial level considered (e.g., "glo" for global, "loc" for local).
#' @param nesting_name A character string indicating the nesting strategy for which predicted values are provided.
#' @param scenar_name A character string indicating the scenario for which predicted values are provided.
#' @param period_name A character string indicating the time period for which predicted values are provided.
#'
#' @return A list of rasters filled with predicted values (for ESM models) or a single raster (for other models).
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.map <- function(template, nona_ix, pred, species_name, model_name, level = NA, nesting_name = NA, scenar_name = NA, period_name = NA) {

  pred_fit_f <- list()  # List to store results

  # Ensure pred is not empty before looping
  if (length(pred) == 0) {
    warning("Input 'pred' is empty. Returning NULL.")
    return(NULL)
  }

  for (m in seq_along(pred)) {  # Safe loop indexing

    # Retrieve model fit
    if (model_name == "esm") {
      if (!is.null(pred[[m]]$fit)) {
        fit <- pred[[m]]$fit
      } else {
        warning(paste("Missing 'fit' in ESM model", m))
        next  # Skip this iteration if 'fit' is missing
      }
    } else {
      if (!is.null(pred$fit)) {
        fit <- pred$fit
      } else {
        stop("Error: 'fit' not found in 'pred' for non-ESM model.")
      }
    }

    # Determine model-specific naming
    model_nameu <- if (model_name == "esm" && !is.null(names(pred)[m])) names(pred)[m] else model_name

    # Fill template raster
    pred_fit <- template
    pred_fit[] <- NA
    pred_fit[nona_ix] <- round(fit * 100)

    # Construct raster name safely
    pred_fit_name <- gsub("_NA", "", paste(species_name, model_nameu, level, nesting_name, scenar_name, period_name, sep = "_"))
    names(pred_fit) <- pred_fit_name

    storage.mode(pred_fit[]) <- "integer"  # Ensure integer format

    # Store result
    pred_fit_f[[m]] <- pred_fit
    names(pred_fit_f)[m] <- model_nameu
  }

  # Return full list for ESM models, single output for others
  if (model_name == "esm") {
    return(pred_fit_f)
  } else {
    return(pred_fit_f[[1]])
  }
}