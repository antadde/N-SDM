#' nsdm.covinfo
#'
#' Generate Covariate Information Table
#'
#' This function generates a covariate information table based on the provided covariate paths and names.
#'
#' @param cov_path A character string specifying the directory where covariates are stored.
#' @param save_path A character string specifying the destination path to save the covariate information table.
#'
#' @return A `.psv` file containing the covariate information table.
#'
#' @export
#' @author Antoine Adde (antoine.adde@eawag.ch)

nsdm.covinfo <- function(
  cov_path,
  save_path){ 
  
 # List all available layers and extract their information
setwd(cov_path)
full_pred_list <- system(paste("find -L", cov_path, "-type f -name '*.tif'"), intern = TRUE)
pred_list <- gsub(paste0("^", cov_path, "/?"), "", full_pred_list)

# Initialize pred_table as a list
pred_table <- list(
  level     = rep(NA, length(pred_list)),
  category  = rep(NA, length(pred_list)),
  dataset   = rep(NA, length(pred_list)),
  cada      = rep(NA, length(pred_list)),
  year      = rep(NA, length(pred_list)),
  scenario  = rep(NA, length(pred_list)),
  variable  = rep(NA, length(pred_list)),
  attribute = rep(NA, length(pred_list)),
  focal     = rep(NA, length(pred_list)),
  file      = full_pred_list
)

# Extract general structure
splits <- strsplit(pred_list, "/")
pred_table$level    <- sapply(splits, `[[`, 1)
pred_table$category <- sapply(splits, `[[`, 2)
pred_table$dataset  <- sapply(splits, `[[`, 3)
pred_table$cada     <- paste(pred_table$category, pred_table$dataset, sep = "_")

# Process each covariate
for (i in seq_along(pred_list)) {
  split <- splits[[i]]
  file_path <- pred_list[i]

  # If scenario folder exists
  if (grepl("scenario", file_path)) {
    pred_table$scenario[i] <- split[length(split) - 2]
    pred_table$year[i]     <- as.numeric(split[length(split) - 3])

  } else if (!is.na(suppressWarnings(as.numeric(split[4])))) {
    pred_table$year[i] <- as.numeric(split[4])
  }

  # Extract variable (2nd last folder) and attribute (filename part after variable_)
  pred_table$variable[i] <- split[length(split) - 1]
  attribute_match <- regmatches(file_path, regexpr(paste0(pred_table$variable[i], "_[^/]+\\.tif$"), file_path))
  if(length(attribute_match)>0) pred_table$attribute[i] <- gsub(paste0(pred_table$variable[i], "_|\\.tif$"), "", attribute_match)

  # Extract focal radius if present (and remove attribute)
  focal_radius <- sub(".*_f([0-9]+)\\.tif$", "\\1", file_path)
  if (!is.na(suppressWarnings(as.numeric(focal_radius)))) {
  pred_table$focal[i] <- suppressWarnings(as.numeric(focal_radius))
  pred_table$attribute[i] <- NA
  }
}

  # Clean and format the table
  pred_table <- data.frame(pred_table, stringsAsFactors = FALSE)
  pred_table[is.na(pred_table)] <- "NA"
  
  # Save covariate table
  row.names(pred_table) <- NULL
  fwrite(pred_table, file.path(save_path, "predictors_available.psv"), sep = "|")
}