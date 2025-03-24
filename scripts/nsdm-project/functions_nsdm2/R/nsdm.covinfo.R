#' nsdm.covinfo
#'
#' Generate Covariate Information Table
#'
#' This function generates a covariate information table in `.xlsx` format based on the provided covariate paths and names.
#'
#' @param cov_path A character string specifying the root directory where covariates are stored.
#' @param save_path A character string specifying the destination path to save the `.xlsx` covariate information table.
#' @param time_cov A character vector specifying the categories or datasets that include temporally dynamic covariates (with scenarios).
#' @param focal_cat A character vector specifying the categories or datasets that include covariates with focal statistics.
#'
#' @return An `.xlsx` file containing the covariate information table.
#' @details This function processes the covariate paths to identify dynamic and focal covariates, and saves the resulting information table in Excel format.
#'
#' @export
#' @author Antoine Adde (antoine.adde@eawag.ch)

nsdm.covinfo <- function(
  cov_path,
  save_path,
  time_cov, 
  focal_cov
) { 
  # Helper function to extract start and end years
  extract_years <- function(sub_period) {
    if (grepl("_", sub_period)) {
      years <- strsplit(sub_period, split = "_")[[1]]
      list(start_year = years[1], end_year = years[2])
    } else {
      list(start_year = sub_period, end_year = sub_period)
    }
  }
  
  # List all available layers and extract their information
  setwd(cov_path)
  full_pred_list <- system(paste("find -L", cov_path, "-type f -name '*.tif'"), intern = TRUE)
  pred_list <- gsub(cov_path, "", full_pred_list)
  pred_list <- gsub("^/", "", pred_list)

    # Initialize pred_table as a list with the same length as pred_list
  pred_table <- list(
    level = rep(NA, length(pred_list)),
    category = rep(NA, length(pred_list)),
    dataset = rep(NA, length(pred_list)),
    cada = rep(NA, length(pred_list)),
    period = rep(NA, length(pred_list)),
    variable = rep(NA, length(pred_list)),
    scenario = rep(NA, length(pred_list)),
    sub_period = rep(NA, length(pred_list)),
    attribute = rep(NA, length(pred_list)),
    start_year = rep(NA, length(pred_list)),
    end_year = rep(NA, length(pred_list)),
    focal = rep(NA, length(pred_list)),
    file = rep(NA, length(pred_list))
  )  
  
  # Extract generic information
  pred_table$level <- sapply(strsplit(pred_list, split = '/'), `[[`, 1)
  pred_table$category <- sapply(strsplit(pred_list, split = '/'), `[[`, 2)
  pred_table$dataset <- sapply(strsplit(pred_list, split = '/'), `[[`, 3)
  pred_table$cada <- paste(pred_table$category, pred_table$dataset, sep = "_")
  
    
  # Process each covariate
  for (i in seq_along(pred_list)) {
    split <- strsplit(pred_list[i], split = '/')[[1]]
    cada <- pred_table$cada[i]
    
    # Determine if covariate is available for scenarions
    if (cada %in% time_cov) {
      if (grepl("scenario", pred_list[i])) {
        pred_table$period[i] <- "scenario"
        pred_table$scenario[i] <- split[length(split) - 2]
        pred_table$sub_period[i] <- split[length(split) - 3]
      } else if (grepl("present", pred_list[i])) {
        pred_table$period[i] <- "present"
        pred_table$sub_period[i] <- split[length(split) - 2]
      } else {
        pred_table$period[i] <- "present"
        pred_table$sub_period[i] <- split[length(split) - 2]
      }
      
      # Extract start and end years
      years <- extract_years(pred_table$sub_period[i])
      pred_table$start_year[i] <- years$start_year
      pred_table$end_year[i] <- years$end_year
      
      # Extract variable and attribute
      pred_table$variable[i] <- split[length(split) - 1]
      pred_table$attribute[i] <- str_match(pred_list[i], paste0(pred_table$variable[i], "_\\s*(.*?)\\s*.tif"))[2]
      
    } else {
      # Static covariates
      pred_table$period[i] <- "present"
      pred_table$variable[i] <- split[length(split) - 1]
      pred_table$attribute[i] <- str_match(pred_list[i], paste0("_", pred_table$variable[i], "_\\s*(.*?)\\s*.tif"))[2]
    }
    
    # Add focal statistics if applicable
    if (cada %in% focal_cov) {
      focal <- sub("^.+_", "", tools::file_path_sans_ext(full_pred_list[i]))
      suppressWarnings(true_focal <- as.numeric(focal))
      pred_table$focal[i] <- true_focal
    }
    
    # Add file path
    pred_table$file[i] <- full_pred_list[i]
  }
  
  # Clean and format the table
  pred_table <- data.frame(pred_table, stringsAsFactors = FALSE)
  pred_table$attribute <- sapply(seq_len(nrow(pred_table)), function(i) {
    gsub(paste0("_", pred_table$focal[i]), "", pred_table$attribute[i])
  })
  pred_table$attribute[pred_table$focal == pred_table$attribute] <- NA
  pred_table[is.na(pred_table)] <- "NA"
  
  # Save covariate table
  row.names(pred_table) <- NULL
fwrite(pred_table, file.path(save_path, "predictors_available.psv"), sep = "|")
}