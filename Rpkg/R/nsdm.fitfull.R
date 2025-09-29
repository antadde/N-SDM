#' nsdm.fitfull
#'
#' Model fitting (full dataset)
#'
#' @md
#' @param x An `nsdm.pseudoabsences` object used as the master table for feature values and SIDs
#' @param mod_args A list of class `multi.input`, one element per model family, each with slots `@mod` a function name, `@args` a list of arguments, `@weight` logical, `@tag` a short label
#' @param level Character string, `"glo"` or `"reg"`
#'
#' @return An `nsdm.fit` object with fitted models stored in `@fits`, grouped by `@tag`, each tag holds a list of per replicate fits
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.fitfull <- function(x,
                         mod_args,
                         level) {

  # container for all model families
  modis <- list()

  # loop over model families
  for (j in seq_along(mod_args)) {

    mod_args_j <- mod_args[[j]]
    
      # build train frame
      df_train <- cbind(
        data.frame(Presence = x@pa),
        x@env_vars
      )

      # weights for presence upweighting if requested
      wt <- NULL
      if (isTRUE(mod_args_j@weight)) {
        wi <- which(df_train$Presence == 1)
        wt <- rep(1, nrow(df_train))
        if (length(wi) > 0 && length(wi) < length(wt)) {
          wt[wi] <- max(1L, round((length(wt) - length(wi)) / length(wi)))
        }
      }

      # choose model engine
      if (identical(mod_args_j@mod, "maxnet")) {
        # Maxent like with maxnet
        mod_args_j@args$data <- df_train[, colnames(x@env_vars), drop = FALSE]
        mod_args_j@args$p    <- df_train$Presence
        fit <- try(do.call(mod_args_j@mod, mod_args_j@args), silent = TRUE)
        if (inherits(fit, "try-error")) {
          mod_args_bis <- mod_args_j@args
          mod_args_bis$addsamplestobackground <- TRUE
          fit <- do.call(mod_args_j@mod, mod_args_bis)
        }

      } else if (identical(mod_args_j@mod, "lgb.train")) {
        # LightGBM
        if (!is.null(wt)) {
          mod_args_j@args$data <- lightgbm::lgb.Dataset(
            data  = as.matrix(df_train[, colnames(x@env_vars), drop = FALSE]),
            label = df_train$Presence,
            weight = wt
          )
        } else {
          mod_args_j@args$data <- lightgbm::lgb.Dataset(
            data  = as.matrix(df_train[, colnames(x@env_vars), drop = FALSE]),
            label = df_train$Presence
          )
        }

        fit <- do.call(mod_args_j@mod, mod_args_j@args)

      } else if (mod_args_j@mod %in% c("glm", "mgcv_gam", "mgcv_fx", "esm")) {
        # GLM or GAM family
        if (!is.null(wt)) mod_args_j@args$weights <- wt
        mod_args_j@args$data <- df_train
        fit <- do.call(mod_args_j@mod, mod_args_j@args)

      } else if (identical(mod_args_j@mod, "randomForest")) {
        # randomForest
        if (!is.null(wt)) mod_args_j@args$weights <- wt
        mod_args_j@args$data <- df_train[, colnames(x@env_vars), drop = FALSE]
        mod_args_j@args$data$Presence <- as.factor(df_train$Presence)
        fit <- do.call(mod_args_j@mod, mod_args_j@args)

      } else {
        stop(sprintf("Unknown model engine: %s", mod_args_j@mod))
      }

   # store under model tag
    modis[[mod_args_j@tag]] <- fit
  }

  # wrap into nsdm.fit object
  out <- nsdm.fit()
  out@fits <- modis
  out@meta$env_vars = paste(colnames(x@env_vars), collapse = ", ")
  return(out)
}
