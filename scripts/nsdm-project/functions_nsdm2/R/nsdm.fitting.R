#' nsdm.fitting
#'
#' Model fitting
#'
#' @md
#' @param x An `nsdm.pseudoabsences` object used as the master table for feature values and SIDs
#' @param sets A named list of replicates produced upstream, each replicate contains `reg_train`, `glo_train`
#' @param mod_args A list of class `multi.input`, one element per model family, each with slots `@mod` a function name, `@args` a list of arguments, `@weight` logical, `@tag` a short label
#' @param level Character string, `"glo"` or `"reg"`, selects which subset of `sets` to use for SIDs
#' @param ncores Integer, number of cores for `mclapply`
#' @param tmp_path Character, directory for temporary files such as LightGBM saves
#'
#' @return An `nsdm.fit` object with fitted models stored in `@fits`, grouped by `@tag`, each tag holds a list of per replicate fits
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.fitting <- function(x,
                         sets,
                         mod_args,
                         level,
                         ncores = 1L,
                         tmp_path) {

  reps  <- length(sets)

  # container for all model families across replicates
  modis <- list()

  # loop over model families
  for (j in seq_along(mod_args)) {

    mod_args_j <- mod_args[[j]]
    
    # fit in parallel across replicates
    modit <- parallel::mclapply(seq_len(reps), function(rep_id) {

      # fetch train SIDs for this level and replicate
      sub_sets  <- sets[[rep_id]][grep(level, names(sets[[rep_id]]))]
      train_sid <- sub_sets[[paste0(level, "_train")]]@sid

      # index rows in x
      i_tr <- match(train_sid, x@sid); i_tr <- i_tr[!is.na(i_tr)]

      # build train frames
      train <- cbind(
        data.frame(sid = x@sid[i_tr], split = "train", Presence = x@pa[i_tr]),
        x@env_vars[i_tr, , drop = FALSE]
      )

      df_train <- train[, c("Presence", colnames(x@env_vars)), drop = FALSE]

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
        mod_args_j@args$data <- train[, colnames(x@env_vars), drop = FALSE]
        mod_args_j@args$p    <- df_train$Presence
        fit <- try(do.call(mod_args_j@mod, mod_args_j@args), silent = TRUE)
        if (inherits(fit, "try-error")) {
          mod_args_bis <- mod_args_j@args
          mod_args_bis$addsamplestobackground <- TRUE
          fit <- do.call(mod_args_j@mod, mod_args_bis)
        }

      } else if (identical(mod_args_j@mod, "lgb.train")) {
        # LightGBM
        taxon <- tryCatch(x@meta$taxon, error = function(e) "unknown_taxon")
        tmp_path_gbm <- file.path(tmp_path, "gbm")
        if (!dir.exists(tmp_path_gbm)) dir.create(tmp_path_gbm, recursive = TRUE)

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

        # save LightGBM model to disk
        out_file <- paste0(tmp_path_gbm,"/",taxon,"_rep",rep_id,"_mod",j,"_",level,".rds")
        lightgbm::lgb.save(fit, out_file)

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

      } else if (identical(mod_args_j@mod, "ranger")) {
        # ranger
        if (!is.null(wt)) mod_args_j@args$case.weights <- wt
        mod_args_j@args$data <- df_train[, colnames(x@env_vars), drop = FALSE]
        mod_args_j@args$data$Presence <- as.factor(df_train$Presence)
        fit <- do.call(mod_args_j@mod, mod_args_j@args)

      } else {
        stop(sprintf("Unknown model engine: %s", mod_args_j@mod))
      }

      return(fit)
    }, mc.cores = ncores)

    # name replicates
    names(modit) <- sprintf("replicate_%02d", seq_len(reps))

    # store under model tag
    modis[[mod_args_j@tag]] <- modit
  }

  # wrap into nsdm.fit object
  out <- nsdm.fit()
  out@fits <- modis
  out@meta$env_vars = paste(colnames(x@env_vars), collapse = ", ")
  return(out)
}
