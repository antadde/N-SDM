#' nsdm.eval
#'
#' Parallel evaluation of fitted models with common SDM metrics.
#'
#' @md
#' @param x An `nsdm.pseudoabsences` object that holds the master covariates, SIDs, and presence values
#' @param models An `nsdm.fit` object produced by your fitting step, models are in `models@fits`
#' @param sets A list of replicates, each replicate contains `reg_train`, `reg_test`, `glo_train`, `glo_test`
#' @param level Character string, either `"glo"` or `"reg"`, selects which SIDs to use for testing
#' @param ncores Integer, number of cores used by `parallel::mclapply`, default is `1`
#' @param tmp_path Character path for LightGBM model files, default is `tempdir()`
#'
#' @return An object of class `nsdm.evaluation` with performance per model tag and replicate
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.eval <- function(x, models, sets, level,
                            ncores = 1L, tmp_path){

  tmp_path_gbm <- file.path(tmp_path, "gbm")

  ## ------------------------
  ## Build testing covariates and presence vectors per replicate
  ## ------------------------
  cols <- colnames(x@env_vars)
  rep_ids <- seq_along(sets)

  testa <- vector("list", length(rep_ids))   # covariates only
  papa  <- vector("list", length(rep_ids))   # presence vector

  for (k in rep_ids) {
    sub_sets  <- sets[[k]][grep(level, names(sets[[k]]))]
    test_sid  <- sub_sets[[paste0(level, "_test")]]@sid

    i_te <- match(test_sid, x@sid)
    i_te <- i_te[!is.na(i_te)]

    df_test <- cbind(
      data.frame(Presence = x@pa[i_te]),
      x@env_vars[i_te, , drop = FALSE]
    )

    testa[[k]] <- df_test[, cols, drop = FALSE]
    papa[[k]]  <- df_test$Presence
  }
  names(testa) <- names(papa) <- sprintf("replicate_%02d", rep_ids)

  ## ------------------------
  ## Prepare evaluation container
  ## ------------------------
  out <- preva.meta(type = "evaluation") 

  taxon <- tryCatch(x@meta$taxon, error = function(e) "unknown_taxon")
  safe_taxon <- gsub("[^A-Za-z0-9_]+", "_", taxon)

  ## ------------------------
  ## Evaluate each model family across replicates
  ## ------------------------
  lisa <- list()

  for (j in seq_along(models@fits)) {
    tag_j <- names(models@fits)[j]
    fit_j <- models@fits[[j]]

    scores <- parallel::mclapply(seq_along(fit_j), function(g) {
      
	   # inside your mclapply over g
		if ("lgb.Booster" %in% class(fit_j[[g]])) {
		  fit_j[[g]] <- lightgbm::lgb.load(
			paste0(tmp_path_gbm, "/", taxon, "_rep", g, "_mod", j, "_", level, ".rds")
		  )
		}

      ## predict and evaluate, guard with try
	  fit <- fit_j[[g]]
      pr <- try(nsdm.prd(fit, testa[[g]]), silent = TRUE)
      if (inherits(pr, "try-error")) return(pr)

      sc <- try(nsdm.ceval2(f = pr, pa = papa[[g]]), silent = TRUE)
      return(sc)
    }, mc.cores = ncores)

    ## drop failures
    ok <- !vapply(scores, inherits, logical(1), what = "try-error")
    scores <- scores[ok]
    names(scores) <- sprintf("replicate_%02d", seq_len(length(scores)))

    lisa[[tag_j]] <- scores
  }

  out@performance <- lisa
  return(out)
}
