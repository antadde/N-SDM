#' nsdm.eval
#'
#' Parallel evaluation of fitted models with common SDM metrics.
#'
#' @md
#' @param x An `nsdm.pseudoabsences` object that holds the master covariates, SIDs, and presence values
#' @param models An `nsdm.fit` object produced by your fitting step, models are in `models@fits`
#' @param sets A list of replicates, each replicate contains `reg_train`, `reg_test`, `glo_train`, `glo_test`
#' @param level Character string, selects which SIDs to use for testing
#' @param ncores Integer, number of cores used by `parallel::mclapply`, default is `1`
#' @param tmp_path Character path for LightGBM model files, default is `tempdir()`
#'
#' @return An object of class `nsdm.evaluation` with performance per model tag and replicate
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.eval <- function(x, models, sets, level,
                            ncores = 1L, tmp_path){

## ------------------------
## Build testing covariates, presence, and row ids per replicate
## ------------------------
tmp_path_gbm <- file.path(tmp_path, "gbm")

cols <- colnames(x@env_vars)
rep_ids <- seq_along(sets)

testa   <- vector("list", length(rep_ids))
papa    <- vector("list", length(rep_ids))
row_ids <- vector("list", length(rep_ids))

lev <- substr(level, 1, 3)

for (k in rep_ids) {

  sub_sets <- sets[[k]][grep(lev, names(sets[[k]]))]
  test_sid <- sub_sets[[paste0(lev, "_test")]]@sid

  i_te <- match(test_sid, x@sid)
  i_te <- i_te[!is.na(i_te)]

  row_ids[[k]] <- i_te

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
lisa    <- list()
lisa_pr <- list()    ## for esm families

for (j in seq_along(models@fits)) {

  tag_j <- names(models@fits)[j]
  fit_j <- models@fits[[j]]

  scores <- parallel::mclapply(seq_along(fit_j), function(g) {

    ## reload booster if needed
    if ("lgb.Booster" %in% class(fit_j[[g]])) {
      fit_j[[g]] <- lightgbm::lgb.load(
        paste0(tmp_path_gbm, "/", taxon, "_rep", g, "_mod", j, "_", level, ".rds")
      )
    }

    fit <- fit_j[[g]]
    pr  <- try(nsdm.prd(fit, testa[[g]]), silent = TRUE)

    ## evaluation
    if (grepl("esm", tag_j)) {
      sc <- NA
    } else {
      sc <- try(nsdm.ceval2(f = pr, pa = papa[[g]]), silent = TRUE)
    }

    ## always return row ids to parent
    list(
      score = sc,
      pred  = pr,
      rid   = row_ids[[g]]
    )

  }, mc.cores = ncores)

  ## extract scores
  scs <- lapply(scores, function(z) z$score)

  ok <- !vapply(scs, inherits, logical(1), "try-error")
  scs <- scs[ok]
  names(scs) <- sprintf("replicate_%02d", seq_len(length(scs)))

  lisa[[tag_j]] <- scs

  ## collect esm predictions + row ids (parent side)
  if (grepl("esm", tag_j)) {

    preds <- lapply(scores, function(z) z$pred)
    rids  <- lapply(scores, function(z) z$rid)

    preds <- preds[ok]
    rids  <- rids[ok]

    lisa_pr[[tag_j]] <- list(
      pred = preds,
      rid  = rids
    )
  }
}

## ------------------------
## Pooled evaluation for esm model families
## ------------------------
esm_scores <- list()

if (length(lisa_pr) > 0) {

  for (tag_j in names(lisa_pr)) {

    pred_list <- lisa_pr[[tag_j]]$pred
    rid_list  <- lisa_pr[[tag_j]]$rid

    pooled <- list()

    for (g in seq_along(pred_list)) {

      pooled[[g]] <- data.frame(
        row_id = rid_list[[g]],
        pa     = papa[[g]],
        pred   = pred_list[[g]]
      )
    }

    pooled_df <- do.call(rbind, pooled)

    ## aggregate by row id
    uids <- unique(pooled_df$row_id)
    nids <- length(uids)

    pa_final   <- numeric(nids)
    pred_final <- numeric(nids)

    for (i in seq_len(nids)) {
      id <- uids[i]
      take <- pooled_df$row_id == id
      pa_final[i]   <- pooled_df$pa[take][1]
      pred_final[i] <- mean(pooled_df$pred[take])
    }

    sc <- nsdm.ceval2(
      f  = pred_final,
      pa = pa_final
    )

    esm_scores[[tag_j]] <- sc
  }

  ## replace esm entries in lisa with pooled results
  for (tag_j in names(esm_scores)) {
    lisa[[tag_j]] <- list(pooled = esm_scores[[tag_j]])
  }
}

## ------------------------
## Store results and return
## ------------------------
out@performance <- lisa
return(out)
}
