#' nsdm.preps3
#'
#' Build spatially stratified train/test sets for regional and global pseudoabsence data.
#' Presences are split by k-means spatial clustering. Each replicate selects one cluster
#' as the test set and uses the others as training. Absences are randomly sampled to
#' balance the sets according to `ratio_abs`.
#'
#' @md
#' @param pseu_reg An `nsdm.pseudoabsences` object for the regional level (can be NULL)
#' @param pseu_glo An `nsdm.pseudoabsences` object for the global level
#' @param n_reps Integer, number of replicates (also number of spatial clusters)
#' @param ratio_abs Numeric, absences per presence, set to `0` to skip absences
#' @return A named list of length `n_reps`. Each element has `reg_train`, `reg_test`,
#' `glo_train`, `glo_test`, each an `nsdm.pseudoabsences` object or `NULL`
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.preps3 <- function(pseu_reg = NULL,
                        pseu_glo,
                        n_reps = 10,
                        ratio_abs = 1) {

  ## --- helpers ---

  .subset_nsdm_by_idx <- function(obj, rows) {
    if (is.null(obj) || length(rows) == 0) return(NULL)
    out <- obj
    out@pa    <- obj@pa[rows]
    out@years <- obj@years[rows]
    out@xy    <- obj@xy[rows, , drop = FALSE]
    out@sid   <- obj@sid[rows]
    if (is.data.frame(obj@env_vars) && nrow(obj@env_vars) > 0) {
      out@env_vars <- obj@env_vars[rows, , drop = FALSE]
    }
    out
  }

  ## --- 1. Build spatial clusters on presences ---

  # Regional
  if (!is.null(pseu_reg)) {
    reg_pres_idx <- which(pseu_reg@pa == 1)
    reg_coords   <- pseu_reg@xy[reg_pres_idx, , drop = FALSE]
    if (nrow(reg_coords) < n_reps) stop("Not enough regional presences for the requested number of clusters")
    set.seed(1)
    reg_km <- stats::kmeans(reg_coords, centers = n_reps)
    reg_clusters <- reg_km$cluster
  } else {
    reg_pres_idx <- integer(0)
    reg_clusters <- NULL
  }

  # Global
  glo_pres_idx <- which(pseu_glo@pa == 1)
  glo_coords   <- pseu_glo@xy[glo_pres_idx, , drop = FALSE]
  if (nrow(glo_coords) < n_reps) stop("Not enough global presences for the requested number of clusters")
  set.seed(1)
  glo_km <- stats::kmeans(glo_coords, centers = n_reps)
  glo_clusters <- glo_km$cluster

  ## --- 2. Create splits: each replicate = one cluster as test ---

  splits <- lapply(seq_len(n_reps), function(k) {
    list(
      reg_train_idx = if (!is.null(reg_clusters)) reg_pres_idx[reg_clusters != k] else integer(0),
      reg_test_idx  = if (!is.null(reg_clusters)) reg_pres_idx[reg_clusters == k] else integer(0),
      glo_train_idx = glo_pres_idx[glo_clusters != k],
      glo_test_idx  = glo_pres_idx[glo_clusters == k]
    )
  })
  names(splits) <- paste0("rep_", seq_len(n_reps))

  ## --- 3. Add random absences to each split ---

  reg_abs_pool_all <- if (!is.null(pseu_reg)) which(pseu_reg@pa == 0) else integer(0)
  glo_abs_pool_all <- which(pseu_glo@pa == 0)

  splits_abs <- lapply(splits, function(sp) {
    # Regional
    if (!is.null(pseu_reg)) {
      need_reg_train <- round(ratio_abs * length(sp$reg_train_idx))
      need_reg_test  <- round(ratio_abs * length(sp$reg_test_idx))
      reg_pool <- reg_abs_pool_all
      reg_train_abs_idx <- if (need_reg_train > 0) sample(reg_pool, min(need_reg_train, length(reg_pool))) else integer(0)
      reg_pool <- setdiff(reg_pool, reg_train_abs_idx)
      reg_test_abs_idx  <- if (need_reg_test  > 0) sample(reg_pool, min(need_reg_test,  length(reg_pool))) else integer(0)
    } else {
      reg_train_abs_idx <- integer(0)
      reg_test_abs_idx  <- integer(0)
    }

    # Global
    need_glo_train <- round(ratio_abs * length(sp$glo_train_idx))
    need_glo_test  <- round(ratio_abs * length(sp$glo_test_idx))
    glo_pool <- glo_abs_pool_all
    glo_train_abs_idx <- if (need_glo_train > 0) sample(glo_pool, min(need_glo_train, length(glo_pool))) else integer(0)
    glo_pool <- setdiff(glo_pool, glo_train_abs_idx)
    glo_test_abs_idx  <- if (need_glo_test  > 0) sample(glo_pool, min(need_glo_test,  length(glo_pool))) else integer(0)

    list(
      reg_train_idx = sp$reg_train_idx,
      reg_test_idx  = sp$reg_test_idx,
      glo_train_idx = sp$glo_train_idx,
      glo_test_idx  = sp$glo_test_idx,
      reg_train_abs_idx = reg_train_abs_idx,
      reg_test_abs_idx  = reg_test_abs_idx,
      glo_train_abs_idx = glo_train_abs_idx,
      glo_test_abs_idx  = glo_test_abs_idx
    )
  })

  ## --- 4. Build final datasets ---

  out <- setNames(
    lapply(seq_along(splits_abs), function(j) {
      sp <- splits_abs[[j]]
      list(
        reg_train = .subset_nsdm_by_idx(pseu_reg, c(sp$reg_train_idx, sp$reg_train_abs_idx)),
        reg_test  = .subset_nsdm_by_idx(pseu_reg, c(sp$reg_test_idx,  sp$reg_test_abs_idx)),
        glo_train = .subset_nsdm_by_idx(pseu_glo,  c(sp$glo_train_idx, sp$glo_train_abs_idx)),
        glo_test  = .subset_nsdm_by_idx(pseu_glo,  c(sp$glo_test_idx,  sp$glo_test_abs_idx))
      )
    }),
    names(splits_abs)
  )

  return(out)
}