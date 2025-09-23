#' nsdm.preps2
#'
#' Build matched train and test sets for regional and global pseudoabsence data.
#' Presences are split first, with global presences inheriting the regional assignment
#' when their SID contains a regional tag. Remaining global presences are sampled
#' independently. Absences are then added to each set to reach a target ratio.
#'
#'
#' @md
#' @param pseu_reg An `nsdm.pseudoabsences` object for the regional level
#' @param pseu_glo An `nsdm.pseudoabsences` object for the global level
#' @param prop_test Numeric in (0, 1), fraction of presences assigned to test
#' @param n_reps Integer, number of replicates
#' @param ratio_abs Numeric, absences per presence, set to `0` to skip absences
#'
#' @return A named list of length `n_reps`. Each element has `reg_train`, `reg_test`,
#' `glo_train`, `glo_test`, each an `nsdm.pseudoabsences` object or `NULL`
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.preps2 <- function(pseu_reg = NULL, pseu_glo,
                        prop_test = 0.3,
                        n_reps = 10,
                        ratio_abs = 1) {

  ## ----- helpers -----

  .reg_tag_from_glo <- function(glo_sid) {
    m <- regexec("(reg_[^_]+_[^_]+(?:_[A-Za-z]+)?)", glo_sid)
    hit <- regmatches(glo_sid, m)
    sapply(hit, function(x) if (length(x)) x[1] else NA_character_)
  }

  .to_core <- function(sid) sub("_[A-Za-z]+$", "", sid)

  .subset_nsdm_by_idx <- function(obj, rows) {
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

  ## ----- 1. split presences only, mirror REG to GLO when tagged -----
  
  if (!is.null(pseu_reg)) {

  .make_presence_indices <- function(pseu_reg, pseu_glo, prop_test, n_reps) {
    reg_sid <- pseu_reg@sid
    glo_sid <- pseu_glo@sid
    reg_pa  <- pseu_reg@pa
    glo_pa  <- pseu_glo@pa

    i_reg_pres <- which(reg_pa == 1)
    i_glo_pres <- which(glo_pa == 1)

    glo_reg_tag <- .reg_tag_from_glo(glo_sid)
    has_reg_tag <- !is.na(glo_reg_tag)

    reps <- lapply(seq_len(n_reps), function(i) {
      reg_core_pres <- .to_core(reg_sid[i_reg_pres])
      u_reg_core <- unique(reg_core_pres)
      reg_core_test <- runif(length(u_reg_core)) < prop_test
      reg_map <- data.frame(core = u_reg_core, test = reg_core_test, stringsAsFactors = FALSE)

      reg_is_test <- reg_map$test[match(reg_core_pres, reg_map$core)]
      reg_test_idx  <- i_reg_pres[reg_is_test]
      reg_train_idx <- i_reg_pres[!reg_is_test]

      i_glo_pres_with_tag <- intersect(which(has_reg_tag), i_glo_pres)
      glo_core_from_reg <- .to_core(glo_reg_tag[i_glo_pres_with_tag])
      copied_flag <- reg_map$test[match(glo_core_from_reg, reg_map$core)]

      matched <- !is.na(copied_flag)
      glo_test_idx_from_copy  <- i_glo_pres_with_tag[ matched & copied_flag]
      glo_train_idx_from_copy <- i_glo_pres_with_tag[ matched & !copied_flag]

      i_glo_pres_remaining <- setdiff(i_glo_pres, c(glo_test_idx_from_copy, glo_train_idx_from_copy))
      need_test <- runif(length(i_glo_pres_remaining)) < prop_test
      glo_test_idx_indep  <- i_glo_pres_remaining[need_test]
      glo_train_idx_indep <- i_glo_pres_remaining[!need_test]

      list(
        reg_train_idx = reg_train_idx,
        reg_test_idx  = reg_test_idx,
        glo_train_idx = c(glo_train_idx_from_copy, glo_train_idx_indep),
        glo_test_idx  = c(glo_test_idx_from_copy,  glo_test_idx_indep)
      )
    })

    names(reps) <- paste0("rep_", seq_len(n_reps))
    reps
  }

  ## ----- 2. add absences to balance each set -----

  .add_balanced_absences <- function(splits_idx, pseu_reg, pseu_glo, ratio) {
    out <- splits_idx
    reps <- names(splits_idx)

    reg_abs_pool_all <- which(pseu_reg@pa == 0)
    glo_abs_pool_all <- which(pseu_glo@pa == 0)

    for (j in seq_along(reps)) {
      i <- splits_idx[[reps[j]]]

      need_reg_train <- round(ratio * length(i$reg_train_idx))
      need_reg_test  <- round(ratio * length(i$reg_test_idx))
      need_glo_train <- round(ratio * length(i$glo_train_idx))
      need_glo_test  <- round(ratio * length(i$glo_test_idx))

      reg_abs_pool <- reg_abs_pool_all
      reg_train_abs_idx <- if (need_reg_train > 0) {
        sample(reg_abs_pool, size = min(need_reg_train, length(reg_abs_pool)))
      } else integer(0)
      reg_abs_pool <- setdiff(reg_abs_pool, reg_train_abs_idx)

      reg_test_abs_idx <- if (need_reg_test > 0) {
        sample(reg_abs_pool, size = min(need_reg_test, length(reg_abs_pool)))
      } else integer(0)

      glo_abs_pool <- glo_abs_pool_all
      glo_train_abs_idx <- if (need_glo_train > 0) {
        sample(glo_abs_pool, size = min(need_glo_train, length(glo_abs_pool)))
      } else integer(0)
      glo_abs_pool <- setdiff(glo_abs_pool, glo_train_abs_idx)

      glo_test_abs_idx <- if (need_glo_test > 0) {
        sample(glo_abs_pool, size = min(need_glo_test, length(glo_abs_pool)))
      } else integer(0)

      out[[reps[j]]]$reg_train_abs_idx <- reg_train_abs_idx
      out[[reps[j]]]$reg_test_abs_idx  <- reg_test_abs_idx
      out[[reps[j]]]$glo_train_abs_idx <- glo_train_abs_idx
      out[[reps[j]]]$glo_test_abs_idx  <- glo_test_abs_idx
    }

    out
  }

  ## ----- 3. build sets for all replicates -----

  .make_sets_with_absences <- function(pseu_reg, pseu_glo, splits_idx_with_abs, rep_id) {
    i <- splits_idx_with_abs[[rep_id]]
    reg_train_rows <- c(i$reg_train_idx, i$reg_train_abs_idx)
    reg_test_rows  <- c(i$reg_test_idx,  i$reg_test_abs_idx)
    glo_train_rows <- c(i$glo_train_idx, i$glo_train_abs_idx)
    glo_test_rows  <- c(i$glo_test_idx,  i$glo_test_abs_idx)
    list(
      reg_train = if (length(reg_train_rows)) .subset_nsdm_by_idx(pseu_reg, reg_train_rows) else NULL,
      reg_test  = if (length(reg_test_rows))  .subset_nsdm_by_idx(pseu_reg, reg_test_rows)  else NULL,
      glo_train = if (length(glo_train_rows)) .subset_nsdm_by_idx(pseu_glo, glo_train_rows) else NULL,
      glo_test  = if (length(glo_test_rows))  .subset_nsdm_by_idx(pseu_glo, glo_test_rows)  else NULL
    )
  }

  .build_all_sets <- function(pseu_reg, pseu_glo, splits_idx_with_abs) {
    setNames(
      lapply(names(splits_idx_with_abs), function(rep_id) {
        .make_sets_with_absences(pseu_reg, pseu_glo, splits_idx_with_abs, rep_id)
      }),
      names(splits_idx_with_abs)
    )
  }

  ## ----- pipeline -----

  splits_idx <- .make_presence_indices(pseu_reg, pseu_glo, prop_test, n_reps)
  splits_idx_bal <- .add_balanced_absences(splits_idx, pseu_reg, pseu_glo, ratio_abs)

  all_sets <- .build_all_sets(pseu_reg, pseu_glo, splits_idx_bal)
  
 ## ----- handle glo-only situations -----
 
 } else {

  .make_presence_indices_glo <- function(pseu_glo, prop_test, n_reps) {
    glo_sid <- pseu_glo@sid
    glo_pa  <- pseu_glo@pa
    i_glo_pres <- which(glo_pa == 1)

    reps <- lapply(seq_len(n_reps), function(i) {

      # split presences at global level only
      glo_core_pres <- .to_core(glo_sid[i_glo_pres])
      u_glo_core <- unique(glo_core_pres)
      glo_core_test <- runif(length(u_glo_core)) < prop_test
      glo_map <- data.frame(core = u_glo_core, test = glo_core_test, stringsAsFactors = FALSE)

      glo_is_test <- glo_map$test[match(glo_core_pres, glo_map$core)]
      glo_test_idx  <- i_glo_pres[glo_is_test]
      glo_train_idx <- i_glo_pres[!glo_is_test]

      list(
        reg_train_idx = integer(0),
        reg_test_idx  = integer(0),
        glo_train_idx = glo_train_idx,
        glo_test_idx  = glo_test_idx
      )
    })

    names(reps) <- paste0("rep_", seq_len(n_reps))
    reps
  }

  .add_balanced_absences_glo <- function(splits_idx, pseu_glo, ratio) {
    out <- splits_idx
    reps <- names(splits_idx)

    glo_abs_pool_all <- which(pseu_glo@pa == 0)

    for (j in seq_along(reps)) {
      i <- splits_idx[[reps[j]]]

      need_glo_train <- round(ratio * length(i$glo_train_idx))
      need_glo_test  <- round(ratio * length(i$glo_test_idx))

      glo_abs_pool <- glo_abs_pool_all
      glo_train_abs_idx <- if (need_glo_train > 0) {
        sample(glo_abs_pool, size = min(need_glo_train, length(glo_abs_pool)))
      } else integer(0)
      glo_abs_pool <- setdiff(glo_abs_pool, glo_train_abs_idx)

      glo_test_abs_idx <- if (need_glo_test > 0) {
        sample(glo_abs_pool, size = min(need_glo_test, length(glo_abs_pool)))
      } else integer(0)

      out[[reps[j]]]$reg_train_abs_idx <- integer(0)
      out[[reps[j]]]$reg_test_abs_idx  <- integer(0)
      out[[reps[j]]]$glo_train_abs_idx <- glo_train_abs_idx
      out[[reps[j]]]$glo_test_abs_idx  <- glo_test_abs_idx
    }

    out
  }

  .make_sets_with_absences_glo <- function(pseu_glo, splits_idx_with_abs, rep_id) {
    i <- splits_idx_with_abs[[rep_id]]
    glo_train_rows <- c(i$glo_train_idx, i$glo_train_abs_idx)
    glo_test_rows  <- c(i$glo_test_idx,  i$glo_test_abs_idx)
    list(
      reg_train = NULL,
      reg_test  = NULL,
      glo_train = if (length(glo_train_rows)) .subset_nsdm_by_idx(pseu_glo, glo_train_rows) else NULL,
      glo_test  = if (length(glo_test_rows))  .subset_nsdm_by_idx(pseu_glo, glo_test_rows)  else NULL
    )
  }

  .build_all_sets_glo <- function(pseu_glo, splits_idx_with_abs) {
    setNames(
      lapply(names(splits_idx_with_abs), function(rep_id) {
        .make_sets_with_absences_glo(pseu_glo, splits_idx_with_abs, rep_id)
      }),
      names(splits_idx_with_abs)
    )
  }

  ## ----- pipeline -----

  splits_idx <- .make_presence_indices_glo(pseu_glo, prop_test, n_reps)
  splits_idx_bal <- .add_balanced_absences_glo(splits_idx, pseu_glo, ratio_abs)

  all_sets <- .build_all_sets_glo(pseu_glo, splits_idx_bal)
 
  }
  return(all_sets)
  }
