#' nsdm.ensembleeval
#'
#' Evaluate ensemble predictions for a given species, level, and nesting strategy.
#'
#' This function loads the best model(s) for each algorithm, applies them to test data,
#' and evaluates ensemble performance across replicates using specified ensemble criteria.
#'
#' @md
#' @param x An `nsdm.pseudoabsences` object that holds the master covariates, SIDs, 
#'   and presence values.
#' @param sets A list of replicates; each replicate contains `reg_train`, `reg_test`, 
#'   `glo_train`, and `glo_test`.
#' @param level Character string, either `"glo"` or `"reg"`, selecting which SIDs to use 
#'   for testing.
#' @param model_names Character vector of modeling algorithms to ensemble 
#'   (e.g. `c("glm", "rf", "esm")`).
#' @param species_name Character string, the name of the species or taxon being evaluated.
#' @param nesting_name Character string (default = `NA`), indicating the nesting strategies evaluated 
#'   `"posthoc"` or `"covariate"`.
#' @param posthoc_nesting_name Character string (default = `NA`), indicating the posthoc_nesting strategies evaluated 
#'   `"multiply"`, `"average"`, etc.
#'
#' @return A named list of ensemble evaluation scores, averaged across replicates. 
#'
#'#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.ensembleeval <- function(sets, level, model_names, species_name, nesting_name=NA, posthoc_nesting_name=NA, scratch_path){

pred_all <- list()

### =========================================================================
### GLO
### =========================================================================
  ## ------------------------
  ## Build testing covariates and presence vectors per replicate
  ## ------------------------
  x=readRDS(file.path(scratch_path, "outputs", "d1_covsels", "glo", species_name, paste0(species_name, ".rds")))$pseu.abs_i
  
  cols <- colnames(x@env_vars)
  rep_ids <- seq_along(sets)

  testa_glo <- vector("list", length(rep_ids))   # covariates only
  papa_glo  <- vector("list", length(rep_ids))   # presence vector

  for (k in rep_ids) {
    sub_sets  <- sets[[k]][grep("glo", names(sets[[k]]))]
    test_sid  <- sub_sets[[paste0("glo", "_test")]]@sid

    i_te <- match(test_sid, x@sid)
    i_te <- i_te[!is.na(i_te)]

    df_test <- cbind(
      data.frame(Presence = x@pa[i_te]),
      x@env_vars[i_te, , drop = FALSE]
    )

    testa_glo[[k]] <- df_test[, cols, drop = FALSE]
    papa_glo[[k]]  <- df_test$Presence
  }
names(testa_glo) <- names(papa_glo) <- sprintf("replicate_%02d", rep_ids)

for (model in model_names) {
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = species_name, format = "psv",
                             read_path = file.path(scr_path, "outputs", "d3_evals", "glo"))
  
   metric_row <- eval_list[eval_list$Metric == best_met, -1, drop = FALSE]
   
   modinp_top <- names(metric_row)[which.max(as.numeric(metric_row))]

  if (model == "esm") {
    modinp_top <- names(metric_row)[as.numeric(metric_row) > best_thre_esm]
  }
  
  # Load model(s)
  mod_m <- nsdm.loadthis(species_name = species_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "glo"))$model
  
  # Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
  
    modinp_top_n <- modinp_top[n]   
 
    # Predict
    outerloop <- length(testa_glo)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
    
    for (rep_id in seq_len(outerloop)) { 
	  
      model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", rep_id, "_mod", gsub(".*-", "", modinp_top_n), "_glo.rds"))
      
      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[rep_id]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", rep_id, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_glo[[rep_id]]))
      } else {
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[rep_id]], testa_glo[[rep_id]])
      }
      
      pred[[rep_id]] <- pred_i
    }
    
	  pop[[n]] <- pred
  }
  
  if (length(pop) > 1) {
   pop <- lapply(1:outerloop, function(i) rowMeans(do.call(cbind, lapply(pop, `[[`, i)), na.rm=TRUE))
  } else {
    pop <- pop[[1]]
  }  
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
 # Rearrange pred_all so it's replicate → model → vector
n_models <- length(pred_all)
n_reps <- length(pred_all[[1]])

# Make list for each replicate
GLO_preds <- vector("list", n_reps)
names(GLO_preds) <- paste0("rep", seq_len(n_reps))

for (rep_idx in seq_len(n_reps)) {
  GLO_preds[[rep_idx]] <- lapply(pred_all, `[[`, rep_idx)
  names(GLO_preds[[rep_idx]]) <- names(pred_all)
}
 
} else {
  warning("No predictions were made.")
  GLO_preds <- NULL
}

if(level == "reg"){

  ## ------------------------
  ## Build testing covariates and presence vectors per replicate
  ## ------------------------
  x=readRDS(file.path(scratch_path, "outputs", "d1_covsels", "reg", species_name, paste0(species_name, ".rds")))
  
  cols <- colnames(x$env_vars)
  rep_ids <- seq_along(sets)

  testa_reg <- vector("list", length(rep_ids))   # covariates only
  papa_reg  <- vector("list", length(rep_ids))   # presence vector
  coord_reg <- vector("list", length(rep_ids))   # coordinates vector

  for (k in rep_ids) {
    sub_sets  <- sets[[k]][grep("reg", names(sets[[k]]))]
    test_sid  <- sub_sets[[paste0("reg", "_test")]]@sid

    i_te <- match(test_sid, x$pseu.abs_i@sid)
    i_te <- i_te[!is.na(i_te)]

    df_test <- cbind(
      data.frame(Presence = x$pseu.abs_i@pa[i_te]),
      x$env_vars[i_te, , drop = FALSE],
	  x$pseu.abs_i@xy[i_te, , drop = FALSE]
    )

    testa_reg[[k]] <- df_test[, cols, drop = FALSE]
    papa_reg[[k]]  <- df_test$Presence
	coord_reg[[k]] <- df_test[, c("X","Y"), drop = FALSE]
  }
  
names(testa_reg) <- names(coord_reg) <- names(papa_reg) <- sprintf("replicate_%02d", rep_ids)

### =========================================================================
### REG POSTHOC
### =========================================================================
if (any(c("posthoc") %in% nesting_name)) {
## Loop on target algorithms
pred_all <- list()
for (model in model_names) {
  
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = species_name, format = "psv",
                             read_path = file.path(scr_path, "outputs", "d3_evals", "reg", "posthoc"))
  
   metric_row <- eval_list[eval_list$Metric == best_met, -1, drop = FALSE]
   
   modinp_top <- names(metric_row)[which.max(as.numeric(metric_row))]

  if (model == "esm") {
    modinp_top <- names(metric_row)[as.numeric(metric_row) > best_thre_esm]
  }
  
  # Load best model
  mod_m <- nsdm.loadthis(species_name = ispi_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "reg", "posthoc"))$model
  
  # Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
    modinp_top_n <- modinp_top[n]
    
    # Predict
    outerloop <- length(testa_reg)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
	   
   for (rep_id in seq_len(outerloop)) { 
	  
    model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", rep_id, "_mod", gsub(".*-", "", modinp_top_n), "_", level, "_", "posthoc", ".rds"))

      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[rep_id]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", rep_id, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_reg[[rep_id]]))
      } else {
	    vars_list <- strsplit(mod_m@meta$env_vars, ", ")[[1]]
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[rep_id]], testa_reg[[rep_id]][, c(vars_list), drop = FALSE])
      }
      pred[[rep_id]] <- pred_i
    }
    
    pop[[n]] <- pred
  }
  
  if (length(pop) > 1) {
   pop <- lapply(1:outerloop, function(i) rowMeans(do.call(cbind, lapply(pop, `[[`, i)), na.rm=TRUE))
  } else {
    pop <- pop[[1]]
  }
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
 # Rearrange pred_all so it's replicate → model → vector
n_models <- length(pred_all)
n_reps <- length(pred_all[[1]])

# Make list for each replicate
REG_posthoc_preds <- vector("list", n_reps)
names(REG_posthoc_preds) <- paste0("rep", seq_len(n_reps))

for (rep_idx in seq_len(n_reps)) {
  REG_posthoc_preds[[rep_idx]] <- lapply(pred_all, `[[`, rep_idx)
  names(REG_posthoc_preds[[rep_idx]]) <- names(pred_all)
}
 
} else {
  warning("No predictions were made.")
  REG_posthoc_preds <- NULL
}
}

### =========================================================================
### REG COVARIATE
### =========================================================================
if (any(c("covariate") %in% nesting_name)) {
## Loop on target algorithms
pred_all <- list()
for (model in model_names) {
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = species_name, format = "psv",
                             read_path = file.path(scr_path, "outputs", "d3_evals", "reg", "covariate"))
  
   metric_row <- eval_list[eval_list$Metric == best_met, -1, drop = FALSE]
   
   modinp_top <- names(metric_row)[which.max(as.numeric(metric_row))]

  if (model == "esm") {
    modinp_top <- names(metric_row)[as.numeric(metric_row) > best_thre_esm]
  }
  
  # Load best model
  mod_m <- nsdm.loadthis(species_name = ispi_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "reg", "covariate"))$model
  
# Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
      modinp_top_n <- modinp_top[n]   
 
    # Predict
    outerloop <- length(testa_reg)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
    
    for (rep_id in seq_len(outerloop)) { 
	  
    model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", rep_id, "_mod", gsub(".*-", "", modinp_top_n), "_", level, "_", "covariate", ".rds"))

      
      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[rep_id]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[rep_id]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", rep_id, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_reg[[rep_id]]))
      } else {
	    vars_list <- strsplit(mod_m@meta$env_vars, ", ")[[1]]
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[rep_id]], testa_reg[[rep_id]][, c(vars_list), drop = FALSE])
      }
      
      pred[[rep_id]] <- pred_i
    }
    
	  pop[[n]] <- pred
  }
  
  if (length(pop) > 1) {
   pop <- lapply(1:outerloop, function(i) rowMeans(do.call(cbind, lapply(pop, `[[`, i)), na.rm=TRUE))
  } else {
    pop <- pop[[1]]
  }
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
 # Rearrange pred_all so it's replicate → model → vector
n_models <- length(pred_all)
n_reps <- length(pred_all[[1]])

# Make list for each replicate
REG_covariate_preds <- vector("list", n_reps)
names(REG_covariate_preds) <- paste0("rep", seq_len(n_reps))

for (rep_idx in seq_len(n_reps)) {
  REG_covariate_preds[[rep_idx]] <- lapply(pred_all, `[[`, rep_idx)
  names(REG_covariate_preds[[rep_idx]]) <- names(pred_all)
}
 
} else {
  warning("No predictions were made.")
  REG_covariate_preds <- NULL
}
}
}

### =========================================================================
### D- Evaluate ensemble predictions
### =========================================================================
scores_ensemble <- list()

# GLO-level
target <- GLO_preds
scores <- list()
for (z in seq_len(outerloop)) {
  z_target <- target[[z]]
  score <- nsdm.ceval2(f = rowMeans(as.data.frame(z_target), na.rm = TRUE),
                       pa = papa_glo[[z]])
  scores[[z]] <- score
}
scores_ensemble[["GLO"]] <- simplify2array(scores)
w_glo <- rowMeans(scores_ensemble[["GLO"]])[weight_metric]

if(level == "reg"){
# REG-Covariate
if ("covariate" %in% nesting_name) {
    target <- REG_covariate_preds
    scores <- list()    
    for (z in seq_len(outerloop)) {
	  z_target <- target[[z]]
      score <- nsdm.ceval2(
        f = rowMeans(as.data.frame(z_target), na.rm = TRUE),
        pa = papa_reg[[z]])
      scores[[z]] <- score
    }   
    scores_ensemble[["COV"]] <- simplify2array(scores)
  }
  
if (any(c("posthoc") %in% nesting_name)) {
# REG without nesting
    target <- REG_posthoc_preds
    scores <- list()  
    for (z in seq_len(outerloop)) {
	  z_target <- target[[z]]
      score <- nsdm.ceval2(
        f = rowMeans(as.data.frame(z_target), na.rm = TRUE),
        pa = papa_reg[[z]])
      scores[[z]] <- score
    }
    scores_ensemble[["REG"]] <- simplify2array(scores)
    w_reg <- rowMeans(scores_ensemble[["REG"]])[weight_metric]

# GLO probabilities
   	glo_out <- rast(list.files(file.path(scr_path, "outputs", "d8_ensembles", "glo", ispi_name), pattern = ".tif", full.names = TRUE))
	glo_prob <- list()
	for (rep_id in seq_len(outerloop)) {
    glo_prob[[rep_id]] <- terra::extract(glo_out, coord_reg[[rep_id]], ID=FALSE) / 100
	}   
	glo_prob2 <- lapply(glo_prob, function(eux) {
      eux[eux < 0] <- 0
      return(eux)
    })
	
# Define weight sum
w_sum <- w_reg + w_glo
	
# Posthoc nesting methods
if (any(c("multiply", "multiplyw", "average", "averagew") %in% posthoc_nesting_name)){
	scores_m <- list()
	scores_mw <- list()
	scores_a <- list()
	scores_aw <- list()

target <- REG_posthoc_preds
    for (z in seq_len(outerloop)) {
    good_ix <- which(!is.na(glo_prob2[[z]]))
	z_target <- target[[z]]
	z_target <- lapply(z_target, function(g){g[good_ix]})
	papa_reg_c <- lapply(papa_reg, function(g){g[good_ix]})
	testa_reg_c <- lapply(testa_reg, function(g){g[good_ix,]})
	z_glo_prob2 <- glo_prob2[[z]]

if (any(c("multiply") %in% posthoc_nesting_name)){  
# Multiply
    score_m <- nsdm.ceval2(
    f = sqrt(rowMeans(as.data.frame(z_target), na.rm = TRUE) * z_glo_prob2)[[1]],
    pa = papa_reg[[z]])
	scores_m[[z]] <- score_m
	}

if (any(c("multiplyw") %in% posthoc_nesting_name)){  
# Multiply weigthed
	score_mw <- nsdm.ceval2(
    f =(((rowMeans(as.data.frame(z_target), na.rm = TRUE) ^ w_reg) * (z_glo_prob2[[1]] ^ w_glo)) ^ (1 / w_sum)),
    pa = papa_reg[[z]])
	scores_mw[[z]] <- score_mw
	}

if (any(c("average") %in% posthoc_nesting_name)) {  
# Average
   score_a <- nsdm.ceval2(
   f = ((rowMeans(as.data.frame(z_target), na.rm = TRUE) + z_glo_prob2[[1]]) / 2),
   pa = papa_reg[[z]])
   scores_a[[z]] <- score_a
   }

if (any(c("averagew") %in% posthoc_nesting_name)) {  
# Average weighted
   score_aw <- nsdm.ceval2(
   f = ((w_reg * rowMeans(as.data.frame(z_target), na.rm = TRUE)) + (w_glo * z_glo_prob2[[1]])) / (w_reg + w_glo),
   pa = papa_reg[[z]])
   scores_aw[[z]] <- score_aw
   }
}
}
}
}

if (any(c("posthoc") %in% nesting_name)) {
  if (exists("scores_m")  && length(scores_m)  > 0) scores_ensemble[["MUL"]]  <- simplify2array(scores_m)
  if (exists("scores_mw") && length(scores_mw) > 0) scores_ensemble[["MULW"]] <- simplify2array(scores_mw)
  if (exists("scores_a")  && length(scores_a)  > 0) scores_ensemble[["AVG"]]  <- simplify2array(scores_a)
  if (exists("scores_aw") && length(scores_aw) > 0) scores_ensemble[["AVGW"]] <- simplify2array(scores_aw)
}

#### Return
return(scores_ensemble)
}
