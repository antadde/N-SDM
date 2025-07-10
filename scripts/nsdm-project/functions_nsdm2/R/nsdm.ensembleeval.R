#' nsdm.ensembleeval
#'
#' Evaluate ensemble predictions for a given species, level, and nesting strategy.
#'
#' This function loads the best model(s) for each algorithm, applies them to test data, 
#' generates predictions, and evaluates ensemble performance using a specified criterion.
#'
#' @param level Character. Either `"glo"` or `"reg"`, indicating the spatial modeling level.
#' @param model_names Character vector. The names of the modeling algorithms to ensemble (e.g., `c("glm", "rf", "esm")`).
#' @param species_name Character. Name of the species/taxon being evaluated.
#' @param nesting_name Character (default = NA). Name of the nesting strategy to evaluate, such as `"multiply"`, `"multiplyw"`, or `"covariate"`.
#'
#' @return A named list of ensemble evaluation scores (mean values across replicates), with one element per ensemble type (e.g., `"GLO"`, `"REG"`, `"MUL"`, `"COV"`).
#'
#' @details This function assumes the model predictions and evaluation metrics have been previously computed and stored.
#' Ensemble strategies include:
#' - GLO: global-level ensemble.
#' - REG: regional ensemble.
#' - MUL: multiply nested ensemble.
#' - MULW: multiplyw nested ensemble.
#' - COV: covariate nested ensemble.
#'
#' Evaluation is performed using `nsdm.ceval()`, and results are returned as averaged scores over replicates.
#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.ensembleeval <- function(level, model_names, species_name, nesting_name=NA){
# Retrieve predictions
### =========================================================================
### GLO
### =========================================================================
## Loop on target algorithms
pred_all <- list()
for (model in model_names) {
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = species_name,
                             read_path = file.path(scr_path, "outputs", "d3_evals", "glo"))
  
  if (!is.null(eval_list) && best_met %in% rownames(eval_list)) {
   modinp_top <- colnames(eval_list)[which.max(eval_list[best_met, ])]
  } else {
    warning("best_met is not found in eval_list for model: ", model)
    next
  }
  
  if (model == "esm") {
    modinp_top <- colnames(eval_list)[eval_list[best_met, ] > best_thre_esm]
  }
  
  # Load model(s)
  mod_m <- nsdm.loadthis(species_name = species_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "glo"))$model
  
  # Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
    modinp_top_n <- modinp_top[n]
    
    # Retrieve test and training data
    testa_GLO <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, !colnames(x) %in% "Presence", drop = FALSE]
      } else {
        warning("Presence column not found in test dataset")
        x
      }
    })
    
    papa_GLO <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, "Presence"]
      } else {
        warning("Presence column not found in training dataset")
        rep(NA, nrow(x))
      }
    })
    
    # Predict
    outerloop <- length(mod_m@tesdat)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
    
    for (k in seq_len(outerloop)) {
      model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", k, "_mod", gsub(".*-", "", modinp_top_n), "_glo.rds"))
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[k]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", k, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_GLO[[k]]))
      } else {
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[k]], testa_GLO[[k]][, !(names(testa_GLO[[k]]) %in% c("X","Y","sid"))])
      }
      
      pred[[k]] <- pred_i
    }
    
    pop_n <- do.call(cbind, pred)
    pop[[n]] <- pop_n
  }
  
  if (length(pop) > 1) {
    pop <- simplify2array(pop)
    pop <- apply(pop, c(1, 2), mean)
  } else {
    pop <- pop[[1]]
  }
  
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
  GLO_preds <- simplify2array(pred_all)
} else {
  warning("No predictions were made.")
  GLO_preds <- NULL
}

if(level == "reg"){
### =========================================================================
### REG MULTIPLY
### =========================================================================
if (any(c("multiply") %in% nesting_name)) {
## Loop on target algorithms
pred_all <- list()
for (model in mod_algo) {
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = ispi_name,
                             read_path = file.path(scr_path, "outputs", "d3_evals", "reg", nesting_name))
  
  if (!is.null(eval_list) && best_met %in% rownames(eval_list)) {
   modinp_top <- colnames(eval_list)[which.max(eval_list[best_met, ])]
  } else {
    warning("best_met is not found in eval_list for model: ", model)
    next
  }
  
  if (model == "esm") {
    modinp_top <- colnames(eval_list)[eval_list[best_met, ] > best_thre_esm]
  }
  
  # Load best model
  mod_m <- nsdm.loadthis(species_name = ispi_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "reg", nesting_name))$model
  
  # Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
    modinp_top_n <- modinp_top[n]
    
    # Retrieve test and training data
    testa_REG_multiply <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, !colnames(x) %in% "Presence", drop = FALSE]
      } else {
        warning("Presence column not found in test dataset")
        x
      }
    })
    
    papa_REG_multiply <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, "Presence"]
      } else {
        warning("Presence column not found in training dataset")
        rep(NA, nrow(x))
      }
    })
    
    # Predict
    outerloop <- length(mod_m@tesdat)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
    glo_prob <- list()
    glo_out <- rast(list.files(file.path(scr_path, "outputs", "d8_ensembles", "glo", ispi_name), 
                                  pattern = ".tif", full.names = TRUE))
    
    for (k in seq_len(outerloop)) {
	  model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", k, "_mod", gsub(".*-", "", modinp_top_n), "_", level, "_", nesting_name, ".rds"))
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[k]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", k, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_REG_multiply[[k]]))
      } else {
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[k]], testa_REG_multiply[[k]][, !(names(testa_REG_multiply[[k]]) %in% c("X","Y","sid"))])

      }
      
      pred[[k]] <- pred_i
      glo_prob[[k]] <- extract(glo_out, data.frame(mod_m@tesdat[[k]]$X, mod_m@tesdat[[k]]$Y)) / 100
    }
    
    pop_n <- do.call(cbind, pred)
    pop[[n]] <- pop_n
  }
  
  if (length(pop) > 1) {
    pop <- simplify2array(pop)
    pop <- apply(pop, c(1, 2), mean)
  } else {
    pop <- pop[[1]]
  }
  
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
  REG_multiply_preds <- simplify2array(pred_all)
} else {
  warning("No predictions were made.")
  REG_multiply_preds <- NULL
}
}

### =========================================================================
### REG COVARIATE
### =========================================================================
if (any(c("covariate") %in% nesting_name)) {
# Load covariate matrix to rescale GLO predictions
mat <- nsdm.loadthis(species_name = ispi_name, read_path = file.path(scr_path, "outputs", "d1_covsels", "reg"))$env_vars

## Loop on target algorithms
pred_all <- list()
for (model in mod_algo) {
  # Retrieve evaluation table and identify best model
  eval_list <- nsdm.loadthis(model_name = model, species_name = ispi_name,
                             read_path = file.path(scr_path, "outputs", "d3_evals", "reg", nesting_name))
  
  if (!is.null(eval_list) && best_met %in% rownames(eval_list)) {
   modinp_top <- colnames(eval_list)[which.max(eval_list[best_met, ])]
  } else {
    warning("best_met is not found in eval_list for model: ", model)
    next
  }
  
  if (model == "esm") {
    modinp_top <- colnames(eval_list)[eval_list[best_met, ] > best_thre_esm]
  }
  
  # Load best model
  mod_m <- nsdm.loadthis(species_name = ispi_name, model_name = model,
                          tag = paste0(model, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "reg", nesting_name))$model
  
  # Loop over models (one in regular cases; more for ESMs)
  pop <- list()
  for (n in seq_along(modinp_top)) {
    modinp_top_n <- modinp_top[n]
    
    # Retrieve test and training data
    testa_REG_covariate <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, !colnames(x) %in% "Presence", drop = FALSE]
      } else {
        warning("Presence column not found in test dataset")
        x
      }
    })
    
    papa_REG_covariate <- lapply(mod_m@tesdat, function(x) {
      if ("Presence" %in% colnames(x)) {
        x[, "Presence"]
      } else {
        warning("Presence column not found in training dataset")
        rep(NA, nrow(x))
      }
    })
    
    # Predict
    outerloop <- length(mod_m@tesdat)
    tmp_path_gbm <- file.path(scr_path, "tmp", "gbm")
    pred <- list()
    
    for (k in seq_len(outerloop)) {
      model_path <- file.path(tmp_path_gbm, paste0(ispi_name, "_rep", k, "_mod", gsub(".*-", "", modinp_top_n), "_", level, "_", nesting_name, ".rds"))
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "lgb.Booster")) {
        if (file.exists(model_path)) {
          mod_m@fits[[modinp_top_n]][[k]] <- lgb.load(model_path)
        } else {
          warning("LightGBM model file not found: ", model_path)
        }
      }
      
      if (inherits(mod_m@fits[[modinp_top_n]][[k]], "try-error")) {
        warning("Model fit for ", modinp_top_n, " rep ", k, " is a try-error.")
        pred_i <- rep(NA, nrow(testa_REG_covariate[[k]]))
      } else {
		pred_i <- nsdm.prd(mod_m@fits[[modinp_top_n]][[k]], testa_REG_covariate[[k]][, !(names(testa_REG_covariate[[k]]) %in% c("X","Y","sid"))])
      }
      
      pred[[k]] <- pred_i
    }
    
    pop_n <- do.call(cbind, pred)
    pop[[n]] <- pop_n
  }
  
  if (length(pop) > 1) {
    pop <- simplify2array(pop)
    pop <- apply(pop, c(1, 2), mean)
  } else {
    pop <- pop[[1]]
  }
  
  # List results
  pred_all[[model]] <- pop
}

if (length(pred_all) > 0) {
  REG_covariate_preds <- simplify2array(pred_all)
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

# GLO-level ensemble
target <- GLO_preds
scores <- list()
for (z in seq_len(outerloop)) {
  score <- nsdm.ceval(f = rowMeans(as.data.frame(target[, z, ]), na.rm = TRUE),
                       pa = papa_GLO[[z]],
                       tesdat = testa_GLO[[z]],
                       crit = eval_crit)
  scores[[z]] <- score
}
scores_ensemble[["GLO"]] <- scores
scores <- simplify2array(scores)
w_glo <- rowMeans(scores)[weight_metric]

if(level == "reg"){
  # REG-level ensemble without nesting (only possible in multiply model)
if (any(c("multiply") %in% nesting_name)) {
    target <- REG_multiply_preds
    scores <- list()
    
    for (z in seq_len(outerloop)) {
      score <- nsdm.ceval(
        f = rowMeans(as.data.frame(target[, z, ]), na.rm = TRUE),
        pa = papa_REG_multiply[[z]],
        tesdat = testa_REG_multiply[[z]],
        crit = eval_crit
      )
      scores[[z]] <- score
    }
    scores_ensemble[["REG"]] <- scores
	scores <- simplify2array(scores)
    w_reg <- rowMeans(scores)[weight_metric]
  }	

  # Covariate nested ensemble
  if ("covariate" %in% nesting_name) {
    target <- REG_covariate_preds
    scores <- list()
    
    for (z in seq_len(outerloop)) {
      score <- nsdm.ceval(
        f = rowMeans(as.data.frame(target[, z, ]), na.rm = TRUE),
        pa = papa_REG_covariate[[z]],
        tesdat = testa_REG_covariate[[z]],
        crit = eval_crit
      )
      scores[[z]] <- score
    }
    
    scores_ensemble[["COV"]] <- scores
  }

  # Multiply nested ensemble
  if ("multiply" %in% nesting_name) {
    target <- REG_multiply_preds
 
   # Preprocess GLO probabilities
    glo_prob2 <- lapply(glo_prob, function(eux) {
      eux[eux < 0] <- 0
      return(eux[, 2])
    })
	
    scores <- list()
    
    for (z in seq_len(outerloop)) {
      good_ix <- which(!is.na(glo_prob2[[z]]))
      target_z <- target[good_ix, z, ]
      glo_prob2_z <- glo_prob2[[z]][good_ix]
      papa_REG_multiply_z <- papa_REG_multiply[[z]][good_ix]
      testa_REG_multiply_z <- testa_REG_multiply[[z]][good_ix, ]
	  
      score <- nsdm.ceval(
        f = sqrt(rowMeans(as.data.frame(target_z), na.rm = TRUE) * glo_prob2_z),
        pa = papa_REG_multiply_z,
        tesdat = testa_REG_multiply_z,
        crit = eval_crit
      )
      scores[[z]] <- score
    }
    
    scores_ensemble[["MUL"]] <- scores
  
  # Multiply weighted nested ensemble
if (multiply_weighted == TRUE) {
  target <- REG_multiply_preds
  
  # Preprocess GLO probabilities
  glo_prob2 <- lapply(glo_prob, function(eux) {
    eux[eux < 0] <- 0
    return(eux[, 2])
  })
  
  scores <- list()
  
  # Define weights
  w_sum <- w_reg + w_glo
  
  for (z in seq_len(outerloop)) {
    good_ix <- which(!is.na(glo_prob2[[z]]))
    target_z <- target[good_ix, z, ]
    glo_prob2_z <- glo_prob2[[z]][good_ix]
    papa_REG_multiply_z <- papa_REG_multiply[[z]][good_ix]
    testa_REG_multiply_z <- testa_REG_multiply[[z]][good_ix, ]
    
    # Compute normalized weighted geometric mean
    target_mean <- rowMeans(as.data.frame(target_z), na.rm = TRUE)
    weighted_gmean <- (target_mean ^ w_reg) * (glo_prob2_z ^ w_glo)
    weighted_gmean <- weighted_gmean ^ (1 / w_sum)
    
    # Evaluate
    score <- nsdm.ceval(
      f = weighted_gmean,
      pa = papa_REG_multiply_z,
      tesdat = testa_REG_multiply_z,
      crit = eval_crit
    )
    
    scores[[z]] <- score
  }
  
  scores_ensemble[["MULW"]] <- scores
}
  }
}

#### Return
scores_array <- lapply(scores_ensemble, simplify2array)
scores_array <- lapply(scores_array, rowMeans)

if (exists("w_reg") && exists("w_glo")) {
  return(list(scores_array = scores_array, w_reg = w_reg, w_glo = w_glo))
} else {
return(list(scores_array = scores_array))}
}

