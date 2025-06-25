#' Preparation function: check input data, collect meta information, take care of data subsetting. Called by model fitting functions.
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (antoine.adde@eawag.ch)
#' @export
preps=function(env=parent.frame(),call){

  env <- as.list(env)

  ### -------------------
  ### generate nsdm.fit obj
  ### -------------------
  out <- nsdm.fit()
  out@call <- call

  ### -------------------
  ### add meta.info
  ### -------------------
  m.i <- list(
    author = Sys.info()[["user"]],
    date = Sys.time(),
    replicatetype = env$replicatetype,
    replicates = env$reps,
    taxon = env$taxon,
    env_vars = paste(colnames(env$env_vars), collapse = ", "),
    model_tag = env$mod_tag
  )

  # Add pseudoabsence info if exists
  if (inherits(env$x, "nsdm.pseudoabsences")) {
    m.i$pseudoabsence_type <- env$x@meta$type
  }

  out@meta <- m.i

### ----------------------
### partition observations
### ----------------------

dat <- cbind(
  data.frame(Presence = env$pa),
  env$env_vars,
  env$xy
)

# 1. Tag points as regional or global
if (env$evaluationdomain == "regionalonly") {
  rsts_ref_file <- file.path(w_path, "tmp", "settings", "ref_rasters.rds")
  rsts_ref <- readRDS(rsts_ref_file)
  rst_reg_gloproj <- unwrap(rsts_ref$rst_reg_gloproj)

  coords <- dat[, c("X", "Y")]
  cells_reg <- cellFromXY(rst_reg_gloproj, coords)
  in_reg <- !is.na(cells_reg)
  dat$level <- ifelse(in_reg, "reg", "glo")
} else {
  dat$level <- "glo"
}

# 2. Define training and validation pools
train_pool <- dat
valid_pool <- dat %>% dplyr::filter(level == "reg")

# 3. Initialize output lists
obschoice <- list()
testing <- list()

# 4. Sampling logic
if (env$replicatetype == "none") {
  obschoice[[1]] <- train_pool

} else if (env$replicatetype == "splitsample") {
  for (i in seq_len(env$reps)) {
    set.seed(i)
    train_pool$sid <- seq_len(nrow(train_pool))
    valid_pool$sid <- seq_len(nrow(valid_pool))

    # Training: random sample from all points
    train_sample <- train_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.7)

    obschoice[[i]] <- train_pool %>%
      dplyr::filter(sid %in% train_sample$sid) %>%
      dplyr::select(-sid)

    # Validation: random sample only from regional points
    test_sample <- valid_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.3)

    testing[[i]] <- test_sample %>% dplyr::select(-sid)
  }

} else if (env$replicatetype == "clustered_splitsample") {
  # Training: no clustering, random sample from all data
  train_pool <- train_pool %>% dplyr::mutate(sid = dplyr::row_number())

  # Validation: clustered sampling only from regional points
  pres_reg <- valid_pool %>% dplyr::filter(Presence == 1)
  abs_reg  <- valid_pool %>% dplyr::filter(Presence == 0)

  n_clusters_reg <- ceiling(sqrt(nrow(pres_reg)))

  set.seed(42)
  kmp_reg <- kmeans(pres_reg[, c("X", "Y")], centers = n_clusters_reg, nstart = 10)
  pres_reg$cluster_id <- kmp_reg$cluster

  abs_reg$cluster_id <- apply(abs_reg[, c("X", "Y")], 1, function(xy) {
    which.min(colSums((t(kmp_reg$centers) - xy)^2))
  })

  valid_clustered <- dplyr::bind_rows(pres_reg, abs_reg) %>%
    dplyr::mutate(sid = dplyr::row_number())

  for (i in seq_len(env$reps)) {
    set.seed(i)

    # Training: random sample from all points
    obschoice[[i]] <- train_pool %>%
      dplyr::slice_sample(prop = 0.7) %>%
      dplyr::select(-sid)

    # Validation: sample clusters only from regional subset
    sampled_clusters <- sample(unique(valid_clustered$cluster_id),
                               size = ceiling(0.3 * n_clusters_reg),
                               replace = FALSE)

    testing[[i]] <- valid_clustered %>%
      dplyr::filter(cluster_id %in% sampled_clusters) %>%
      dplyr::select(-sid, -cluster_id)
  }
}

# 5. Add testing data to model object if it exists
if (length(testing) > 0 && exists("out")) {
  out@tesdat <- testing
}

# 6. Return
list(nsdm.i = out, train = obschoice)
}