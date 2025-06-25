#' Preparation function: check input data, collect meta information, take care of data subsetting. Called by model fitting functions.
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (antoine.adde@eawag.ch)
#' @export
preps=function(env=parent.frame(), call, evaluationdomain=character()){

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
	evaluationdomain = env$evaluationdomain,
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

# Regional-only points for validation

if (evaluationdomain == "regionalonly") {
# 1. Tag points as regional or global
  rsts_ref_file <- file.path(w_path, "tmp", "settings", "ref_rasters.rds")
  rsts_ref <- readRDS(rsts_ref_file)
  rst_reg_gloproj <- unwrap(rsts_ref$rst_reg_gloproj)

  coords <- dat[, c("X", "Y")]
  cells_reg <- cellFromXY(rst_reg_gloproj, coords)
  in_reg <- !is.na(cells_reg)
  dat$level <- ifelse(in_reg, "reg", "glo")

# 2. Define pools
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

    # Sample validation from regional pool
    test_sample <- valid_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.3)

    # Remove validation coordinates from training pool
    val_coords <- test_sample %>% dplyr::select(X, Y)
    training_pool <- train_pool %>%
      dplyr::anti_join(val_coords, by = c("X", "Y"))

    # Sample training set from remaining pool
    train_sample <- training_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.7)

    obschoice[[i]] <- train_sample
    testing[[i]]   <- test_sample
  }

} else if (env$replicatetype == "clustered_splitsample") {
  for (i in seq_len(env$reps)) {
    set.seed(i)

    # ---- VALIDATION: Clustered sampling from regional points ----
    pres_reg <- valid_pool %>% dplyr::filter(Presence == 1)
    abs_reg  <- valid_pool %>% dplyr::filter(Presence == 0)

    n_clusters_reg <- ceiling(sqrt(nrow(pres_reg)))

    kmp_reg <- kmeans(pres_reg[, c("X", "Y")], centers = n_clusters_reg, nstart = 10)
    pres_reg$cluster_id <- kmp_reg$cluster

    abs_reg$cluster_id <- apply(abs_reg[, c("X", "Y")], 1, function(xy) {
      which.min(colSums((t(kmp_reg$centers) - xy)^2))
    })

    valid_clustered <- dplyr::bind_rows(pres_reg, abs_reg)

    sampled_clusters <- sample(unique(valid_clustered$cluster_id),
                               size = ceiling(0.3 * n_clusters_reg),
                               replace = FALSE)

    test_sample <- valid_clustered %>%
      dplyr::filter(cluster_id %in% sampled_clusters) %>%
      dplyr::select(-cluster_id)

    # ---- TRAINING: Remove validation coords and sample from remainder ----
    val_coords <- test_sample %>% dplyr::select(X, Y)
    training_pool <- train_pool %>%
      dplyr::anti_join(val_coords, by = c("X", "Y"))

    train_sample <- training_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.7)

    obschoice[[i]] <- train_sample
    testing[[i]]   <- test_sample
  }
}

# 5. Add testing data to model object if applicable
if (length(testing) > 0 && exists("out")) {
  out@tesdat <- testing
}

# 6. Return
return(list(nsdm.i = out, train = obschoice))
}

# All points for validation

if (evaluationdomain == "all") {
	
obschoice <- list()
testing <- list()

if (env$replicatetype == "none") {

obschoice[[1]] <- dat

} else if (env$replicatetype == "splitsample") {

for (i in seq_len(env$reps)) {
set.seed(i)
dat$sid <- seq_len(nrow(dat))

chc <- dat %>%
dplyr::group_by(Presence) %>%
dplyr::slice_sample(prop = 0.7)

obschoice[[i]] <- dat %>% dplyr::filter(sid %in% chc$sid) %>% dplyr::select(-sid)
testing[[i]]   <- dat %>% dplyr::filter(!sid %in% chc$sid) %>% dplyr::select(-sid)
}

} else if (env$replicatetype == "clustered_splitsample") {

pres_points <- dat %>% dplyr::filter(Presence == 1)
abs_points  <- dat %>% dplyr::filter(Presence == 0)

n_clusters <- ceiling(sqrt(nrow(pres_points)))

set.seed(42)
kmp <- kmeans(pres_points[, c("X", "Y")], centers = n_clusters, nstart = 10)
pres_points$cluster_id <- kmp$cluster

abs_points$cluster_id <- apply(abs_points[, c("X", "Y")], 1, function(xy) {
which.min(colSums((t(kmp$centers) - xy)^2))
})

dat_clustered <- dplyr::bind_rows(pres_points, abs_points) %>%
dplyr::mutate(sid = dplyr::row_number())

for (i in seq_len(env$reps)) {
set.seed(i)

sampled_clusters <- sample(unique(dat_clustered$cluster_id),
					   size = ceiling(0.7 * n_clusters),
					   replace = FALSE)

obschoice[[i]] <- dat_clustered %>%
dplyr::filter(cluster_id %in% sampled_clusters) %>%
dplyr::select(-sid, -cluster_id)

testing[[i]] <- dat_clustered %>%
dplyr::filter(!cluster_id %in% sampled_clusters) %>%
dplyr::select(-sid, -cluster_id)
}
}

if (length(testing) > 0 && exists("out")) {
out@tesdat <- testing
}
return(list(nsdm.i = out, train = obschoice))
}
}
