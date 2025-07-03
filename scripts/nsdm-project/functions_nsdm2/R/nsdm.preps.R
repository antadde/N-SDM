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

#########
##### Note ::: there is something wrong with the double 0.3 then 0.7 sampling. check.

dat <- cbind(
  data.frame(Presence = env$pa),
  env$env_vars,
  env$xy,
  data.frame(sid = env$sid)
)

# Regional-only points for validation

if (env$evaluationdomain == "regionalonly") {
# 1. Define pools
train_pool <- dat
valid_pool <- dat %>%
  dplyr::filter(
    stringr::str_starts(sid, "reg") | stringr::str_starts(sid, "NA_")
  )
  
# 2. Initialize output lists
obschoice <- list()
testing <- list()

# 3. Sampling logic
if (env$replicatetype == "none") {
  obschoice[[1]] <- train_pool
} else if (env$replicatetype == "splitsample") {
  for (i in seq_len(env$reps)) {
    set.seed(i)
 
    # 1. Sample validation set from regional points
    test_sample <- valid_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.3) %>%
      dplyr::ungroup()

    # 2. Remove validation rows using row ID to avoid ambiguity
	val_ids <- test_sample$sid
    valid_ids_df <- valid_pool
    reg_left <- valid_ids_df %>% dplyr::filter(!sid %in% val_ids)

    # 3. Combine global data with remaining regional points
    train_candidates <- dplyr::bind_rows(
      train_pool %>% dplyr::filter(!sid %in% val_ids),
      reg_left
    )
	
	# 4. Sample balanced training set
    train_sample <- train_candidates %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.7) %>%
      dplyr::ungroup()

obschoice[[i]] <- train_sample %>%
  dplyr::select(-X, -Y) %>%
  as.data.frame()

testing[[i]] <- test_sample %>%
  dplyr::select(-X, -Y) %>%
  as.data.frame()
  }

} else if (env$replicatetype == "clustered_splitsample") {
  for (i in seq_len(env$reps)) {
    set.seed(i)

    # ---- VALIDATION: Clustered sampling from regional points ----
    pres_reg <- valid_pool %>% dplyr::filter(Presence == 1)
    abs_reg  <- valid_pool %>% dplyr::filter(Presence == 0)

    n_clusters_reg <- ceiling(sqrt(nrow(valid_pool)))

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

    train_sample <- data.frame(training_pool %>%
      dplyr::group_by(Presence) %>%
      dplyr::slice_sample(prop = 0.7))

    obschoice[[i]] <- train_sample %>% dplyr::select(-X, -Y)
    testing[[i]]   <- test_sample %>% dplyr::select(-X, -Y)
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

if (env$evaluationdomain == "all") {
	
obschoice <- list()
testing <- list()

if (env$replicatetype == "none") {

obschoice[[1]] <- dat

} else if (env$replicatetype == "splitsample") {

for (i in seq_len(env$reps)) {
set.seed(i)

chc <- dat %>%
dplyr::group_by(Presence) %>%
dplyr::slice_sample(prop = 0.7)

obschoice[[i]] <- dat %>% dplyr::filter(sid %in% chc$sid) %>% dplyr::select(-X, -Y)
testing[[i]]   <- dat %>% dplyr::filter(!sid %in% chc$sid) %>% dplyr::select(-X, -Y)
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

dat_clustered <- dplyr::bind_rows(pres_points, abs_points)

for (i in seq_len(env$reps)) {
set.seed(i)

sampled_clusters <- sample(unique(dat_clustered$cluster_id),
					   size = ceiling(0.7 * n_clusters),
					   replace = FALSE)

obschoice[[i]] <- dat_clustered %>%
dplyr::filter(cluster_id %in% sampled_clusters) %>%
dplyr::select(-cluster_id, -X, -Y)

testing[[i]] <- dat_clustered %>%
dplyr::filter(!cluster_id %in% sampled_clusters) %>%
dplyr::select(-cluster_id, -X, -Y)
}
}

if (length(testing) > 0 && exists("out")) {
out@tesdat <- testing
}
return(list(nsdm.i = out, train = obschoice))
}
}
