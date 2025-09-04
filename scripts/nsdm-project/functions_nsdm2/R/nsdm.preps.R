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
  env$xy,
  data.frame(sid = env$sid)
)

#
# Global level strategies
#
if(level == "glo"){
 # ============================================================
  # CASE 1: Regional-only points for validation
  # ============================================================
  if (env$evaluationdomain == "regionalonly") {

    # ------------------------------------------------------------
    # Step 1: Define training and validation pools
    # ------------------------------------------------------------
    train_pool <- dat
    valid_pool <- dat %>% dplyr::filter(stringr::str_ends(sid, "_r"))

    # ------------------------------------------------------------
    # Step 2: Initialize output lists
    # ------------------------------------------------------------
    obschoice <- list()
    testing   <- list()

    # ------------------------------------------------------------
    # Step 3: Sampling logic
    # ------------------------------------------------------------
    if (env$replicatetype == "none") {
      # No replication, all data used for training
      obschoice[[1]] <- train_pool

    } else if (env$replicatetype == "splitsample") {
      # Random split: 70% training, 30% validation (only regional)
      for (i in seq_len(env$reps)) {
        set.seed(i)

        # Validation set from regional points
        test_sample <- valid_pool %>%
          dplyr::group_by(Presence) %>%
          dplyr::slice_sample(prop = 0.3)
        val_ids <- test_sample$sid

        # Training set = all except validation IDs
        train_sample <- train_pool %>% dplyr::filter(!sid %in% val_ids)

        # Store results
        obschoice[[i]] <- train_sample %>% dplyr::select(-X, -Y) %>% as.data.frame()
        testing[[i]]   <- test_sample %>% dplyr::select(-X, -Y) %>% as.data.frame()
      }

    } else if (env$replicatetype == "clustered_splitsample") {
      # Clustered split: preserve spatial grouping
      for (i in seq_len(env$reps)) {
        set.seed(i)

        # ---- Step 3.1: Separate presences and absences ----
        pres_reg <- valid_pool %>% dplyr::filter(Presence == 1)
        abs_reg  <- valid_pool %>% dplyr::filter(Presence == 0)

        # ---- Step 3.2: Cluster presence points ----
        n_clusters_reg <- ceiling(sqrt(nrow(pres_reg)))
        kmp_reg <- kmeans(pres_reg[, c("X", "Y")], centers = n_clusters_reg, nstart = 10)
        pres_reg$cluster_id <- kmp_reg$cluster

        # ---- Step 3.3: Assign absences to nearest cluster ----
        abs_reg$cluster_id <- apply(abs_reg[, c("X", "Y")], 1, function(xy) {
          which.min(colSums((t(kmp_reg$centers) - xy)^2))
        })

        # ---- Step 3.4: Combine and sample clusters ----
        valid_clustered <- dplyr::bind_rows(pres_reg, abs_reg)
        sampled_clusters <- sample(unique(valid_clustered$cluster_id),
                                   size = ceiling(0.3 * n_clusters_reg),
                                   replace = FALSE)

        # Validation set = points in sampled clusters
        test_sample <- valid_clustered %>%
          dplyr::filter(cluster_id %in% sampled_clusters) %>%
          dplyr::select(-cluster_id)
        val_ids <- test_sample$sid

        # Training set = remaining points
        train_sample <- train_pool %>% dplyr::filter(!sid %in% val_ids)

        # Store results
        obschoice[[i]] <- train_sample %>% dplyr::select(-X, -Y)
        testing[[i]]   <- test_sample %>% dplyr::select(-X, -Y)
      }
    }

    # ------------------------------------------------------------
    # Step 4: Attach testing data to model object if applicable
    # ------------------------------------------------------------
    if (length(testing) > 0 && exists("out")) {
      out@tesdat <- testing
    }

    # ------------------------------------------------------------
    # Step 5: Return
    # ------------------------------------------------------------
    return(list(nsdm.i = out, train = obschoice))
  }


  # ============================================================
  # CASE 2: All points for validation
  # ============================================================
  if (env$evaluationdomain == "all") {

    # ------------------------------------------------------------
    # Step 1: Initialize output lists
    # ------------------------------------------------------------
    obschoice <- list()
    testing   <- list()

    # ------------------------------------------------------------
    # Step 2: Sampling logic
    # ------------------------------------------------------------
    if (env$replicatetype == "none") {
      # All points used for training
      obschoice[[1]] <- dat

    } else if (env$replicatetype == "splitsample") {
      # Random split: 70% training, 30% validation
      for (i in seq_len(env$reps)) {
        set.seed(i)

        # Training selection
        chc <- dat %>%
          dplyr::group_by(Presence) %>%
          dplyr::slice_sample(prop = 0.7)

        obschoice[[i]] <- dat %>%
          dplyr::filter(sid %in% chc$sid) %>%
          dplyr::select(-X, -Y)

        # Remaining = validation set
        testing[[i]] <- dat %>%
          dplyr::filter(!sid %in% chc$sid) %>%
          dplyr::select(-X, -Y)
      }

    } else if (env$replicatetype == "clustered_splitsample") {
      # Clustered split: preserve spatial grouping

      # ---- Step 2.1: Separate presences and absences ----
      pres_points <- dat %>% dplyr::filter(Presence == 1)
      abs_points  <- dat %>% dplyr::filter(Presence == 0)

      # ---- Step 2.2: Cluster presence points ----
      n_clusters <- ceiling(sqrt(nrow(pres_points)))
      set.seed(42)
      kmp <- kmeans(pres_points[, c("X", "Y")], centers = n_clusters, nstart = 10)
      pres_points$cluster_id <- kmp$cluster

      # ---- Step 2.3: Assign absences to nearest cluster ----
      abs_points$cluster_id <- apply(abs_points[, c("X", "Y")], 1, function(xy) {
        which.min(colSums((t(kmp$centers) - xy)^2))
      })

      dat_clustered <- dplyr::bind_rows(pres_points, abs_points)

      # ---- Step 2.4: Replicated clustered sampling ----
      for (i in seq_len(env$reps)) {
        set.seed(i)

        # Select clusters for training
        sampled_clusters <- sample(unique(dat_clustered$cluster_id),
                                   size = ceiling(0.7 * n_clusters),
                                   replace = FALSE)

        # Training set
        obschoice[[i]] <- dat_clustered %>%
          dplyr::filter(cluster_id %in% sampled_clusters) %>%
          dplyr::select(-cluster_id, -X, -Y)

        # Validation set
        testing[[i]] <- dat_clustered %>%
          dplyr::filter(!cluster_id %in% sampled_clusters) %>%
          dplyr::select(-cluster_id, -X, -Y)
      }
    }

    # ------------------------------------------------------------
    # Step 3: Attach testing data to model object if applicable
    # ------------------------------------------------------------
    if (length(testing) > 0 && exists("out")) {
      out@tesdat <- testing
    }

    # ------------------------------------------------------------
    # Step 4: Return
    # ------------------------------------------------------------
    return(list(nsdm.i = out, train = obschoice))
  }
  }
  
#
# Regional level strategies
#  
if(level != "glo"){
# Load model(s)
 mod_glo <- nsdm.loadthis(species_name = species, model_name = model_name,
                          tag = paste0(model_name, "_tune"),
                          read_path = file.path(scr_path, "outputs", "d2_models", "glo"))$model@tesdat
						  
# forced_ids_list: list of character vectors, one per replicate, with reg_*_r SIDs to force in validation
forced_ids_list <- lapply(mod_glo, function(df) {
  df %>%
    dplyr::filter(
      stringr::str_ends(sid, "_r"),
      !stringr::str_starts(sid, "glo")
    ) %>%
    dplyr::pull(sid) %>%
    unique()
})

# resampling function
# Apply the companion sampling on dat, forcing specific reg_*_r SIDs into validation per replicate
# forced_ids_list: list of character vectors, one per replicate
# env must provide env$replicatetype and env$reps
apply_companion_sampling <- function(dat, env, forced_ids_list, prop_train = 0.7, seed_clusters = 42) {

  # ===== Step 0: Checks =====
  stopifnot(is.list(forced_ids_list))
  reps <- if (!is.null(env$reps)) env$reps else length(forced_ids_list)
  if (length(forced_ids_list) != reps) stop("forced_ids_list length must equal env$reps")

  # ===== Step 1: Pools and constants =====
  train_pool <- dat
  valid_pool <- dat
  prop_val   <- 1 - prop_train
  take_cols  <- function(df) df %>% dplyr::select(-X, -Y)

  # ===== Step 2: Outputs =====
  obschoice <- vector("list", reps)
  testing   <- vector("list", reps)

  # ===== Step 3: No replication =====
  if (identical(env$replicatetype, "none")) {
    obschoice[[1]] <- take_cols(train_pool)
    if (exists("out")) out@tesdat <- list()
    return(list(nsdm.i = if (exists("out")) out else NULL, train = obschoice))
  }

  # ===== Step 4: Split sample with forced ids =====
  if (identical(env$replicatetype, "splitsample")) {
    for (i in seq_len(reps)) {
      set.seed(i)

      forced_ids <- intersect(forced_ids_list[[i]], valid_pool$sid)

      # Start validation with all forced ids
      test_seed <- valid_pool %>% dplyr::filter(sid %in% forced_ids)

      # Target by class for validation
      target_by_class <- valid_pool %>%
        dplyr::count(Presence, name = "n") %>%
        dplyr::mutate(target = ceiling(n * prop_val)) %>%
        dplyr::select(Presence, target)

      # Already taken by class
      taken_by_class <- test_seed %>% dplyr::count(Presence, name = "taken")

      # Remaining need per class
      to_get <- target_by_class %>%
        dplyr::left_join(taken_by_class, by = "Presence") %>%
        dplyr::mutate(taken = dplyr::coalesce(taken, 0L),
                      need  = pmax(target - taken, 0L)) %>%
        dplyr::select(Presence, need)

      # Sample remaining per class from non forced pool
      remain_pool <- valid_pool %>% dplyr::filter(!sid %in% forced_ids)

# make an empty template without Presence for empty returns
.empty_no_presence <- remain_pool[0, setdiff(names(remain_pool), "Presence"), drop = FALSE]

add_samples <- to_get %>%
  dplyr::group_by(Presence) %>%
  dplyr::group_modify(function(df, key) {
    n_need <- df$need[1]

    pool <- remain_pool %>% dplyr::filter(Presence == key$Presence)

    if (n_need <= 0 || nrow(pool) == 0) {
      return(.empty_no_presence)
    }

    k <- min(n_need, nrow(pool))
    pool %>%
      dplyr::slice_sample(n = k) %>%
      dplyr::select(-Presence)  # drop grouping column before returning
  }) %>%
  dplyr::ungroup()

      test_sample <- dplyr::bind_rows(test_seed, add_samples) %>%
        dplyr::distinct(sid, .keep_all = TRUE)

      val_ids      <- test_sample$sid
      train_sample <- train_pool %>% dplyr::filter(!sid %in% val_ids)

      obschoice[[i]] <- take_cols(train_sample) %>% as.data.frame()
      testing[[i]]   <- take_cols(test_sample) %>% as.data.frame()
    }

  # ===== Step 5: Clustered split with forced ids at cluster level =====
  } else if (identical(env$replicatetype, "clustered_splitsample")) {

    # Build clusters once
    pres_pts <- valid_pool %>% dplyr::filter(Presence == 1)
    abs_pts  <- valid_pool %>% dplyr::filter(Presence == 0)

    n_clusters <- max(1L, ceiling(sqrt(nrow(pres_pts))))
    set.seed(seed_clusters)
    kmp <- kmeans(pres_pts[, c("X", "Y")], centers = n_clusters, nstart = 10)
    pres_pts$cluster_id <- kmp$cluster

    abs_pts$cluster_id <- apply(abs_pts[, c("X", "Y")], 1, function(xy) {
      which.min(colSums((t(kmp$centers) - xy)^2))
    })

    clustered <- dplyr::bind_rows(pres_pts, abs_pts)

    for (i in seq_len(reps)) {
      set.seed(i)

      forced_ids <- intersect(forced_ids_list[[i]], clustered$sid)

      # Any cluster containing forced ids is placed in validation
      forced_clusters <- clustered %>%
        dplyr::filter(sid %in% forced_ids) %>%
        dplyr::pull(cluster_id) %>%
        unique()

      target_val_clusters <- ceiling(prop_val * n_clusters)
      remaining_clusters  <- setdiff(unique(clustered$cluster_id), forced_clusters)
      n_more <- max(0L, target_val_clusters - length(forced_clusters))
      sampled_extra <- if (n_more > 0) {
        sample(remaining_clusters, size = min(n_more, length(remaining_clusters)), replace = FALSE)
      } else integer(0)

      val_clusters <- unique(c(forced_clusters, sampled_extra))

      test_sample <- clustered %>%
        dplyr::filter(cluster_id %in% val_clusters) %>%
        dplyr::select(-cluster_id)

      val_ids      <- test_sample$sid
      train_sample <- train_pool %>% dplyr::filter(!sid %in% val_ids)

      obschoice[[i]] <- take_cols(train_sample)
      testing[[i]]   <- take_cols(test_sample)
    }

  } else {
    stop("env$replicatetype must be 'none', 'splitsample', or 'clustered_splitsample'")
  }

    # ------------------------------------------------------------
    # Step 3: Attach testing data to model object if applicable
    # ------------------------------------------------------------
    if (length(testing) > 0 && exists("out")) {
      out@tesdat <- testing
    }

    # ------------------------------------------------------------
    # Step 4: Return
    # ------------------------------------------------------------
    return(list(nsdm.i = out, train = obschoice))
  }

}  
}