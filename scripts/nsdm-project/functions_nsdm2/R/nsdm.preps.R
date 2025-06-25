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