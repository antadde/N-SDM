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
  dat <- cbind(data.frame(Presence = env$pa), env$env_vars, env$xy, level=env$level)

  obschoice <- list()
  testing <- list()

  if (env$replicatetype == "none") {
    obschoice[[1]] <- dat
	
} else if (env$replicatetype == "splitsample") {
  for (i in 1:env$reps) {
    dat$sid <- 1:nrow(dat)
    set.seed(i)

    # Stratified sampling: group by Presence and level
    chc <- dat %>%
      dplyr::group_by(Presence, level) %>%
      dplyr::slice_sample(prop = 0.7, replace = FALSE)

    obschoice[[i]] <- dat[chc$sid, ]
    testing[[i]]   <- dat[-chc$sid, ]

    obschoice[[i]]$sid <- NULL
    testing[[i]]$sid   <- NULL
  }
}

  if (length(testing) > 0) {
    out@tesdat <- testing
  }

  return(list(nsdm.i = out, train = obschoice))
}


