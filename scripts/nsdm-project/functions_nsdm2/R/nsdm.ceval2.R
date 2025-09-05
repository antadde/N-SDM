#' nsdm.ceval2
#'
#' Model evaluation core function.
#'
#' @md
#' @param f Numeric vector of predicted scores for the test data, higher means more suitable
#' @param pa Integer or logical vector of presence–absence labels, must be 0 or 1, same length as `f`
#'
#' @return Named numeric vector with `AUC`, `AUC_S`, `RMSE`, `CBI`, `Score`, `threshold`,
#'   plus all elements returned by `nsdm.all.metrics`
#' @details The threshold is chosen as the prediction value that maximizes TSS
#'   along the ranked predictions. Ties are broken by the first maximum.
#' @seealso [ROCR::prediction()], [ROCR::performance()], [ecospat::ecospat.boyce()],
#'   `nsdm.all.metrics()`, `tss()`
#' @importFrom ROCR prediction performance
#' @importFrom ecospat ecospat.boyce
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export
nsdm.ceval2 <- function(f, pa) {

  ## basic hygiene
  if (length(f) != length(pa)) stop("f and pa must have the same length")
  keep <- !is.na(f) & !is.na(pa)
  f  <- f[keep]
  pa <- as.integer(pa[keep])

  ## order by score, highest first
  ordf <- order(f, decreasing = TRUE)
  rdf  <- f[ordf]
  opa  <- pa[ordf]
  cpa  <- cumsum(opa)

  ## drop redundant tail after last presence
  nsns <- which(cpa == max(cpa))
  nsns <- nsns[-1]
  if (length(nsns)) {
    rdf <- rdf[-nsns]
    opa <- opa[-nsns]
    cpa <- cpa[-nsns]
  }

  ## scan TSS over ranks
  tsss <- apply(cbind(seq_along(cpa), cpa), 1, function(x, y, z) {
    tss(x[2], x[1] - x[2], z - x[2], y - x[1] - (z - x[2]))
  }, y = length(pa), z = sum(pa))

  tre <- rdf[which.max(tsss)]

  ## threshold dependent metrics
  bina <- ifelse(f < tre, 0, 1)
  tb <- table(factor(bina, levels = c("0", "1")),
              factor(pa,   levels = c("0", "1")))
  tdep <- unlist(nsdm.all.metrics(tb[2, 2], tb[2, 1], tb[1, 2], tb[1, 1]))

  ## Boyce index
  boy <- try(ecospat::ecospat.boyce(fit = f, obs = f[pa == 1],
                                    PEplot = FALSE, rm.duplicate = TRUE,
                                    method = "spearman"),
             silent = TRUE)
  cbi <- if (!inherits(boy, "try-error") && is.numeric(boy$cor)) boy$cor else NA_real_

  ## AUC and RMSE
  prd  <- ROCR::prediction(f, pa)
  auc  <- ROCR::performance(prd, measure = "auc")@y.values[[1]]
  rmse <- ROCR::performance(prd, measure = "rmse")@y.values[[1]]

  ## Somers scaled AUC
  aucS <- 2 * auc - 1

  ## simple consensus score
  score_vec <- na.omit(c(aucS, cbi, tdep["maxTSS"]))
  score <- if (length(score_vec)) mean(score_vec) else NA_real_

  ## return
  weg <- c(auc, aucS, rmse, cbi, score, tre, tdep)
  names(weg)[1:6] <- c("AUC", "AUC_S", "RMSE", "CBI", "Score", "threshold")
  return(unlist(weg))
}
