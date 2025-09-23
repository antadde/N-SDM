#' nsdm.partialPlot.gbm
#'
#' Partial dependency plot for xgboost, lightgbm, or ranger object
#'
#' @author Michael Mayer, \email{mayer@consultag.ch}
#' @param obj model object
#' @param pred.data Matrix to be used in prediction (no xgb.DMatrix, no lgb.Data)
#' @param xname Name of column in \code{pred.data} according to that dependency plot is calculated
#' @param n.pt Evaluation grid size (used only if x is not discrete)
#' @param x.discrete If TRUE, the evaluation grid is set to the unique values of x
#' @param subsample Fraction of random lines in pred.data to be used in prediction
#' @param which.class Which class if objective is "multi:softprob" (value from 0 to num_class - 1)
#' @param xlab, ylab, main, type, ... Parameters passed to \code{plot}

#' @return List of prepared objects

nsdm.partialPlot.gbm <- function(obj, pred.data, xname, n.pt = 19, discrete.x = FALSE, 
                        subsample = pmin(1, n.pt * 100 / nrow(pred.data)), which.class = NULL,
                        xlab = deparse(substitute(xname)), ylab = "", type = if (discrete.x) "p" else "b",
                        main = "", rug = TRUE, ...) {
  stopifnot(dim(pred.data) >= 1)
  
  if (subsample < 1) {
    n <- nrow(pred.data)
    picked <- sample(n, trunc(subsample * n))
    pred.data <- pred.data[picked, , drop = FALSE]
  }
  xv <- pred.data[, xname]
  
  if (discrete.x) {
    x <- unique(xv)
  } else {
    x <- quantile(xv, seq(0.03, 0.97, length.out = n.pt), names = FALSE)
  }
  y <- numeric(length(x))
  
  isRanger <- inherits(obj, "ranger")
  isLm <- inherits(obj, "lm") | inherits(obj, "lmrob") | inherits(obj, "lmerMod")

  for (i in seq_along(x)) {
   pred.data[, xname] <- x[i]

    if (isRanger) {
      if (!is.null(which.class)) {
        if (obj$treetype != "Probability estimation") {
          stop("Choose probability = TRUE when fitting ranger multiclass model") 
        }
        preds <- predict(obj, pred.data)$predictions[, which.class]
      }
      else {
        preds <- predict(obj, pred.data)$predictions
      }
    } else if (isLm) {
      preds <- predict(obj, pred.data) 
    } else {
      if (!is.null(which.class)) {
        preds <- predict(obj, pred.data, reshape = TRUE)[, which.class + 1] 
      } else {
        preds <- predict(obj, pred.data)
      }
    }
    
    y[i] <- mean(preds)
  }
  
  data.frame(x = x, y = y)
}
