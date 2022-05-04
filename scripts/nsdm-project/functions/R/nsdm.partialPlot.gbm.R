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
#' @param seed Random seed used if \code{subsample < 1}

#' @return List of prepared objects

nsdm.partialPlot.gbm <- function(obj, pred.data, xname, n.pt = 19, discrete.x = FALSE, 
                        subsample = pmin(1, n.pt * 100 / nrow(pred.data)), which.class = NULL,
                        xlab = deparse(substitute(xname)), ylab = "", type = if (discrete.x) "p" else "b",
                        main = "", rug = TRUE, seed = NULL, ...) {
  stopifnot(dim(pred.data) >= 1)
  
  if (subsample < 1) {
    if (!is.null(seed)) {
      set.seed(seed)
    } 
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
  
  #plot(x, y, xlab = xlab, ylab = ylab, main = main, type = type, ...)
  data.frame(x = x, y = y)
}

# Requires h2o package 
# h2o.partialPlot <- function(obj, pred.data, xname, n.pt = 19, discrete.x = FALSE, 
#                             subsample = pmin(1, n.pt * 100 / h2o.nrow(pred.data))) {
#   if (subsample < 1) {
#     pred.data <- h2o.splitFrame(pred.data, ratios = subsample)[[1]]
#   }
#   
#   xv <- pred.data[, xname]
#   
#   if (discrete.x) {
#     x <- h2o.unique(xv)
#   } else {
#     x <- as.h2o(quantile(xv, seq(0.03, 0.97, length.out = n.pt)))
#   }
#   
#   y <- numeric(h2o.nrow(x))
#   xout <- as.data.frame(x)[, 1]
#     
#   for (i in seq_along(xout)) {
#     pred.data[, xname] <- h2o.rep_len(x[i], h2o.nrow(pred.data))
#     y[i] <- mean(as.data.frame(predict(obj, pred.data))$predict)
#   }
#   data.frame(x = xout, y = y)
# }
# 
# 
# 
# partialPlot2 <- function(fit, data, x_names, x_gridsize = c(19, 5), x_discrete = NULL, n_max = 1000, seed = NULL) {
#   stopifnot((n <- nrow(data)) >= 1, (m <- length(x_names)) %in% 1:2)
#   if (is.null(x_discrete)) {
#     x_discrete <- !vapply(data[, x_names], is.numeric, TRUE, USE.NAMES = FALSE)
#   }
#   if (!is.null(seed)) {
#       set.seed(seed)
#   } 
#   if (n_max < n) {
#     data <- data[sample(n, n_max), , drop = FALSE]
#   }
#   
#   x_grid <- setNames(vector(mode = "list", length = m), x_names)
#   
#   for (i in seq_along(x_names)) {
#     xv <- data[, x_names[i]]
#     x_grid[[i]] <- sort(if (x_discrete[i]) unique(xv) else quantile(xv, seq(0.03, 0.97, length.out = x_gridsize[i]), names = FALSE))
#   }
#   
#   full <- merge(expand.grid(x_grid), data[, setdiff(colnames(data), x_names)], by = NULL)
#   out <- aggregate(list(prediction = predict(fit, full)), by = full[, x_names, drop = FALSE], FUN = mean)
#   
#   plt <- ggplot(data = out, aes_string(x = x_names[1], y = "prediction", 
#                                        color = if (m == 2L) x_names[m],
#                                        group = if (m == 2L) x_names[m]))
#   print(plt + if (x_discrete[1]) geom_point() else geom_line())
#   out
# }

# out <- partialPlot2(fit, iris, x_names = c("Sepal.Width", "Petal.Length"))
# out <- partialPlot2(fit, iris, x_names = c("Species", "Sepal.Width"))
# out <- partialPlot2(fit, iris, x_names = "Species")
# out <- partialPlot2(fit, iris, x_names = "Sepal.Width")