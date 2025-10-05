#' ecospat.boyce
#'
#' Updated ecospat.boyce (ecospat package) function for calculating Boyce index as in Hirzel et al. 2006
#'
#' @param fit A vector containing the predicted suitability values.
#' @param obs A vector containing the predicted suitability values of validation points.
#' @param nclass Number of classes or a vector with class thresholds.
#'               If nclass = 0, a moving window approach is used.
#' @param window.w Width of the moving window (default = 1/10 of suitability range).
#' @param res Resolution of the moving window (default = 100 focal points).
#' @param rm.duplicate Logical. If TRUE, remove successive duplicated values.
#' @param method Correlation method used for the Boyce index ('spearman' or 'pearson').
#'
#' @return A list with F.ratio (Predicted/Expected), cor (Boyce index), and HS (suitability classes).
#' @export

ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                          rm.duplicate = TRUE, method = 'spearman') {
  
  # Internal function to compute predicted-to-expected ratio
  boycei <- function(interval, obs, fit) {
    pi <- sum(obs >= interval[1] & obs <= interval[2]) / length(obs)
    ei <- sum(fit >= interval[1] & fit <= interval[2]) / length(fit)
    return(round(pi / ei, 10))
  }
   
  mini <- min(fit, obs)
  maxi <- max(fit, obs)
  
  if (length(nclass) == 1) {
    if (nclass == 0) {
      # Moving window approach
      if (window.w == "default") {
        window.w <- (max(fit) - min(fit)) / 10
      }
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w) / res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  # Avoid closed interval error
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else {
      # Fixed number of bins
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini) / nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else {
    # User-defined breaks
    vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  # Compute P/E ratios
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(!is.nan(f))
  f <- f[to.keep]
  
  # Compute correlation (Boyce index)
  if (length(f) < 2) {
    b <- NA  # Need at least 2 bins for correlation
    r <- NULL
  } else {
    r <- seq_along(f)
    if (rm.duplicate) {
      r <- which(f != c(f[-1], TRUE))
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)
  }
  
  # Calculate average suitability in each bin
  HS <- rowMeans(interval)
  if (length(nclass) == 1 && nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  # Correct "trick"
  }
  HS <- HS[to.keep]
   
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}