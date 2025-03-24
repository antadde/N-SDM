#' ecospat.boyce
#'
#' Updated ecospat.boyce (ecospat package) function for calculating Boyce index as in Hirzel et al. 2006
#'
#' @param fit A vector or Raster-Layer containing the predicted suitability values
#' @param obs A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
#' @param nclass Number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
#' @param window.w width of the moving window (by default 1/10 of the suitability range)
#' @param res resolution of the moving window (by default 100 focals)
#' @param PEplot if True, plot the predicted to expected ratio along the suitability class
#' @param rm.duplicate if TRUE, the correlation exclude successive duplicated values
#' @param method correlation method used to compute the boyce index
#'
#' @author Blaise Petitpierre (blaise.petitpierre@unil.ch)
#' @export

ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                          PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman') {
  
  #### internal function calculating predicted-to-expected ratio for each class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }
  
  if (class(fit) == "RasterLayer") {
    if (class(obs) == "data.frame" || class(obs) == "matrix") {
      obs <- extract(fit, obs)
    }
    fit <- getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  if(length(nclass)==1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}
