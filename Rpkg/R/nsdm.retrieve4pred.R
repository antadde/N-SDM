#' nsdm.retrieve4pred
#'
#' Prepare covariate data for prediction.
#'
#' @param covstk A RasterStack containing the covariate layers used for prediction.
#' @param scaleparam A data frame or matrix containing the scaling parameters to be reapplied.
#' @param mask A character string specifying the path to an `.rds` file containing a raster where cells with a value of `0` will be masked during prediction.
#' @param observational A character vector specifying the names of "observational" covariates, which will be assigned a constant value (`obsval`) for prediction.
#' @param obsval A character string indicating the strategy for setting the common value for observational covariates. Options: `"zero"`, `"median"`, or `"max"`.
#'
#' @return A list containing:
#' \item{data}{A `data.frame` containing the covariate values for prediction.}
#' \item{valid_cells}{An index vector identifying non-NA cells in the original raster stack.}
#'
#' @author Antoine Adde (\email{antoine.adde@eawag.ch})
#' @export


nsdm.retrieve4pred <- function(covstk, scaleparam, mask=NULL, observational=NULL, obsval=NULL){
# Transform as data.frame and retrieve non-na cells indices
covdf<-as.data.frame(covstk, row.names=T)
covdf<-covdf[complete.cases(covdf), ]
covdf_ix<-as.numeric(row.names(covdf))

# If a mask is provided, mask predictions
if(length(mask)>0){
m<-rast(mask)
tbm<-which(m[]==0)
tbr<-which(covdf_ix %in% tbm)
if(length(tbr)>1){
covdf<-covdf[-tbr,]
covdf_ix<-covdf_ix[-tbr]
}
}

# If any, set observational covariates to specified value at all locations
if(length(observational)>0){
for(i in 1:length(observational)){
if(obsval=="zero") covdf[,observational[i]]<-0
if(obsval=="median") covdf[,observational[i]]<-median(covdf[,observational[i]])
if(obsval=="max") covdf[,observational[i]]<-max(covdf[,observational[i]])
}
}

# Transform covariates with desired scaling parameters
mat_sc<-list()
for(i in 1:ncol(covdf)){
col<-covdf[,i] 
col_sc<-(col - scaleparam$'scaled:center'[names(covdf[i])])/scaleparam$'scaled:scale'[names(covdf[i])]
mat_sc<-append(mat_sc, data.frame(col_sc))
names(mat_sc)[i]<-names(covdf[i])
}
covdf<-data.frame(do.call(cbind,mat_sc))

return(list(covdf=covdf, covdf_ix=covdf_ix))
}