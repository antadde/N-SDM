#' nsdm.retrieve4pred
#'
#' Prepare covariate data for prediction
#'
#' @param covstk A raster stack containing the covariate layers used for prediction
#' @param scaleparam Attribute part of a scaled matrix containing the scaling parameters to be reapplied
#' @param mask A character string for the path to .rds file containing the raster with cells to be masked (value of 0) when predicting
#' @param observational A character vector containing the names of "observational" covariates for which the same value (obsval) will be used for prediction
#' @param obsval A character string indicating the strategy used to set the common value for observational covariates ("zero", "median", or "max"imum)
#'
#' @return A data.frame object on which predictions will be applied and a companion vector containing the indices for non-na cells on the original raster stack
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.retrieve4pred <- function(covstk, scaleparam, mask=NULL, observational=NULL, obsval=NULL){
# Transform as data.frame and retrieve non-na cells indices
covdf<-raster::as.data.frame(covstk, row.names=T)
covdf<-covdf[complete.cases(covdf), ]
covdf_ix<-as.numeric(row.names(covdf))

# If a mask is provided, mask predictions
if(length(mask)>0){
m<-readRDS(mask)
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