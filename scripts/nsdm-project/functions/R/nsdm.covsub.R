#' nsdm.covsub
#'
#' Final subsetting of covariate set after selection and ranking procedures
#'
#' @param covdata A data.frame containing covariate data
#' @param rasterdata A raster stack containing corresponding original covdata
#' @param ranks A nsdm.covselrk output
#' @param thre The target number (numeric) of covariates in a model
#' @param max.thre The maximum possible number (numeric) of covariates in a model
#' @param glo.out If mainGLO included as covariate, GLO model output
#' @param glo.xy If mainGLO included as covariate, two column matrix of for XY GLO occurences coordinates
#' @param glo.buff If mainGLO included as covariate, buffer radius used for mainGLO extraction
#'
#' @author Antoine Adde
#' @export

nsdm.covsub <- function(covdata, rasterdata, ranks, thre, max.thre, glo.out=NULL, glo.xy=NULL, glo.buff=NULL){
# set rules
if(thre>max.thre) thre<-max.thre
if(thre>nrow(cov.rk_i))thre<-nrow(cov.rk_i)
# subset
covdata<-covdata[sub('.*\\.', '', unlist(ranks["var"]))[1:thre]]
rasterdata<-rasterdata[[sub('.*\\.', '', unlist(ranks["var"]))[1:thre]]]
# Add mainGLO to covariate set if wanted
if(!is.null(glo.out)){
names(glo.out)<-"mainGLO"
rasterdata<-raster::stack(rasterdata, glo.out) # add mainGLO to covariate stack
mainGLO<-scale(raster::extract(glo.out, glo.xy, buffer=glo.buff, fun=mean))
covdata<-cbind(covdata, mainGLO) 
}
# Return results
if(is.null(glo.out)) {return(list(covdata=covdata, rasterdata=rasterdata))}
if(!is.null(glo.out)){return(list(covdata=covdata, rasterdata=rasterdata, mainGLO=mainGLO))}
}


