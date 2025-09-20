#' covsel.filter
#'
#' Apply the collinearity filtering algorithm at each target level (i=variable level; ii=category level; iii= all remainders)
#'
#' @param covdata data.frame containing covariate data (continuous values) extracted at 'pa' locations
#' @param pa numeric vector of species presences (1) and absences (0)
#' @param weights numeric vector containing the weights for each value in 'pa' (of length 'pa')
#' @param force optional character vector indicating the name(s) of the covariate(s) to be forced in the final set
#' @param corcut numeric value for the correlation coefficient threshold used for identifying collinearity
#' @param variables character vector of length ncol(covdata) containing variable-level info
#' @param categories character vector of length ncol(covdata) containing category-level info
#'
#' @return A data.frame of "non-colinear" candidate covariates
#' @author Antoine Adde (antoine.adde@unil.ch)
#' @examples
#' library(covsel)
#' covdata<-data_covsel$env_vars
#' dim(covdata)
#' covdata_filter<-covsel.filter(covdata,
#'                               pa=data_covsel$pa,
#'                               variables=data_covsel$catvar$variable,
#'                               categories=data_covsel$catvar$category)
#' dim(covdata_filter)
#' @export

covsel.filter <- function(covdata, pa, weights=NULL, force=NULL, corcut=0.7, variables=NULL, categories=NULL){
# i-variable level (if available, filtering per variable)
if(length(variables)>0){
covdata.candidates <- split.default(covdata, variables)
covdata.variable.filter<-lapply(covdata.candidates, covsel.filteralgo, pa=pa, weights=weights, force=force, corcut=0) # corcut=0 select one covariate per variable
covdata.variable.filter<-Filter(Negate(is.null), covdata.variable.filter)
covdata<-do.call("cbind", covdata.variable.filter)
names(covdata)<-gsub("^.*\\.","", names(covdata))
if(length(categories)>0) categories<-categories[match(names(covdata.variable.filter), variables)]
}

# ii-category level (if available, filtering per category)
if(length(categories)>0){
covdata.candidates <- split.default(covdata, categories)
covdata.category.filter<-lapply(covdata.candidates, covsel.filteralgo, pa=pa, weights=weights, force=force, corcut=corcut)
covdata<-do.call("cbind", covdata.category.filter)
names(covdata)<-gsub("^.*\\.","", names(covdata))
}

# iii-all remainders
covdata.filter<-covsel.filteralgo(covdata, pa, weights=weights, force=force, corcut=corcut)
names(covdata.filter)<-gsub("^.*\\.","", names(covdata.filter))

# return results
return(covdata.filter)
}


