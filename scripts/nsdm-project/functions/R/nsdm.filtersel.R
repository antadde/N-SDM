#' nsdm.filtersel
#'
#' Apply the colinearity filtering algorithm at each target level (i=covariates with focals; ii=within each category; iii= all remainders)
#'
#' @param covdata A data.frame containing covariate data
#' @param pa A vector of presences (1) and absences (0)
#' @param datasets A character vector indicating the names of the datasets (for each candidate covariate)
#' @param varnames A character vector indicating the names of the variables. Only useful if one variable is supplied with multiple attributes (e.g. min, max, mean) and focal sizes
#' @param weights A numeric vector of weights of length pa
#' @param focals A character vector indicating which datasets include focals
#' @param method A character string indicating the target method for covariate ranking (default=quadratic glm p-values)
#' @param force A character vector indicating the name of the covariates to be forced in the final set (e.g.: mainEU)
#' @param corcut The value (numeric) of the correlation coefficient threshold for identifying colinearity
#'
#' @return A thinned set of "non-colinear" candidate covariates
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.filtersel <- function(covdata, pa, datasets, varnames=NULL, weights=NULL, focals=NULL, method, force=NULL, corcut=0){
# Split candidate covariates by datasets
covdata.candidates <- split.default(data.frame(covdata), datasets)
if(length(varnames)>0) variables.candidates<- split.default(varnames, datasets)

## Level 0 for focals, if any
if(length(focals>0)){
for(f in 1:length(focals)){
foc<-focals[f]
covdata.focals<-split.default(covdata.candidates[[foc]], variables.candidates[[foc]])
covdata.focals.filter<-lapply(covdata.focals, nsdm.covfilter, method=method, pa=pa, weights=weights, corcut=0)
covdata.candidates[[foc]]<-do.call("cbind",covdata.focals.filter)
}
}

## Level 1 for each dataset
covdata.candidates.L1<-lapply(covdata.candidates, nsdm.covfilter, method=method, pa=pa, weights=weights, corcut=corcut)

## Level 2 for all selected
covdata.candidates.L2 <-nsdm.covfilter(covdata=do.call("cbind",covdata.candidates.L1), method=method, pa=pa, force=force, weights=weights, corcut=corcut)

# return results
names(covdata.candidates.L2)<-gsub("^.*\\.","", names(covdata.candidates.L2 ))
return(covdata.candidates.L2)
}


