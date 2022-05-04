#' nsdm.covfilter
#'
#' Colinearity filtering algorithm function
#'
#' @param covdata A data.frame containing covariate data
#' @param pa A vector of presences (1) and absences (0)
#' @param method A character string indicating the target method for covariate ranking (default=quadratic glm p-values)
#' @param weights A vector of weights (of length pa)
#' @param force A character vector indicating the name(s) of covariates to be forced in final set (e.g.: mainEU)
#' @param corcut The value (numeric) of the correlation coefficient threshold used for identifying colinearity
#'
#' @return A thinned set of "non-colinear" candidate covariates
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.covfilter <- function(covdata, pa, method, weights, force=NULL, corcut=0){

# If only one covariate in the candidate set, don't do anything 
  if(ncol(covdata)==1){
  covdata.filter<-covdata
  return(covdata.filter)}
  
# If covariates to force, eliminate colinear ones
if(is.character(force)){
force_dat<-covdata[,force]
ix<-which(abs(cor(covdata)[,force])<cor_cut)
covdata<-covdata[,ix]
}
  
# Bind covariates and PA data
  covdata<-cbind(covdata, pa=as.factor(pa))
  
# Rank candidate covariates using selected method
  if(method=="IQR"){
    covranked<-data.frame(iqr=sort(sapply(subset(covdata[covdata$pa==1,], select=-c(pa)), IQR)))
  }
  
  if(method=="IQR.M"){
  covranked<-data.frame(iqr=sort(abs(sapply(subset(covdata[covdata$pa==1,], select=-c(pa)), IQR)/ sapply(subset(covdata[covdata$pa==1,], select=-c(pa)), median))))
  }
  
  if(method=="COR.P"){
  covranked<-data.frame(cor=abs(t(cor(pa, as.matrix(subset(covdata, select=-c(pa))), method="pearson", use="complete.obs"))))
  covranked<-covranked[order(-covranked$cor), , drop = FALSE]
  }
  
  if(method=="COR.S"){
    covranked<-data.frame(cor=abs(t(cor(pa, as.matrix(subset(covdata, select=-c(pa))), method="spearman", use="complete.obs"))))
    covranked<-covranked[order(-covranked$cor), , drop = FALSE]
  }
  
  if(method=="WIL"){
    res <- sort(apply(subset(covdata, select=-c(pa)), 2, function(x){
    wilcox.test(x ~ covdata$pa)$p.value
}))
   covranked<-data.frame(wil=res, row.names=names(res))
}

  if(method=="GLM"){
    res <- sort(apply(subset(covdata, select=-c(pa)), 2, function(x){
	 min(summary(glm(covdata$pa ~ poly(x, degree=2), family="binomial", weights=weights))$coefficients[2:3,4])
}))
   covranked<-data.frame(pval=res, row.names=names(res))
}

   if(method=="MAD.M"){
    covranked<-data.frame(iqr=sort(abs(sapply(subset(covdata[covdata$pa==1,], select=-c(pa)), mad)/ sapply(subset(covdata[covdata$pa==1,], select=-c(pa)), median))))
  }
  
# Reorder covdata accordingly and compute initial correlation matrix  
  covdata.ranked <- data.frame(covdata[, row.names(covranked)])
  cor.mat<-abs(cor(covdata.ranked,use="pairwise.complete.obs"))

# Thin candidate covariate set until no pairwise correlation > corcut
    sanctuaire<-NULL
	
    if (all(cor.mat[cor.mat!=1] < corcut)) {
    sanctuaire <-c(sanctuaire,row.names(cor.mat))
    } else {
    while(length(cor.mat) > 1) {
    tmp_cm<-data.frame(cor.mat[,1])
    ix<-which(tmp_cm[,1]<=corcut)                      
    var.name<-row.names(cor.mat)[1]
    sanctuaire<-c(sanctuaire,var.name)
    cor.mat<-cor.mat[ix,ix]
    if (all(cor.mat[cor.mat!=1] < corcut)) { 
    sanctuaire <-c(sanctuaire, row.names(cor.mat))
    break
    }
    }
    }

# Return dataframe of selected covariates  
covdata.filter<-as.data.frame(covdata[,sanctuaire])

# If covariates to force, add it
if(is.character(force)){
covdata.filter<-cbind(covdata.filter, force_dat)
names(covdata.filter)[ncol(covdata.filter)]<-force
}
	
    if(length(sanctuaire)==1){
    colnames(covdata.filter)<-sanctuaire}

	return(covdata.filter)
}