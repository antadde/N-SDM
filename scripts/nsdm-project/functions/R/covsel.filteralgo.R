#' covsel.filteralgo
#'
#' Collinearity filtering algorithm
#'
#' @param covdata data.frame containing covariate data (continuous values) extracted at 'pa' locations
#' @param pa numeric vector of species presences (1) and absences (0)
#' @param weights numeric vector containing the weights for each value in 'pa' (of length 'pa')
#' @param force optional character vector indicating the name(s) of the covariate(s) to be forced in the final set
#' @param corcut numeric value for the correlation coefficient threshold used for identifying collinearity
#'
#' @return A data.frame of "non-collinear" candidate covariates
#' @author Antoine Adde (antoine.adde@unil.ch)
#' @examples
#' library(covsel)
#' covdata<-data_covsel$env_vars
#' dim(covdata)
#' covdata_filter<-covsel.filteralgo(covdata, pa=data_covsel$pa)
#' dim(covdata_filter)
#' @export

covsel.filteralgo <- function(covdata, pa, weights=NULL, force=NULL, corcut=0.7){

# Retrieve covariate names
candidates<-names(covdata)

# If only one covariate in the candidate set, don't do anything
  if(ncol(covdata)==1){
  covdata.filter<-covdata
  return(covdata.filter)}
  
# Remove covariates with less than 10 unique points (required for embedding part), but keep forced ones
pointless10<-which(apply(covdata, 2, function(x) length(unique(x)))<10)
if(length(pointless10)>0){
pointless<-pointless10[!names(pointless10)%in% force]
if(length(pointless)>0){
print(paste0("Covariate '", names(covdata)[pointless], "' has less than 10 unique points and will be discarded"))
covdata<-covdata[,-pointless]
}
pointless_f<-pointless10[names(pointless10)%in% force]
if(length(pointless_f)>0){
print(paste0("Warning: forced covariate '", names(covdata)[pointless_f], "' has less than 10 unique points"))
}
}

# If only zero or one remaining covariate in the candidate set, stop
  if(class(covdata) != "data.frame"){
  covdata.filter<-data.frame(covdata)
  names(covdata.filter)<-candidates[-pointless10]
  return(covdata.filter)}

  if(class(covdata) == "data.frame" && ncol(covdata)==0) return(NULL)

# If covariate(s) to force, eliminate collinear ones
force_dat<-data.frame()
if(is.character(force)){
if("TRUE" %in% (force %in% names(covdata))){
force_dat<-data.frame(covdata[, force])
names(force_dat)<-force
if(length(force)>1)  ix<-unique(rownames(which(abs(cor(covdata)[,force])>corcut, arr.ind=TRUE)))
if(length(force)==1) ix<-names(which(abs(cor(covdata)[,force])>corcut))
covdata<-covdata[, !(colnames(covdata) %in% ix)]
}
}

# Bind covariates and PA data
covdata<-cbind(covdata, pa=as.factor(pa))

# Rank candidate covariates using selected method
res <- sort(apply(subset(covdata, select=-c(pa)), 2, function(x){
min(summary(suppressWarnings(glm(covdata$pa ~ poly(x, degree=2), family="binomial", weights=weights)))$coefficients[2:3,4])
}))
covranked<-data.frame(pval=res, row.names=names(res))

# Reorder covdata accordingly and compute initial correlation matrix
covdata.ranked <- data.frame(covdata[, row.names(covranked)])
cor.mat<-abs(cor(covdata.ranked,use="pairwise.complete.obs"))

# Thin candidate covariate set until no correlation > corcut
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
if(ncol(force_dat)>0) covdata.filter<-cbind(covdata.filter, force_dat)

if(length(sanctuaire)==1) colnames(covdata.filter)<-sanctuaire

return(covdata.filter)
}
