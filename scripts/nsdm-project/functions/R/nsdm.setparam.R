#' nsdm.setparam
#'
#' Define model parameter settings
#'
#' @param model_name A character string indicating the name of the target modelling algorithm
#' @param param_grid A character string indicating the complete path where the .xlsx document containing the parameter settings is stored
#' @param covariate_names A character vector indicating covariate names
#' @param tmp_path A character string indicating the path where to store temporary outputs
#' @param ncov.esm If the ESM framework is used, a numeric value indicating the total number of covariates evaluated
#' @param comb.esm If the ESM framework is used, a numeric value indicating the number of covariates combined in each individual small model
#' @param force.esm If the ESM framework is used, a character string indicating the name of covariate(s) to be forced (typically mainGLO)
#' @param weights A numeric vector of weights to be applied
#' @param nthreads Number (numeric) of cores to be used during parallel operations
#'
#' @return List with elements of class 'multi.input' which specify models to be fitted
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.setparam <- function(model_name, param_grid, covariate_names, tmp_path, ncov.esm=NULL, comb.esm=NULL, force.esm=NULL, weights=1, nthreads){
# Import parameter grid
grid.list<-lapply(excel_sheets(param_grid), read_excel, path = param_grid)
names(grid.list)<-excel_sheets(param_grid)

# Weighting?
if(length(weights)==1){
  weighting=FALSE
} else {
  weighting=TRUE}

# Set model settings
modinp<-list() # list where all model settings will be stored

## GLM
if(model_name=="glm"){
params_glm<-data.frame(grid.list[["glm"]])

multis_glm<-list()
for(p in 1:nrow(params_glm)){
param_glm<-params_glm[p,]
form_glm<-as.formula(paste("Presence~", paste(paste0("poly(", covariate_names,",", as.numeric(param_glm), ", raw=FALSE)"), collapse="+")))
multi_glm<-nsdm.multi("glm", list(formula=form_glm, family="binomial"), tag=paste0("glm-",p), weight=weighting)
modinp<-append(modinp, multi_glm)
}
}

## GAM
if(model_name=="gam"){
mgcv_gam<<-mgcv::bam

## with penalization
if(nrow(grid.list[["gam.auto"]])>1){
  params_gam_auto<-na.omit(expand.grid(data.frame(grid.list[["gam.auto"]])))
} else {
params_gam_auto<-data.frame(grid.list[["gam.auto"]])
}

multis_gam_auto<-list()
for(p in 1:nrow(params_gam_auto)){
param_gam<-params_gam_auto[p,]
form_gam<-as.formula(paste0("Presence~ " , paste(paste0("s(", covariate_names, ",bs='", param_gam$reg.spline,"')"), collapse=" + ")))
multi_gam<-nsdm.multi("mgcv_gam", list(formula=form_gam, family="binomial", method=as.character(param_gam$method), select=FALSE, control=list(nthreads=nthreads)), tag=paste0("gam-",p), weight=weighting)
multis_gam_auto<-append(multis_gam_auto, multi_gam)
}

## pure regression splines without penalization
params_gam_fx<-data.frame(grid.list[["gam.fx"]])

multis_gam_fx<-list()
for(p in 1:nrow(params_gam_fx)){
param_gam<-params_gam_fx[p,]
form_gam<-as.formula(paste0("Presence~ " ,paste(paste0("s(",covariate_names,", k=", as.numeric(param_gam),", fx=TRUE)"),collapse=" + ")))
multi_gam<-nsdm.multi("mgcv_gam", list(formula=form_gam, family="binomial", method="REML", select=FALSE),tag=paste0("gam-",p+nrow(params_gam_auto)), weight=weighting)
multis_gam_fx<-append(multis_gam_fx, multi_gam)
}
modinp<-append(multis_gam_auto,multis_gam_fx)
}

## Maxent; Note for Maxent default weights: "upweighting of background points" (https://par.nsf.gov/servlets/purl/10079053)
if(model_name=="max"){
if(nrow(grid.list[["max"]])>1){
  params_maxent<-na.omit(expand.grid(data.frame(grid.list[["max"]])))
} else {
  params_maxent<-data.frame(grid.list[["max"]])
}

for(p in 1:nrow(params_maxent)){
param_maxent<-params_maxent[p,]
multi_maxent<-nsdm.multi("maxnet",list(classes = as.character(param_maxent$classes), addsamplestobackground=F, regmult = as.numeric(param_maxent$regmult)), tag=paste0("max-",p))
modinp<-append(modinp, multi_maxent)
}
}

## RF
if(model_name=="rf"){
if(nrow(grid.list[["rf"]])>1){
  params_rf<-na.omit(expand.grid(data.frame(grid.list[["rf"]])))
} else {
  params_rf<-data.frame(grid.list[["rf"]])
}

for(p in 1:nrow(params_rf)){
param_rf<-params_rf[p,]
#multi_rf<-nsdm.multi("ranger",list(formula=Presence~.,num.trees=param_rf$num.trees, mtry=param_rf$mtry, min.node.size = param_rf$min.node.size, num.threads=nthreads), tag=paste0("rf-",p))
form_rf<-as.formula(paste("Presence~", paste(covariate_names, collapse=" + ")))
multi_rf<-nsdm.multi("randomForest",list(formula=form_rf,ntree=as.numeric(param_rf$num.trees), mtry=floor(sqrt(length(covariate_names))), nodesize = as.numeric(param_rf$min.node.size), classwt=c("0"=min(weights), "1"=max(weights))), tag=paste0("rf-",p), weight=weighting)
modinp<-append(modinp, multi_rf)
}
}

## (light)GBM
if(model_name=="gbm"){
if(nrow(grid.list[["gbm"]])>1){
  params_gbm<-na.omit(expand.grid(data.frame(grid.list[["gbm"]])))
} else {
  params_gbm<-data.frame(grid.list[["gbm"]])
}
tmp_path_gbm<-paste0(tmp_path, "/gbm")
dir.create(tmp_path_gbm, recursive=TRUE)
tmp_path_gbm_save<-paste0(tmp_path_gbm, "/lightgbm.model")

for(p in 1:nrow(params_gbm)){
param_gbm<-params_gbm[p,]
multi_gbm<-nsdm.multi("lightgbm",list(num_leaves=2^as.numeric(param_gbm$max_depth)-1, min_data_in_leaf=as.numeric(param_gbm$min_data_in_leaf), max_depth=as.numeric(param_gbm$max_depth), objective="binary",
                                   learning_rate=as.numeric(param_gbm$learning_rate), num_iterations=as.numeric(param_gbm$num_iterations), save_name=tmp_path_gbm_save, verbose=-10),
                                   weight=weighting, tag=paste0("gbm-",p))
modinp<-append(modinp, multi_gbm)
}
}

## ESM
if(model_name=="esm"){
covs<-covariate_names[1:ncov.esm]
if(is.character(force.esm)) covs<-c(covs[-length(covs)], force.esm)
combs<-combn(covs, comb.esm)


### Generate list of all possible formulas
formulas<-list()
for (j in 1:ncol(combs)){
  formula <- as.formula(paste0("Presence ~", paste0("poly(", combs[,j],",2)", collapse = "+")))
  formulas <- c(formulas, formula)
}

### Set esm settings
for(p in 1:length(formulas)){
  form_glm<-formulas[[p]]
  multi_esm<-nsdm.multi("glm", list(formula=form_glm, family="binomial"), tag=paste0("esm-",p), weight=weighting)
  modinp<-append(modinp, multi_esm)
}
}

names(modinp)<- unlist(lapply(modinp, function(x){x@tag}))

return(modinp)
}