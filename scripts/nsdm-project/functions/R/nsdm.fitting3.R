#' nsdm.flex3
#'
#' Model fitting (parallel version)
#'
#' @param x A nsdm.pseudoabsences object
#' @param pa A vector containing 0/1 presence/absence values
#' @param env_vars A data.frame containing environmental predictors
#' @param taxon  A character string of the taxon for which models are fitted
#' @param replicatetype A charachter string indicating how should replicates be generated? may be 'none', 'splitsample',
#' 'cv' 'block-cv'
#' @param strata A numeric vector of the same length as observations with integers separating cross validation replicates (used when replicatetype='block-cv')
#' @param reps Numeric number of replicates
#' @param mod_args List of class 'multi.input' containing information on models to be fitted
#' @param timer Logical (TRUE or FALSE) indicating if model computation time is to be saved
#' @param level A charachter string indicating the level evaluated (e.g. glo or loc)
#' @param ncores Number of cores to be used during parallel operations
#' @param save  should the model be saved in a structured way? (not implemented yet)
#' @param project character indicating the name of the project within which the models are run (later used to define saving directories)
#' @param path where to save? (not implemented yet)
#' @param tmp_path Character indicating the path where to store temporary outputs
#'
#' @return A nsdm.fit object including slots for meta info, testing data for evaluation, and model objects
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)
#' @export


nsdm.flex3<-function(x=numeric(),
                   pa=numeric(),
                   env_vars=data.frame(),
                   taxon=character(),
                   replicatetype=character(),
                   reps,
                   mod_args=list(),
				   timer=FALSE,
				   level,
				   ncores,
				   save=FALSE,
                   project=NA,
                   path=NA,
				   tmp_path){

 # Check supplied model types
  for(i in 1:length(mod_args)){
    if(!(mod_args[[i]]@mod%in%c("glm","gam","gbm","maxnet","randomForest","mgcv_gam","mgcv_gam_fx","lightgbm","ranger"))){
      warning(paste(mod_args[[i]]@mod,"not in focal model functions. You might run in to problems when evaluating/predicting..."))
    }
  }

  # Check if pseudo absence object is supplied
  if(class(x)=="nsdm.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
	xy=data.frame(x@xy)
  }

  # check and prepare data and output
  lis=preps(call=match.call())
  modis=list()
  time_res=list()

  # loop over model types
	
  for(j in 1:length(mod_args)){
  
  print(j)
  
  mod_args_j=mod_args[[j]]
  
  # Retrieve computation time
  ptm <- proc.time()
  
  # parallelize over replicates
  modi<-mclapply(1:reps, function(x){
  
  ## Weights
  if(mod_args_j@weight){
  wi=which(lis$train[[x]][,"Presence"]==1)
  wt=rep(1,length(lis$train[[x]][,"Presence"]))
  wt[wi]<-round((length(wt)-length(wi))/length(wi))}
  
  ### MAX
  if(mod_args_j@mod=="maxnet"){
  mod_args_j@args$data<-lis$train[[x]][,-which(colnames(lis$train[[x]]) %in% c("Presence", "X", "Y"))]
  mod_args_j@args$p<-lis$train[[x]][,"Presence"]
  modi<-try(do.call(mod_args_j@mod, mod_args_j@args), TRUE)
  if(class(modi)=="try-error"){
  mod_args_bis<-mod_args_j@args
  mod_args_bis$addsamplestobackground<-TRUE
  modi<-do.call(mod_args_j@mod, mod_args_bis)}}

  ### GBM
 if(mod_args_j@mod=="lightgbm"){
  tmp_path_gbm<-paste0(tmp_path, "/gbm")
  dir.create(tmp_path_gbm, recursive=TRUE)
  if(mod_args_j@weight){
  mod_args_j@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                       label = lis$train[[x]][,"Presence"], weight=wt)
  } else {
  mod_args_j@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                       label = lis$train[[x]][,"Presence"])}  
  modi<-do.call(mod_args_j@mod,mod_args_j@args)
  saveRDS.lgb.Booster(modi, paste0(tmp_path_gbm,"/",taxon,"_rep",x,"_mod",j,"_",level,".rds"))}
  
  ### GLM or GAM
  if(mod_args_j@mod %in% c("glm", "mgcv_gam", "mgcv_fx", "esm")){
  if(mod_args_j@weight) mod_args_j@args$weights=wt
  mod_args_j@args$data=lis$train[[x]]
  modi=do.call(mod_args_j@mod,mod_args_j@args)}
  
  ### RF
  if(mod_args_j@mod=="randomForest"){
  if(mod_args_j@weight) mod_args_j@args$weights=wt
   mod_args_j@args$data=lis$train[[x]]
   mod_args_j@args$data$Presence=as.factor(mod_args_j@args$data$Presence)
   modi=do.call(mod_args_j@mod, mod_args_j@args)}
   
  return(modi)}, mc.cores=ncores)
  
# Rename replicates
names(modi)=paste0("replicate_",sprintf("%02d",1:reps))

# Retrieve computation time
time<-c(proc.time() - ptm)
time_res[[j]]<-time

# Supply fitted replicates
modis[[mod_args_j@tag]]<-modi
}

# supply fitted objects
lis$nsdm.i@fits=modis

if(timer){ return(list(mod=lis$nsdm.i, time=time_res))
} else {
return(lis$nsdm.i)}
}
