#' nsdm.flex3
#'
#' Model fitting (parallel version)
#'
#' @param x A nsdm.pseudoabsences object
#' @param replicatetype A charachter string indicating how should replicates be generated? may be 'none', 'splitsample', or 'clustered_splitsample'
#' @param evaluationdomain A charachter string indicating the spatial domain for evaluation? may be 'regionalonly', or 'all'
#' @param reps Numeric number of replicates
#' @param mod_args List of class 'multi.input' containing information on models to be fitted
#' @param level A charachter string indicating the level evaluated (e.g. glo or reg)
#' @param path where to save? (not implemented yet)
#' @param ncores Number of cores to be used during parallel operations
#' @param tmp_path Character indicating the path where to store temporary outputs
#'
#' @return A nsdm.fit object including slots for meta info, testing data for evaluation, and model objects
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (antoine.adde@eawag.ch)
#' @export


nsdm.flex3<-function(x=numeric(),
                   replicatetype=character(),
				   evaluationdomain=character(),
                   reps,
                   mod_args=list(),
                   path=NA,
				   level,
				   ncores,
				   tmp_path){

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

  # loop over model types
	
  for(j in 1:length(mod_args)){
   
  mod_args_j=mod_args[[j]]
  
  # parallelize over replicates
  modi<-mclapply(1:reps, function(x){
  
  ## Weights
  if(mod_args_j@weight){
  wi=which(lis$train[[x]][,"Presence"]==1)
  wt=rep(1, nrow(lis$train[[x]]))
  wt[wi]<-round((length(wt)-length(wi))/length(wi))
  }
  
  ### MAX
  if(mod_args_j@mod=="maxnet"){
  mod_args_j@args$data<-lis$train[[x]][,-which(colnames(lis$train[[x]]) %in% c("Presence", "X", "Y"))]
  mod_args_j@args$p<-lis$train[[x]][,"Presence"]
  modi<-try(do.call(mod_args_j@mod, mod_args_j@args), TRUE)
  if("try-error" %in% class(modi)){
  mod_args_bis<-mod_args_j@args
  mod_args_bis$addsamplestobackground<-TRUE
  modi<-do.call(mod_args_j@mod, mod_args_bis)}}

  ### GBM
 if(mod_args_j@mod=="lgb.train"){
  tmp_path_gbm<-paste0(tmp_path, "/gbm")
  dir.create(tmp_path_gbm, recursive=TRUE)
  if(mod_args_j@weight){
  mod_args_j@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                       label = lis$train[[x]][,"Presence"], weight=wt)
  } else {
  mod_args_j@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                       label = lis$train[[x]][,"Presence"])}  
  modi<-do.call(mod_args_j@mod,mod_args_j@args)
  lgb.save(modi, paste0(tmp_path_gbm,"/",taxon,"_rep",x,"_mod",j,"_",level,".rds"))}
  
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
   
   ### RF with Ranger
  if(mod_args_j@mod=="ranger"){
  if(mod_args_j@weight) mod_args_j@args$case.weights=wt
   mod_args_j@args$data=lis$train[[x]]
   mod_args_j@args$data$Presence=as.factor(mod_args_j@args$data$Presence)
   modi=do.call(mod_args_j@mod, mod_args_j@args)}
   
  return(modi)}, mc.cores=ncores)
  
# Rename replicates
names(modi)=paste0("replicate_",sprintf("%02d",1:reps))

# Supply fitted replicates
modis[[mod_args_j@tag]]<-modi
}

# supply fitted objects
lis$nsdm.i@fits=modis

return(lis$nsdm.i)
}
