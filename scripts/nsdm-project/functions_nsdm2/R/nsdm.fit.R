#' nsdm.fit
#'
#' Model fitting
#'
#' @param x A nsdm.pseudoabsences object
#' @param sets
#' @param mod_args List of class 'multi.input' containing information on models to be fitted
#' @param level A charachter string indicating the level evaluated (e.g. glo or reg)
#' @param ncores Number of cores to be used during parallel operations
#' @param tmp_path Character indicating the path where to store temporary outputs
#'
#' @return A nsdm.fit object including slots for meta info, testing data for evaluation, and model objects
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.fit<-function(x,
                   sets,
                   mod_args,
				   level,
				   ncores,
				   tmp_path){


  # Check sets
  reps = length(sets)
  
  # loop over model types
  for(j in 1:length(mod_args)){
   
  mod_args_j=mod_args[[j]]
  
  # parallelize over replicates
 modis=list()
 modit <-mclapply(1:reps, function(rep_id){
  
  # Retrieve replicate values
sub_sets  <- sets[[rep_id]][grep(level, names(sets[[rep_id]]))]
train_sid <- sub_sets[[paste0(level, "_train")]]@sid
test_sid  <- sub_sets[[paste0(level,  "_test")]]@sid

## build train and test rows from x
i_tr <- match(train_sid, x@sid); i_tr <- i_tr[!is.na(i_tr)]
i_te <- match(test_sid,  x@sid); i_te <- i_te[!is.na(i_te)]

train <- cbind(
  data.frame(sid = x@sid[i_tr], split = "train", Presence = x@pa[i_tr]),
  x@env_vars[i_tr, , drop = FALSE]
)
test <- cbind(
  data.frame(sid = x@sid[i_te], split = "test", Presence = x@pa[i_te]),
  x@env_vars[i_te, , drop = FALSE]
)

## only presence and env vars
df_train <- train[, c("Presence", colnames(x@env_vars)), drop = FALSE]
df_test <- test[, c("Presence", colnames(x@env_vars)), drop = FALSE]


  ## Weights
  if(mod_args_j@weight){
  wi=which(df_train$Presence==1)
  wt=rep(1, nrow(df_train))
  wt[wi]<-round((length(wt)-length(wi))/length(wi))
  }
  
  ### MAX
  if(mod_args_j@mod=="maxnet"){
  mod_args_j@args$data<-train[, c(colnames(x@env_vars)), drop = FALSE]
  mod_args_j@args$p<-df_train$Presence
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
  mod_args_j@args$data<-lgb.Dataset(as.matrix(df_train[, c(colnames(x@env_vars)), drop = FALSE]), 
                                    label = df_train$Presence, weight=wt)
  } else {
  mod_args_j@args$data<-lgb.Dataset(as.matrix(df_train[, c(colnames(x@env_vars)), drop = FALSE]), 
                                    label = df_train$Presence)}  
									
  modi<-do.call(mod_args_j@mod,mod_args_j@args)
  lgb.save(modi, paste0(tmp_path_gbm,"/",taxon,"_rep",x,"_mod",j,"_",level,".rds"))}
  
  ### GLM or GAM
  if(mod_args_j@mod %in% c("glm", "mgcv_gam", "mgcv_fx", "esm")){
  if(mod_args_j@weight) mod_args_j@args$weights=wt
  mod_args_j@args$data=df_train
  modi=do.call(mod_args_j@mod, mod_args_j@args)
  }
  
  ### RF
  if(mod_args_j@mod=="randomForest"){
  if(mod_args_j@weight) mod_args_j@args$weights=wt
   mod_args_j@args$data=df_train[, c(colnames(x@env_vars)), drop = FALSE]
   mod_args_j@args$data$Presence=as.factor(df_train$Presence)
   modi=do.call(mod_args_j@mod, mod_args_j@args)}
   
   ### RF with Ranger
  if(mod_args_j@mod=="ranger"){
  if(mod_args_j@weight) mod_args_j@args$case.weights=wt
   mod_args_j@args$data=df_train[, c(colnames(x@env_vars)), drop = FALSE]
   mod_args_j@args$data$Presence=as.factor(df_train$Presence)
   modi=do.call(mod_args_j@mod, mod_args_j@args)}
   
 
  return(modi)}, mc.cores=ncores)
  
# Rename replicates
names(modit)=paste0("replicate_",sprintf("%02d",1:reps))

# Supply fitted replicates
modis[[mod_args_j@tag]]<-modit
}

# supply fitted objects
# Initiate out object
  out <- nsdm.fit()
  lis <- list(nsdm.i = out)
  lis$nsdm.i@fits=modis

return(lis$nsdm.i)
}
