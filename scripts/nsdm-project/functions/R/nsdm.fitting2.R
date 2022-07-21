#' nsdm.flex2
#'
#' Model fitting (parallel version)
#'
#' @param x A nsdm.pseudoabsences object
#' @param pa Vector with presence/absence values
#' @param env_vars Data.frame with environmental predictors
#' @param taxon  A character sting of the taxon for which models are fitted
#' @param replicatetype (how) should replicates be generated? may be 'none', 'splitsample',
#' 'cv' 'block-cv'
#' @param strata a numeric vector of the same length as observations with integers separating
#' cross validation replicates (used when replicatetype='block-cv')
#' @param reps Number of replicates
#' @param mod_args List with elements of class 'multi.input' which specify models to be fitted
#' @param level Character indicating the level evaluated (e.g. glo or loc)
#' @param ncores Number of cores to be used during parallel operations
#' @param save  should the model be saved in a structured way? (not implemented yet)
#' @param project character indicating the name of the project within which the models are run
#' (later used to define saving directories)
#' @param path where to save? (not implemented yet)
#' @param tmp_path Character indicating the path where to store temporary outputs
#'
#' @return A nsdm.fit object including slots for meta info, testing data for evaluation, and model objects
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)
#' @export


nsdm.flex2<-function(x=numeric(),
                   pa=numeric(),
                   env_vars=data.frame(),
                   taxon=character(),
                   replicatetype=character(),
                   reps,
                   strata=NA,
                   mod_args=list(),
				   level,
				   ncores,
				   save=FALSE,
                   project=NA,
                   path=NA,
				   tmp_path){

 # Check supplied model types
  for(i in 1:length(mod_args)){
    if(!(mod_args[[i]]@mod%in%c("glm","gam","gbm","maxnet","randomForest","mgcv_gam","lightgbm","ranger"))){
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
  
  # parallelize over replicates
  fits<-mclapply(1:reps, function(x){

    modi=list()
    # loop over models
    for(j in 1:length(mod_args)){

      if(mod_args[[j]]@mod=="maxnet"){

        mod_args[[j]]@args$data<-lis$train[[x]][,-which(colnames(lis$train[[x]]) %in% c("Presence", "X", "Y"))]
        mod_args[[j]]@args$p<-lis$train[[x]][,"Presence"]

        modi[[j]]<-try(do.call(mod_args[[j]]@mod, mod_args[[j]]@args), TRUE)
		if(class(modi[[j]])=="try-error"){
		mod_args_bis<-mod_args[[j]]@args
		mod_args_bis$addsamplestobackground<-TRUE
		modi[[j]]<-do.call(mod_args[[j]]@mod, mod_args_bis)}


      } else if(mod_args[[j]]@mod=="lightgbm"){
	  
	  tmp_path_gbm<-paste0(tmp_path, "/gbm")
      dir.create(tmp_path_gbm, recursive=TRUE)
       
		if(mod_args[[j]]@weight){
		
		# Make weight vector
        wi=which(lis$train[[x]][,"Presence"]==1)
        wt=rep(1,length(lis$train[[x]][,"Presence"]))
        wt[wi]<-round((length(wt)-length(wi))/length(wi))
		
        mod_args[[j]]@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                             label = lis$train[[x]][,"Presence"],
											 weight=wt)
											 
		} else {
		
		       mod_args[[j]]@args$data<-lgb.Dataset(as.matrix(lis$train[[x]][,-which(colnames(lis$train[[x]])%in% c("Presence", "X", "Y"))]), 
                                             label = lis$train[[x]][,"Presence"])
		}  
		
        modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args)

		saveRDS.lgb.Booster(modi[[j]], paste0(tmp_path_gbm,"/",taxon,"_rep",x,"_mod",j,"_",level,".rds"))
        
      } else{

        mod_args[[j]]@args$data=lis$train[[x]]

        # Make weight vector
        wi=which(mod_args[[j]]@args$data$Presence==1)
        wt=rep(1,nrow(mod_args[[j]]@args$data))
        wt[wi]<-round((nrow(mod_args[[j]]@args$data)-length(wi))/length(wi))

        if(mod_args[[j]]@weight){
          mod_args[[j]]@args$weights=wt
        }

        if(mod_args[[j]]@mod=="randomForest"){
          mod_args[[j]]@args$data$Presence=as.factor(mod_args[[j]]@args$data$Presence)
        }
		
        modi[[j]]=do.call(mod_args[[j]]@mod,mod_args[[j]]@args)
      }

      names(modi)[j]=ifelse(mod_args[[j]]@tag=="",paste0("model_",j),mod_args[[j]]@tag)
	  
	  return(modi)

    }}, mc.cores=ncores)

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$nsdm.i@fits=fits

  return(lis$nsdm.i)

}
