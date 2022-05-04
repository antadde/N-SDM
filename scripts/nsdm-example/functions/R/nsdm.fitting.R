#' nsdm.flex
#'
#' Model fitting
#'
#' @param x A nsdm.pseudoabsences object
#' @param pa Vector with presence/absence values
#' @param env_vars Data.frame with environmental predictors
#' @param taxon  A character string of the taxon for which models are fitted
#' @param replicatetype (how) should replicates be generated? may be 'none', 'splitsample',
#' 'cv' 'block-cv'
#' @param strata A numeric vector of the same length as observations with integers separating
#' cross validation replicates (used when replicatetype='block-cv')
#' @param reps Number of replicates (numeric)
#' @param mod_args List with elements of class 'multi.input' which specify models to be fitted
#' @param save  should the model be saved in a structured way? (not implemented yet)
#' @param project character indicating the name of the project within which the models are run
#' (later used to define saving directories)
#' @param path where to save? (not implemented yet)
#' @param tmp_path A character string indicating the path where to store temporary outputs
#'
#' @return A nsdm.fit object including slots for meta info, testing data for evaluation, and model objects
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)
#' @export

nsdm.flex<-function(x=numeric(),
                   pa=numeric(),
                   env_vars=data.frame(),
                   taxon=character(),
                   replicatetype=character(),
                   reps,
                   strata=NA,
                   save=FALSE,
                   project=NA,
                   path=NA,
				   mod_args=list()){

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
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    # loop over models
    for(j in 1:length(mod_args)){

      if(mod_args[[j]]@mod=="maxnet"){

        mod_args[[j]]@args$data<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
        mod_args[[j]]@args$p<-lis$train[[i]][,"Presence"]

        modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args)

      } else if(mod_args[[j]]@mod=="lightgbm"){
	  
       
		if(mod_args[[j]]@weight){
		
		# Make weight vector
        wi=which(lis$train[[i]][,"Presence"]==1)
        wt=rep(1,length(lis$train[[i]][,"Presence"]))
        wt[wi]<-round((length(wt)-length(wi))/length(wi))
		
        mod_args[[j]]@args$data<-lgb.Dataset(as.matrix(lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]), 
                                             label = lis$train[[i]][,"Presence"],
											 weight=wt)
											 
		} else {
		
		       mod_args[[j]]@args$data<-lgb.Dataset(as.matrix(lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]), 
                                             label = lis$train[[i]][,"Presence"])
		}  
		
        modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args)
        
      } else{

        mod_args[[j]]@args$data=lis$train[[i]]

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

    }

    fits[[i]]=modi

  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$nsdm.i@fits=fits


  return(lis$nsdm.i)

}
