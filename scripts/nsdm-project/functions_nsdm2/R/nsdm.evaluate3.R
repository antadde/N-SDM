#' nsdm.evaluate3
#'
#' Parallel evaluation of nsdm.fit models with a suite of model assessment metrics (parallel version)
#'
#' @param x A nsdm.fit object
#' @param tester Data.frame with testing data (only mandatory if replicatetype='none' was chosen when models were fitted)
#' @param thres A character vector of the same length as number of models chosen with custom thresholds for model evaluation. for nsdm.flex outputs the thresholds have to be labelled with the same names provided to models
#' @param crit A character string indicating which threshold criterion should be considered? Currently 'pp=op'(predicted prevalence = observed prevalence), 'maxTSS' (threshold yielding maximum TSS), and 'external' (thresholds manually supplied) are possible
#' @param prevalence_correction Logical (TRUE FALSE) should imbalanced presence/absence data be upsampled to prevalence 0.5 for model evaluation.
#' @param level A character string indicating the evaluated level (e.g. CH or EU)
#' @param ncores Number of cores to be used during parallel operations
#' @param tmp_path A character string indicating the path where to store gbm temporary outputs
#'
#' @return An object of class 'nsdm.evaluation'
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)
#' @export

nsdm.evaluate3<-function(x,tester=data.frame(),thres=numeric(),crit="pp=op",prevalence_correction=FALSE, level, ncores=ncores, tmp_path=NULL){

  ### ------------------------
  ### check tresholds
  ### ------------------------

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
    }
    if(crit!="external"){
      warning("Assuming you want external tresholds to be used - setting crit='external'!")
      crit="external"
    }
  }

  if(!(crit%in%c("pp=op","maxTSS","external"))){
    stop("Invalid threshold criterion chosen!")
  }


  ### ------------------------
  ### Check testing data and prepare for evaluation
  ### ------------------------

  if(x@meta$replicatetype=="none" && nrow(tester)==0){
    stop("External testing data must be supplied for replicatetype 'none'")
  } else if(x@meta$replicatetype%in%c("cv","block-cv","splitsample")) {

    if(prevalence_correction){
      x@tesdat=lapply(x@tesdat,function(y){
        tdpres=y[which(y$Presence==1),]
        tdabs=y[which(y$Presence==0),]
        if(nrow(tdabs)<nrow(tdpres)){
          tdabs=tdabs[sample(1:nrow(tdabs),nrow(tdpres),replace=T),]
        } else if(nrow(tdpres)<nrow(tdabs)){
          tdpres=tdpres[sample(1:nrow(tdpres),nrow(tdabs),replace=T),]
        }
        return(rbind(tdpres,tdabs))
      })
    }

    outerloop<-length(x@tesdat)
    testa<-lapply(x@tesdat,function(x){
      y<-x[,-which(colnames(x)=="Presence"),drop=FALSE]
    })
    papa<-lapply(x@tesdat,function(x){
      y<-x[,"Presence"]
    })

  } else if(x@meta$replicatetype=="none"){

    outerloop<-1
    testa<-list(tester[,-which(colnames(tester)=="Presence"),drop=FALSE])
    papa<-list(tester[,"Presence"])

  }

  ### ------------------------
  ### generate nsdm.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="evaluation")
  
   # Unique identifiers for tmp gbm files
   taxon<-x@meta$taxon
   tmp_path_gbm<-paste0(tmp_path, "/gbm")

  ### -------------------------------------------
  ### Evaluate models
  ### -------------------------------------------
  
    lisa<-list()

  # loop over model types
  for(j in 1:length(x@fits)){
  
  fit_j<-x@fits[[j]]
 
  # parallelize over replicates
  scores<-mclapply(1:length(fit_j), mc.preschedule = FALSE, function(g){
  
  # Make predictions
   if("lgb.Booster" %in% class(fit_j[[g]])){
   fit_j[[g]]<-lgb.load(paste0(tmp_path_gbm, "/", taxon,"_rep",g,"_mod",j,"_",level,".rds"))
   testa[[g]]<- testa[[g]][,-which(colnames(testa[[g]]) %in% c("X","Y"))]}
   
   pred=nsdm.prd(fit_j[[g]], testa[[g]])
   scores<-NULL
	  
  # Evaluate (with external threshold if available)
    if(length(thres)==0){
       scores<-nsdm.ceval(f=pred,
                          pa=papa[[g]],
                          tesdat=testa[[g]],
                          crit=crit)
 
      } else{

        scores<-nsdm.ceval(f=pred,
                      pa=papa[[g]],
                      tesdat=testa[[g]],
                      tre=thres[which(names(thres)==names(fit_j[[g]])[j])],
                      crit=crit)

      }

      if(scores["threshold"]==Inf){
        scores["threshold"]=.5
      }
	  
	
    return(scores)}, mc.cores=ncores)

# Remove replicate failures
scores<-scores[grep("error",lapply(scores, attributes), invert=T)]
  
# Rename replicates
names(scores)=paste0("replicate_",sprintf("%02d",1:length(scores)))

# Supply fitted replicates
lisa[[names(x@fits[j])]]<-scores

}

out@performance<-lisa

return(out)
}
