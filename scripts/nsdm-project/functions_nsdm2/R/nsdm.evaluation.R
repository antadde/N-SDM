#' nsdm.evaluation
#'
#' Parallel evaluation of nsdm.fit models with a suite of model assessment metrics
#'
#' @param x A nsdm.fit object
#' @param level A character string indicating the evaluated level (REG or GLO)
#' @param ncores Number of cores to be used during parallel operations
#' @param tmp_path A character string indicating the path where to store gbm temporary outputs
#'
#' @return An object of class 'nsdm.evaluation'
#' @export

nsdm.evaluation<-function(x, models, sets, level, ncores=ncores, tmp_path=NULL){

  ### ------------------------
  ### Check testing data and prepare for evaluation
  ### ------------------------

cols <- colnames(x@env_vars)
rep_ids <- seq_along(sets)

testa <- vector("list", length(rep_ids))
papa  <- vector("list", length(rep_ids))

for (k in rep_ids) {
  sub_sets  <- sets[[k]][grep(level, names(sets[[k]]))]
  test_sid  <- sub_sets[[paste0(level, "_test")]]@sid

  i_te <- match(test_sid, x@sid)
  i_te <- i_te[!is.na(i_te)]

  df_test <- cbind(
    data.frame(Presence = x@pa[i_te]),
    x@env_vars[i_te, , drop = FALSE]
  )

  testa[[k]] <- df_test[, cols, drop = FALSE]  # covariates only
  papa[[k]]  <- df_test$Presence                          # response
}

names(testa) <- names(papa) <- sprintf("replicate_%02d", rep_ids)

  ### ------------------------
  ### generate nsdm.evaluation and add meta info
  ### ------------------------

   out<-preva.meta(type="evaluation")
  
   # Unique identifiers for tmp gbm files
   taxon<-species
   tmp_path_gbm<-paste0(tmp_path, "/gbm")
   
   #### TO BE RESTART HERE
   
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
   testa[[g]]<- testa[[g]][,-which(colnames(testa[[g]]) %in% c("X","Y","sid"))]}
   
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
