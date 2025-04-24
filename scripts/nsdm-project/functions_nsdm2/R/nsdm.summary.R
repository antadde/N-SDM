#' nsdm.summary
#'
#' Summarize nsdm.evaluation objects
#'
#' @param object A nsdm.evaluation object to be summarized
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.summary=function(object){

  mats=list()
  for(k in 1:length(object@performance)){
  perf_k<-object@performance[[k]]
  mats[[k]]=do.call("cbind",perf_k)
  }
  
  mn=do.call(cbind, lapply(mats, rowMeans, na.rm = TRUE))
  mn=as.matrix(mn)
  colnames(mn)=names(object@performance)

  return(mn)
}
