#' Prediction-evaluation meta.info function: generate meta information for prediction and evaluation
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)
#'
preva.meta=function(env=parent.frame(),type=character()){

  ### ------------------------
  ### generate nsdm.evaluation and add meta info
  ### ------------------------

  m.i=list()
  m.i$author=Sys.info()[["user"]]
  m.i$date=Sys.time()

  # Generate pevaluate object
  if(type=="evaluation"){
    out<-nsdm.evaluation()
    m.i$nsdm.fit=env$x@meta
    m.i$cutoff=env$crit
  } else if(type=="pseudoabsence") {
    out<-nsdm.pseudoabsences()
  } else {
    out<-nsdm.prediction()
    m.i$nsdm.fit=env$x@meta
  }

  #add Meta info
  out@meta<-m.i

  return(out)

}
