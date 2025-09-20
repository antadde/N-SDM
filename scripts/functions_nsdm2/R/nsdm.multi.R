#' nsdm.multi
#'
#' Define settings for function that should be supplied to nsdm.flex
#'
#' @param mod A character with the name of the function to be called. E.g. "gam"
#' @param args A list with arguments to be passed to the function specified in mod
#' @param tag A character with name for model set-up
#' @return A multi.input object that efficiently stores model specifications
#' @author Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.multi=function(mod,args,tag="",step=FALSE,weight=FALSE){

  out=nsdm.multi.input()
  out@tag=tag
  out@args=args
  out@mod=mod
  out@weight=weight

  return(out)
}
