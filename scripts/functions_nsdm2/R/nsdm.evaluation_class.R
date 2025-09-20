#' Generate an S4 class to store evaluation data
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.evaluation<-setClass("nsdm.evaluation",slots=c(meta="list", # Meta information
                                                    thres="numeric", # supply external threshold
                                                    performance="list")) # conserve function call
