#' Generate an S4 class to store fitted models
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.fit<-setClass("nsdm.fit",slots=c(meta="list", # Meta information
                                      tesdat="list", # Test data subset
                                      fits="list", # Model objects
                                      call="call")) # conserve function call
