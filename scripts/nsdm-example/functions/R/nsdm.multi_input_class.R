#' Generate an S4 class to store muti.input data
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)
#'
nsdm.multi.input<-setClass("multi.input",slots=c(mod="character", # Model function
                                            args="list", # Model function arguments
                                            tag="character", # Model set-up name
                                            weight="logical")) # Should observations be weighted?

