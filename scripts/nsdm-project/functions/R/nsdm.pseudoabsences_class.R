#' Generate an S4 class to store pseudoabsences data
#'
#' Not to be called directly by the user
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.pseudoabsences<-setClass("nsdm.pseudoabsences",slots=c(meta="list", # Meta information
                                                            pa="numeric", # store presence/pseudoabsence information
															years="numeric", # store year information
                                                            env_vars="data.frame", # store extracted env variables
                                                            xy="matrix", # store coordiantes
                                                            call="call")) # conserve function call