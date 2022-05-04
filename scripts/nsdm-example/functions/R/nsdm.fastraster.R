#' nsdm.fastraster
#'
#' Parallel loading of raster data
#'
#' @param files A character vector containing paths to individual raster layers
#' @param nsplits Number of cores (numeric) to be used during parallel loading
#'
#' @return A raster stack
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.fastraster <- function(files, nsplits=parallel::detectCores()-1){

ix<-splitIndices(length(files), nsplits)

stk<-mclapply(ix, function(x){lapply(files[x],readRDS)}, mc.cores = nsplits) 

stk2<-raster::stack(unlist(stk))

return(stk2)
}