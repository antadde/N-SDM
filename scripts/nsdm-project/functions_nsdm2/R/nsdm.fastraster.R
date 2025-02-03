#' nsdm.fastraster
#'
#' Parallelized loading of raster datasets using multiple cores.
#'
#' @param files A character vector containing file paths to raster layers (`.tif`).
#' @param nsplits A numeric value specifying the number of CPU cores to use for parallel loading. 
#'   Defaults to `parallel::detectCores() - 1` if NULL.
#'
#' @return A `SpatRaster`.
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.fastraster <- function(files, nsplits=parallel::detectCores()-1){

ix<-splitIndices(length(files), nsplits)

stk<-mclapply(ix, function(x){lapply(files[x], rast)}, mc.cores = nsplits) 

# Extract SpatRasters from nested lists
raster_list <- unlist(stk, recursive = TRUE)

# Stack them into a single SpatRaster
raster_stack <- rast(raster_list)

return(raster_stack)}