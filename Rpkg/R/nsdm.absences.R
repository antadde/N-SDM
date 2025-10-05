#' `nsdm.absences`
#'
#' Prepare background (pseudoabsence) points for species distribution modeling.
#'
#' @param n Integer. Number of background points to be generated.
#' @param pres `sf` object. Spatial points dataset containing species occurrence data.
#' @param taxon Character. Target species name.
#' @param type Character. Specifies the type of species data:
#'   - `"po"`: Presence-only data
#'   - `"pa"`: Presence-absence data
#' @param rst_ref Raster (`SpatRaster`). Reference raster used in subsequent analyses.
#' @param rst_background_weight Raster (`SpatRaster`) or `NULL`. Background weight raster used for weighted pseudoabsence sampling. If `NULL`, sampling is random.
#' @param set_max_npres_to_nabs Logical. If `TRUE`, limits the maximum number of presence points to the number of background points.
#'
#' @return An `nsdm.pseudoabsences` object.
#' @author Antoine Adde (\email{antoine.adde@eawag.ch})
#' @export
nsdm.absences<-function(n=10000,
                       pres,
                       taxon=character(),
					   type="po",
					   rst_ref,
					   rst_background_weight=NULL,
                       set_max_npres_to_nabs=TRUE,
					   rst_reg_gloproj=NULL,
					   level){

  ### ------------------------
  ### generate nsdm.pseudoabsences object and add meta info
  ### ------------------------
  out<-preva.meta(type="pseudoabsence")
  out@meta$type="random"
  out@meta$taxon=taxon
  call=match.call()
  out@call<-call

  ### ------------------------
  ### Prepare presences
  ### ------------------------
if (type == "pa") {
  
  if (!"pa" %in% colnames(pres)) stop("Error: 'pa' column not found in pres dataset.")
 
  set_max_npres_to_nabs=FALSE

  abs <- pres[pres$pa == 0, ]
  pres <- pres[pres$pa == 1, ]
}
	
# Set max_npres_to_nabs if requested
if (set_max_npres_to_nabs && nrow(pres) > n) {
target_n <- n + 100
pres <- pres[sample(seq_len(nrow(pres)), target_n), ]
}

### ------------------------
### Prepare background points
### ------------------------
if (type == "po") {
  
  if (is.null(rst_background_weight)) {
    cat("Performing random background sampling (no background weight layer provided).", "\n")
    
    # Convert raster values to a data.table
    dt_ref <- data.table(cell = 1:ncell(rst_ref), weight = values(rst_ref)[,1])  # Rename weight column
    dt_ref <- dt_ref[!is.na(weight)]  # Remove NA values
    
    # Check if dt_ref is empty
    if (nrow(dt_ref) == 0) stop("Error: No valid background points available for sampling.")
    
    # Sample randomly
    sampled_cells <- dt_ref[sample(.N, size = n * 1.5, replace = TRUE)]
  } else {
    cat("Performing weighted background sampling using the background weight layer.", "\n")
    
    # Using the background weight layer
    rst_bck <- mask(rast(rst_background_weight), rst_ref)
    dt_weights <- data.table(cell = 1:ncell(rst_bck), weight = values(rst_bck)[,1])  # Rename weight column
    dt_weights <- dt_weights[!is.na(weight)]  # Remove NA values

    # Check if dt_weights is empty
    if (nrow(dt_weights) == 0) stop("Error: No valid weighted background points available for sampling.")
    
    # Sample based on weights
    sampled_cells <- dt_weights[sample(.N, size = n * 1.5, prob = dt_weights$weight, replace = TRUE)]
  }

  # Convert sampled cells to coordinates
  coords <- xyFromCell(rst_ref, sampled_cells$cell)

  # Create an sf object
  abs <- st_as_sf(data.table(coords), coords = c("x", "y"), crs = crs(rst_ref))

  # B- Subsample to n+200 (some points will be dropped after NA covariate cleaning)
  if (nrow(abs) > n + 200) {
    abs <- abs[sample(1:nrow(abs), n + 200),]
  }
}

  ### ------------------------
  ### Prepare output
  ### ------------------------
  out@pa<-c(rep(1,nrow(pres)),
            rep(0,nrow(abs)))
  out@years<-as.numeric(c(pres$year, rep(NA, nrow(abs))))
  out@env_vars=data.frame()
  out@xy=as.matrix(rbind(data.frame(X=st_coordinates(pres)[,1], Y=st_coordinates(pres)[,2]), data.frame(X=st_coordinates(abs)[,1], Y=st_coordinates(abs)[,2])))
  sid=c(pres$sid,
            paste("NA", 0, 1:nrow(abs), sep="_"))

  
	### ------------------------
	### Tag points (assign sid)
	### ------------------------
	# Tag points
if (inherits(rst_reg_gloproj, "SpatRaster") && level == "glo") {
	  # Extract raster values at point locations
	  vals <- terra::extract(rst_reg_gloproj, out@xy)

	  # Identify if the point is inside or outside
	  gloorreg <- ifelse(is.na(vals), "g", "r")

	  # Build sid
	  out@sid <- paste(sid, gloorreg, sep="_")
	                
	} else {
	  
	  out@sid <- paste(sid, substr(level, 1, 1), sep="_")
	  
	}

  # return
  return(out)
}
