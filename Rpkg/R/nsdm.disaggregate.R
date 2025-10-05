#' nsdm.disaggregate
#'
#' Spatiotemporal Disaggregation of Species Data
#'
#' This function performs spatiotemporal disaggregation of species occurrence data by ensuring a minimum distance between occurrences in space and time, and filtering species with insufficient occurrences for modeling.
#'
#' @param pres A `data.frame` containing species occurrence data with the following key columns:
#'   - `species`: Taxon name.
#'   - `X` and `Y`: Coordinates.
#'   - `year`: Observation year.
#' @param rst A reference `Raster*` object used in subsequent analyses. Defines the spatial resolution and extent for the disaggregation.
#' @param thindist A numeric value specifying the minimum spatial distance (in the units of `rst`) between two occurrences for spatial disaggregation.
#' @param thinyear A numeric value specifying the minimum temporal separation (in years) between two occurrences at the same raster pixel.
#' @param min_occ A numeric value specifying the minimum number of occurrences required for a species to be included in the modeling dataset.
#' @param ncores A numeric value specifying the number of cores to use for parallel operations during disaggregation.
#'
#' @return A `sf` object containing the updated and disaggregated species occurrence data.
#'
#' @details 
#' The function ensures that species occurrences meet the specified spatial and temporal separation criteria. It also filters out species that do not meet the minimum occurrence threshold. The function is optimized for parallel processing when `ncores` is greater than 1.
#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.disaggregate<-function(pres=numeric(), rst, thindist=0, thinyear=0, min_occ=0, ncores=1){

# Convert presences to sf
pres <- st_as_sf(
  pres,
  coords = c("X", "Y"),  # Specify columns for coordinates
  crs = crs(rst)         # Use CRS from the reference raster
)

# List target species
sps <- unique(pres$species)

# Add "year" column if missing (for static data)
if (!"year" %in% names(pres)) {
  pres$year <- NA
}

# Prepare reference raster with target resolution for spatial disaggregation
terraOptions(progress = 0)  # Suppress progress output
fact=round(thindist / res(rst)[1])
rstthin <- terra::aggregate(rst, fact = fact)

# Remove species with occurrences below the minimum threshold for modeling
sp_count <- as.data.frame(table(pres$species))
colnames(sp_count) <- c("species", "freq")
bad_count <- which(sp_count$freq < min_occ)
if (length(bad_count) > 0) {
  pres <- pres[!pres$species %in% sp_count[bad_count, "species"], ]
}

# Update target species list after filtering
sps <- unique(pres$species)

# Remove presences with no environmental data coverage
po <- mclapply(sps, function(i) {
  # Subset presences for the current species
  pres_po <- pres[pres$species == i, ]
  
  # Extract environmental data for the presence points
  extracted <- terra::extract(rst, pres_po, ID = FALSE)
  
  # Filter out rows with missing environmental data
  valid_rows <- !is.na(extracted)
  pres_po <- pres_po[valid_rows, ]
  
  return(pres_po)
}, mc.cores = ncores)

# Remove species with occurrences below the minimum threshold for modeling
species_counts <- unlist(lapply(po, nrow))  # Get the number of occurrences for each species
bad_count <- which(species_counts < min_occ)  # Identify species with insufficient occurrences
if (length(bad_count) > 0) {
  po <- po[-bad_count]  # Remove these species from the list
}

# Update the presence data and the list of species
pres <- do.call(rbind, po)
sps <- unique(pres$species)

# Spatial thinning
if (is.numeric(thindist) && thindist > 0) {
  pq <- mclapply(sps, function(i) {
    
	# Subset presence data for the current species
    pres_pq <- pres[pres$species == i, ]
    
    # Extract coordinates from sf object
    coords <- sf::st_coordinates(pres_pq)  # Extract X, Y coordinates
    
    # Assign cell IDs based on the thinned raster
    pres_pq$cell <- terra::cellFromXY(rstthin, coords)
	   
	# Identify duplicates based on species × cell × year
	dups_idx <- duplicated(sf::st_drop_geometry(pres_pq)[, c("species", "cell", "year")])

	# Filter out the duplicates, keeping the geometry
	pres_pq <- pres_pq[!dups_idx, ]
    
    return(pres_pq)
  }, mc.cores = ncores)
  
  # Remove species with occurrences below the minimum threshold for modeling
  species_counts <- unlist(lapply(pq, nrow))  # Count occurrences per species
  bad_count <- which(species_counts < min_occ)  # Identify species with insufficient occurrences
  if (length(bad_count) > 0) {
    pq <- pq[-bad_count]  # Remove these species from the list
  }
  
# Update the presence data and the list of species
pres <- do.call(rbind, pq)
sps <- unique(pres$species)
}

# Temporal thinning: Ensure minimum period between years within the same cell
if (is.numeric(thinyear) && thinyear > 0) {
  # Parallelize across species
  pp <- mclapply(sps, function(i) {
    # Subset presence data for the current species
    pres_pp <- pres[pres$species == i, ]
    
    # Assign cell IDs based on the thinned raster
    coords <- sf::st_coordinates(pres_pp)  # Extract coordinates
    pres_pp$cell <- terra::cellFromXY(rstthin, coords)
    pres_pp$id <- seq_len(nrow(pres_pp))  # Create unique IDs for each record
    
    # Initialize list to collect valid IDs
    thin_id <- c()
    
    # Loop across unique cells
    for (c in unique(pres_pp$cell)) {
      # Subset records for the current cell
      sub <- pres_pp[pres_pp$cell == c, ]
      
      # Sort records by year
      years <- sub[order(sub$year), ]$year
      
      # Identify and remove years that violate the thinning condition
      while (TRUE) {
        bad_years <- which(diff(years) < thinyear) + 1
        if (length(bad_years) == 0) break
        years <- years[-bad_years[1]]
      }
      
      # Retain only valid years
      sub <- sub[sub$year %in% years, ]
      thin_id <- c(thin_id, sub$id)
    }
    
    # Return thinned presence data for the species
    pres_pp[pres_pp$id %in% thin_id, ]
  }, mc.cores = ncores)
  
  # Remove species with occurrences below the minimum threshold for modeling
  species_counts <- unlist(lapply(pp, nrow))  # Count occurrences per species
  bad_count <- which(species_counts < min_occ)  # Identify species with insufficient occurrences
  if (length(bad_count) > 0) {
    pp <- pp[-bad_count]  # Remove these species from the list
  }
  
  # Retain only valid rows in the presence data
  pres <- do.call(rbind, pp)
  
  # Clean columns
  pres$id<-NULL
  pres$cell<-NULL
}

return(pres)
}
