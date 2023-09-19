#' nsdm.disaggregate
#'
#' Spatiotemporal disaggregation of species data
#'
#' @param pres Data.frame containing species occurence data with key columns "species" (taxon name), "year" (observation year), "uncertainty" (coordinate uncertainty) 
#' @param rst Reference raster used in subsequent analyses
#' @param thindist Distance (numeric; in rst units) used for spatial disaggregation (minimal distance between two occurences) 
#' @param thinyear Time period (numeric; in years) used for temporal disaggregation (minimal number of years between two occurences at same pixel)
#' @param max_uncertain Maximum coordinate uncertainty (numeric; in "uncertainty unit) for an occurence to be retained in the modelling set
#' @param min_occ Minimal number of occurences (numeric) required for a species to be modelled
#' @param ncores Number of cores (numeric) used during parallel operations
#'
#' @return An updated Spatial point dataframe object with disaggregated species occurence data
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.disaggregate<-function(pres=numeric(), rst, thindist=0, thinyear=0, max_uncertain=9999, min_occ=0, ncores=1){

# SpatialPointDataFrame
pres <- SpatialPointsDataFrame(coords      = pres[,c("X","Y")],
                               data        = pres, 
                               proj4string = crs(rst))
# List target species
sps<-unique(pres$species)

# If "year" column is missing (static data) add one
if(!"year" %in% names(pres)) pres$year<-NA

# Prepare reference raster with target resolution for spatial dissagregation  
rstthin<-aggregate(rst, fact=(thindist/res(rst))[1])

# Remove presences with uncertainty (coordinates uncertainty) > max_uncertain
bad_uncertain<-which(pres$uncertainty>max_uncertain)
if(length(bad_uncertain)>0) pres<-pres[-bad_uncertain,]
# Remove species with n_occ under min_occ (minimal number of occurences for modelling)
sp_count<-count(pres@data, "species")
bad_count<-which(sp_count$freq<min_occ)
if(length(bad_count)>0) pres<-pres[!pres$species %in% sp_count[bad_count,"species"],]
# update sps
sps<-unique(pres$species)

# Remove presences with no environmental data coverage
po<-mclapply(sps, function(i){
pres_po<-pres[pres$species==i,]
extracted<-extract(rst, pres_po)
sna<-!is.na(extracted)
pres_po<-pres_po[sna,]}, mc.cores=ncores)
# Remove species with n_occ under min_occ (minimal number of occurences for modelling)
bad_count<-which(unlist(lapply(po, length))<min_occ)
if(length(bad_count)>0) po<-po[-bad_count]
pd<-unlist(lapply(po, function(p){row.names(p@data)}))
pres<-pres[which(row.names(pres) %in% pd),]
# update sps
sps<-unique(pres$species)

# Spatial thinning (keep one value per species x cell x year combination)
if(is.numeric(thindist) && thindist>0){
pq<-mclapply(sps, function(i){
pres_pq<-pres[pres$species==i,]
pres_pq$cell<-cellFromXY(rstthin, pres_pq)
pres_dups <- data.frame(pres_pq@data[c("species", "cell", "year")])
if(TRUE %in% duplicated(pres_dups)) pres_pq<-pres_pq[!duplicated(pres_dups),]
return(pres_pq)
}, mc.cores=ncores)
# Remove species with n_occ under min_occ (minimal number of occurences for modelling)
bad_count<-which(unlist(lapply(pq, length))<min_occ)
if(length(bad_count)>0) pq<-pq[-bad_count]
pq<-unlist(lapply(pq, function(p){row.names(p@data)}))
pres<-pres[which(row.names(pres) %in% pq),]
}
# update sps
sps<-unique(pres$species)

# Temporal thinning (for cells with multiple years ensure minimum period between years)
if(is.numeric(thinyear) && thinyear>0){
# Parallelize across species
pp<-mclapply(sps, function(i){
pres_pp<-pres[pres$species==i,]
pres_pp$cell<-cellFromXY(rstthin, pres_pp)   
pres_pp$id<-1:length(pres_pp)
thin_id<-c()
# Loop across cells
for(c in pres_pp$cell){
sub<-pres_pp[pres_pp$cell==c,]
years<-sub[order(sub$year),]$year
bad_years<-which(diff(years)<thinyear)+1
# Identify years to be removed for meeting thinyear condition
while (length(bad_years) > 0)
  {
years<-years[-bad_years[1]]
bad_years<-which(diff(years)<thinyear)+1
  }
# Remove years
sub<-sub[sub$year %in% years,]
thin_id<-c(thin_id, sub$id)
}
pres_pp[pres_pp$id %in% thin_id,]
}, mc.cores=ncores)
# Remove species with n_occ under min_occ (minimal number of occurences for modelling)
bad_count<-which(unlist(lapply(pp, length))<min_occ)
if(length(bad_count)>0) pp<-pp[-bad_count]
pp<-unlist(lapply(pp, function(p){row.names(p@data)}))
pres<-pres[which(row.names(pres) %in% pp),]
}

return(pres)
}
