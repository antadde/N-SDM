#' nsdm.savemap
#'
#' Save individual predictive maps
#'
#' @param maps A nsdm.map output (list of rasters filled with predicted values)
#' @param save_path A character string indicating the path where output maps will be saved
#' @param species_name A character string indicating the name of the target taxon
#' @param model_name A character string indicating the name of the target modelling alogirthm
#' @param format A character string indicating the saving format(s) ("rds", "tif" or "both")
#'
#' @return Saved maps in both .rds and/or .tif formats
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.savemap<-function(maps, save_path, species_name, model_name=NULL, format="both"){

if(!is.null(model_name) && model_name=="esm"){
for(m in 1:length(maps)){
map<-raster::readAll(maps[[m]])
model_nameu<-names(maps)[m]

# Save rasters
save_this_path<-paste(save_path, species_name,model_name, model_nameu,sep="/")
suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
f_fit_tif<-gsub("_\\.",".",paste0(save_this_path,"/", names(map),".tif"))
f_fit_rds<-gsub("_\\.",".",paste0(save_this_path,"/", names(map),".rds"))

if(format=="both"){
# TIF
suppressWarnings(raster::writeRaster(map, f_fit_tif, format="GTiff", datatype="INT2U", overwrite=TRUE, options=c("COMPRESS=LZW", "PREDICTOR=2")))
# RDS
saveRDS(map, f_fit_rds, compress = TRUE)
}

if(format=="rds"){
saveRDS(map, f_fit_rds, compress = TRUE)
}

if(format=="tif"){
suppressWarnings(raster::writeRaster(map, f_fit_tif, format="GTiff", datatype="INT2U", overwrite=TRUE, options=c("COMPRESS=LZW", "PREDICTOR=2")))
}

}  
}
  
else {
  
map<-raster::readAll(maps)

# Save rasters
save_this_path<-paste(save_path, species_name, model_name,sep="/")
suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
f_fit_tif<-gsub("_\\.",".",paste0(save_this_path,"/", names(map),".tif"))
f_fit_rds<-gsub("_\\.",".",paste0(save_this_path,"/", names(map),".rds"))

if(format=="both"){
# TIF
suppressWarnings(raster::writeRaster(map, f_fit_tif, format="GTiff", datatype="INT2U", overwrite=TRUE, options=c("COMPRESS=LZW", "PREDICTOR=2")))
# RDS
saveRDS(map, f_fit_rds, compress=TRUE)
}

if(format=="rds"){
saveRDS(map, f_fit_rds, compress=TRUE)
}

if(format=="tif"){
suppressWarnings(raster::writeRaster(map, f_fit_tif, format="GTiff", datatype="INT2U", overwrite=TRUE, options=c("COMPRESS=LZW", "PREDICTOR=2")))
}

}
}