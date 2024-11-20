#' nsdm.map
#'
#' Fill empty template rasters with predicted values
#'
#' @param template An empty raster layer to fill with predicted values
#' @param nona_ix Numeric indices of cells that should be filled with predicted values
#' @param pred List containing vectors of fitted values (one list by modelling algorithm)
#' @param species_name A character string indicating the name of the taxon for which models are fitted
#' @param model_name A character vector indicating the name of the modelling algorithm(s) for which predicted values are provided
#' @param level A character vector indicating the name of the level considered (e.g., glo or loc)
#' @param nesting_name A character vector indicating the name of the nesting strategy for which predicted values are provided
#' @param scenar_name A character vector indicating the name of the scenario for which predicted values are provided
#' @param period_name A character vector indicating the name of the period for which predicted values are provided
#'
#' @return A list of rasters filled with predicted values
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.map<-function(template, nona_ix, pred, species_name, model_name, level=NA, nesting_name=NA, scenar_name=NA, period_name=NA){

pred_fit_f<-list() # list where results will be stored
  
for(m in 1:length(pred)){ # Loop over models (1 in regular cases but much more in esm settings)
    
# Retrieve model fit
if(model_name=="esm")fit<-pred[[m]]$fit
if(model_name!="esm")fit<-pred$fit

# In case of esm model uses its name its index for naming
if(model_name=="esm"){
model_nameu<-names(pred)[m]
}else{
model_nameu<-model_name}    

# Fill template rasters
pred_fit<-template
pred_fit[]<-NA
pred_fit[nona_ix]<-round(fit*100)
names(pred_fit)<-gsub("_NA", "", paste(species_name, model_nameu, level, nesting_name, scenar_name, period_name, sep="_"))
storage.mode(pred_fit[]) = "integer"

pred_fit_f[[m]]<-pred_fit
names(pred_fit_f)[m]<-model_nameu

}
if(model_name=="esm")return(pred_fit_f)
if(model_name!="esm")return(pred_fit_f[[1]])}