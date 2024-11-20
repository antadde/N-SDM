#' nsdm.ensemble
#'
#' Ensemble prediction surfaces obtained by individual algorithms
#'
#' @param model_names A character string indicating the modelling algorithms to ensemble
#' @param species_name A character string indicating the name of the taxon
#' @param map_path A character string indicating the path where algorithm-specific prediction surfaces are stored
#' @param score_path A character string indicating the path where algorithm-specific evaluation metrics are stored
#' @param discthre A numeric value for the threshold for weight_metric under which to discard an algorithm
#' @param weighting Logical. If TRUE, weight individual algorithms by weight_metric
#' @param weight_metric A character string indicating the name of the metric used for weighting or discarding
#' @param level A character vector indicating the name of the level considered (e.g., glo or reg)
#' @param nesting_name A character vector indicating the name of the nesting strategy for which predicted values are provided
#' @param scenar_name A character vector indicating the name of the scenario for which predicted values are provided
#' @param period_name A character vector indicating the name of the period for which predicted values are provided
#'
#' @return A list of two rasters: "ensemble" (average) and "ensemble_cv" (coefficient of variation)
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.ensemble <- function(model_names, species_name, level=NA, nesting_name=NA, scenar_name=NA, period_name=NA, map_path, score_path=NULL, discthre=NULL, weighting=FALSE, weight_metric="Score"){
  
# Retrieve prediction maps
stack_map<-raster::stack()
for(i in 1:length(model_names)){
  model_name<-model_names[i]
  full_map_path<-paste0(paste(map_path, species_name, model_name,sep="/"))
  map_files<-list.files(full_map_path, pattern=".rds", full.names = TRUE, recursive = TRUE)
  map2<-lapply(map_files, readRDS)
  map<-raster::stack(unlist(map2))
  stack_map<-raster::stack(stack_map, map)
}

# Compute ensembling evaluation table
res<-data.frame(matrix(nrow=raster::nlayers(stack_map), ncol = 3))
colnames(res)<-c("model_name", "score", "discard")
for(i in 1:length(model_names)){
  model_name<-model_names[i]
  full_score_path<-paste0(paste(score_path, species_name, model_name,sep="/"), paste0("/",species_name,"_",model_name,".rds"))
  score<-readRDS(full_score_path)
## Identify selected model(s) and retrieve scores
  if(model_name=="esm"){
    esm_ix<-grep("_esm",names(stack_map))
    esm_names<-paste("esm", stri_extract_first_regex(names(stack_map)[esm_ix], "[0-9]+"), sep="-")
    score_val<-score[weight_metric,esm_names]
    ## store outputs
    res[esm_ix,"score"]<-score_val
    res[esm_ix,"model_name"]<-esm_names
  } else {
  if(length(score)>1){
    score_val<-data.frame(sort(score[weight_metric,],decreasing=T)[1])
    ## store outputs
    res[i,"score"]<-as.numeric(score_val)
    res[i,"model_name"]<-model_name
  } else {
   score_val<-score[weight_metric,]
   ## store outputs
   res[i,"score"]<-score_val
   res[i,"model_name"]<-model_name
  }
  }
  }

# Check if models fulfill discard threshold, and remove if needed
## Check
if(!"esm" %in% model_names){
if(discthre!="NULL"){
res[,"discard"]<-res[,"score"]<discthre
} else {
res[,"discard"]<-"FALSE" 
}

## Discard
if(discthre!="NULL"){
stack_map<-raster::dropLayer(stack_map, c(which(as.logical(res[,"discard"]))))
res<-res[-which(as.logical(res[,"discard"])),]
}
}

# do weighted ensemble
if(weighting){
  stack_map<-raster::dropLayer(stack_map, c(which(as.logical(res[,"discard"]))))
  ensemble<-raster::weighted.mean(stack_map, w=as.numeric(res[,"score"]), na.rm=T)
} else {
# or un-weighted
ensemble <- raster::mean(stack_map)
}
ensemble_mn<-ensemble

# Compute coefficient of variation	  
rasterstack_sd_fast <- function(x) {
  s0 <- nlayers(x)
  s1 <- raster(x, layer=1)
  s2 <- s1^2
  for(ri in 2:s0) {
    r <- raster(x, layer=ri)
    s1 <- s1 + r
    s2 <- s2 + r^2
  }
  sqrt((s0 * s2 - s1 * s1)/(s0 * (s0 - 1)))
}

ensemble_sd<-rasterstack_sd_fast(stack_map) 
ensemble_cv<-(ensemble_sd/ensemble_mn)*100

# return
ensemble_cv<-round(ensemble_cv)
ensemble_mn<-round(ensemble_mn)

storage.mode(ensemble_cv[]) = "integer"
storage.mode(ensemble_mn[]) = "integer"

crs(ensemble_cv)<-proj4string(crs(readRDS(map_files[1])))
crs(ensemble_mn)<-proj4string(crs(readRDS(map_files[1])))

names(ensemble_mn)<-gsub("_NA", "", paste(ispi_name, level, nesting_name, scenar_name, period_name, "ensemble", sep="_"))
names(ensemble_cv)<-gsub("_NA", "", paste(ispi_name, level, nesting_name, scenar_name, period_name, "ensemble_cv", sep="_"))

return(list(ensemble=ensemble_mn, ensemble_cv=ensemble_cv))
}
