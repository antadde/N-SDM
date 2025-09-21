#' nsdm.loadthis
#'
#' Easy loading of nsdm.ine objects
#'
#' @param model_name A character string indicating the name of the modelling alogirthm related to the target object to be loaded
#' @param species_name  A character string  indicating the name of the species related to the target object to be loaded
#' @param tag  A character string indicating the tag for the hyperparameter X modelling algorithm combination of the target object to be loaded
#' @param read_path A character string indicating the upper path where the target object is stored (above model and species)
#'
#' @return The desired nsdm.ine object
#' @author Antoine Adde (aadde@unil.ch)
#' @export
nsdm.loadthis<-function(model_name=NULL, species_name=NULL, tag=NULL, read_path){

read_this_path<-paste(read_path, species_name, model_name, sep="/")

if(is.null(tag)){f<-gsub("_\\.",".",paste0(read_this_path,"/", paste(species_name,model_name,sep="_"),".rds"))
}else{f<-gsub("_\\.",".",paste0(read_this_path,"/", paste(species_name,tag,sep="_"),".rds"))}

#if(model_name=="gbm"){
#readRDS.lgb.Booster(f)
#}else{
readRDS(f)
#}

}