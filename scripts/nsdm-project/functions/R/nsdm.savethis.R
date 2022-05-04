#' nsdm.savethis
#'
#' Easy saving of .rds nsdm.ine objects
#'
#' @param object Object to be saved
#' @param model_name A character string indicating the name of the target modelling algorithm
#' @param species_name A character string  indicating the name of the target taxon
#' @param tag A character string  indicating the tag for the target hyperparameter X modelling algorithm combination
#' @param compression Logical (TRUE FALSE) for saveRDS "compress" argument
#' @param save_path A character string indicating the upper path where to save the object
#'
#' @return Saved object in .rds format
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.savethis<-function(object, model_name=NULL, species_name=NULL, tag=NULL, compression=FALSE, save_path){
save_this_path<-paste(save_path, species_name, model_name, sep="/")
suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
if(is.null(tag)){f<-gsub("_\\.",".",paste0(save_this_path,"/", paste(species_name,model_name,sep="_"),".rds"))
}else{f<-gsub("_\\.",".",paste0(save_this_path,"/", paste(species_name,tag,sep="_"),".rds"))}
if(class(object)=="lgb.Booster"){
saveRDS.lgb.Booster(object, file=f)
}else{saveRDS(object, file=f, compress=compression)}
}