#############################################################################
## 2_mainREG
## C: ensembling, nesting and mapping
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/2_mainREG","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","tmp",gsub("/main/2_mainREG","",getwd())),"/settings/nsdm-settings.RData"))

# Set permissions for new files
Sys.umask(mode="000")

# Set your working directory
setwd(w_path)

# Set lib path
.libPaths(lib_path)

# Load nsdm package
require(nsdm)

### =========================================================================
### B- Definitions
### =========================================================================
# SBATCH param
args<-eval(parse(text=args))
arrayID<-eval(parse(text=arrayID))

# Target species
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

# Scale-nesting methods for combining GLO and REG predictions
nesting_methods<-nesting_methods

# SBATCH array
array<-expand.grid(nesting=nesting_methods, species=species)
ispi_name <- array[arrayID,"species"]
nesting_method <- array[arrayID,"nesting"]

cat(paste0('Ready for mapping and ensembling REG predictions obtained for ', ispi_name, ' with ',nesting_method,' nesting method...\n'))

### =========================================================================
### C- Save prediction raster
### =========================================================================
for(i in 1:length(mod_algo)){
# Load raw prediction data
model_name<-mod_algo[i]
pred_path<-paste0(scr_path,"/outputs/",project,"/d6_preds/reg/", nesting_method)
full_pred_path<-paste0(paste(pred_path, ispi_name, model_name, sep="/"))
pred_file<-list.files(full_pred_path, pattern=".rds", full.names=TRUE)
pred<-readRDS(pred_file)

# Predict 
map_i<-nsdm.map(template=pred$template,
                nona_ix=pred$nona_ix,
				level="reg",
                species_name=ispi_name,
				model_name=model_name,
				nesting_name=nesting_method,
                pred=pred$ndata_bck)

# Save
nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="rds", save_path=paste0(scr_path,"/outputs/",project,"/d7_maps/reg/", nesting_method))
cat(paste0(model_name,' predictions saved \n'))
}
cat(paste0('Predictions saved \n'))

### =========================================================================
### D- Ensemble predictions
### =========================================================================
ensemble_reg<-nsdm.ensemble(model_names=mod_algo, 
                           species_name=ispi_name,
						   nesting_name=nesting_method,
						   level="reg",
                           map_path=paste0(scr_path,"/outputs/",project,"/d7_maps/reg/", nesting_method),
                           score_path=paste0(scr_path,"/outputs/",project,"/d3_evals/reg/", nesting_method), 
                           weighting=do_weighting, 
                           weight_metric=weight_metric, 
                           discthre=disc_thre)

nsdm.savemap(maps=ensemble_reg$ensemble, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d8_ensembles/reg/",nesting_method))

nsdm.savemap(maps=ensemble_reg$ensemble_cv, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d9_ensembles-cv/reg/",nesting_method))

### =========================================================================
### E- Combine REG and GLO predictions
### =========================================================================
# E.1.1 "Multiply" nesting
if(nesting_method=="multiply"){
  ## Ensembling
  # response
  ensemble_glo<-readRDS(list.files(paste0(scr_path,"/outputs/",project,"/d8_ensembles/glo/",ispi_name), pattern=".rds", full.names = TRUE))
  ensemble_nested<-sqrt(ensemble_glo*ensemble_reg$ensemble)
  names(ensemble_nested)<-names(ensemble_reg$ensemble)
  # cv
  ensemble_glo_cv<-readRDS(list.files(paste0(scr_path,"/outputs/",project,"/d9_ensembles-cv/glo/",ispi_name), pattern=".rds", full.names = TRUE))
  ensemble_nested_cv<-raster::mean(raster::stack(ensemble_reg$ensemble_cv, ensemble_glo_cv))
  names(ensemble_nested_cv)<-names(ensemble_reg$ensemble_cv)
  # Save
  nsdm.savemap(map=ensemble_nested, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d10_nested-ensembles/",nesting_method))

  
  nsdm.savemap(map=ensemble_nested_cv, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d11_nested-ensembles-cv/",nesting_method))

} 

# # E.1.2 "Covariate" nesting (identical to D)
if(nesting_method=="covariate"){
## Ensembling
ensemble_nested<-ensemble_reg$ensemble
ensemble_nested_cv<-ensemble_reg$ensemble_cv
# Save
 nsdm.savemap(map=ensemble_nested, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d10_nested-ensembles/",nesting_method))

 nsdm.savemap(map=ensemble_nested_cv, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d11_nested-ensembles-cv/",nesting_method))
}

cat(paste0('GLO and REG predictions nested and saved \n'))
cat(paste0('Finished!\n'))