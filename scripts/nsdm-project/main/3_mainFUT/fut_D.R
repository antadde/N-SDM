#############################################################################
## 3_mainFUT
## D: ensembling, nesting and mapping local-level projections
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/3_mainFUT","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","outputs",gsub("/main/3_mainFUT","",getwd())),"/settings/nsdm-settings.RData"))

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
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))

# Target RCP for future predictions
scenars<-proj_scenarios

# Scale-nesting methods for combining GLO and LOC predictions
nesting_methods<-nesting_methods

# SBATCH array
array<-expand.grid(nesting=nesting_methods, species=species, scenarios=scenars)
ispi_name <- array[arrayID,"species"]
nesting_method <- array[arrayID,"nesting"]
scenar<-array[arrayID,"scenarios"]

# Target period for future predictions
pers<-proj_periods

for (per in pers){

cat(paste('Ready for mapping and ensembling', scenar, per, 'future LOC predictions obtained for', ispi_name, 'under', nesting_method, 'nesting method...\n', sep=" "))

### =========================================================================
### C- Save prediction raster
### =========================================================================
for(i in 1:length(mod_algo)){
# Load raw prediction data
model_name<-mod_algo[i]
pred_path<-paste0(scr_path,"/outputs/",project,"/d13_preds-fut/loc/",nesting_method,"/",scenar,"/",per)
full_pred_path<-paste0(paste(pred_path, ispi_name, model_name, sep="/"))
pred_file<-list.files(full_pred_path, pattern=".rds", full.names=TRUE)
pred<-readRDS(pred_file)

# Predict 
map_i<-nsdm.map(template=pred$template,
                nona_ix=pred$nona_ix, # indices for non-na cells (stored in clim_df_loc)
                species_name=ispi_name,
				model_name=model_name,
				level="loc",
				scenar_name=scenar,
				period_name=per,
				nesting_name=nesting_method,
                pred=pred$ndata_bck) 

# Save
nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="rds", save_path=paste0(scr_path,"/outputs/",project,"/d14_maps-fut/loc/",nesting_method,"/",scenar,"/",per))
cat(paste0(model_name,' predictions saved \n'))
}

### =========================================================================
### D- Ensemble predictions
### =========================================================================
ensemble_loc<-nsdm.ensemble(model_names= mod_algo, # models for ensembling
                           species_name=ispi_name,
						   level="loc",
						   scenar_name=scenar,
				           period_name=per,
				           nesting_name=nesting_method,
                           map_path=paste0(scr_path,"/outputs/",project,"/d14_maps-fut/loc/",nesting_method,"/",scenar,"/",per), # path where prediction rasters are stored
                           score_path=paste0(scr_path,"/outputs/",project,"/d3_evals/loc/", nesting_method), # path where model evaluation tables are stored
                           weighting=do_weighting, # use weights when ensembling
                           weight_metric=weight_metric, # evaluation metric for weighting/discarding
                           discthre=disc_thre) # threshold under which to discard a model

nsdm.savemap(maps=ensemble_loc$ensemble, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d15_ensembles-fut/loc/",nesting_method,"/",scenar,"/",per))
nsdm.savemap(maps=ensemble_loc$ensemble_cv, species_name=ispi_name, model_name=NULL, format="rds", save_path=paste0(scr_path,"/outputs/",project,"/d16_ensembles-cv-fut/loc/",nesting_method,"/",scenar,"/",per))

### =========================================================================
### E- Combine LOC and GLO predictions
### =========================================================================
# E.1 "Multiply" nesting
if(nesting_method=="multiply"){
  # response
  ensemble_glo<-readRDS(list.files(paste0(scr_path,"/outputs/",project,"/d15_ensembles-fut/glo/",scenar,"/",per,"/",ispi_name), pattern=".rds", full.names = TRUE))
  ensemble_nested<-sqrt(ensemble_glo*ensemble_loc$ensemble)
  names(ensemble_nested)<-names(ensemble_loc$ensemble)
  # # cv
  ensemble_glo_cv<-readRDS(list.files(paste0(scr_path,"/outputs/",project,"/d16_ensembles-cv-fut/glo/",scenar,"/",per,"/",ispi_name), pattern=".rds", full.names = TRUE))
  ensemble_nested_cv<-raster::mean(raster::stack(ensemble_loc$ensemble_cv, ensemble_glo_cv))
  names(ensemble_nested_cv)<-names(ensemble_loc$ensemble_cv)
  # Save
nsdm.savemap(map=ensemble_nested, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d17_nested-ensembles-fut/",nesting_method,"/",scenar,"/",per))
nsdm.savemap(map=ensemble_nested_cv, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d18_nested-ensembles-cv-fut/",nesting_method,"/",scenar,"/",per))

} 

# E.2 "Covariate" nesting
if(nesting_method=="covariate"){
ensemble_nested<-ensemble_loc$ensemble
ensemble_nested_cv<-ensemble_loc$ensemble_cv
# Save
nsdm.savemap(map=ensemble_nested, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d17_nested-ensembles-fut/",nesting_method,"/",scenar,"/",per))
nsdm.savemap(map=ensemble_nested_cv, species_name=ispi_name, save_path=paste0(scr_path,"/outputs/",project,"/d18_nested-ensembles-cv-fut/",nesting_method,"/",scenar,"/",per))
}

}

cat(paste0('GLO and LOC predictions nested and saved \n'))
cat(paste0('Finished!\n'))