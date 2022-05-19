#############################################################################
## 3_mainFUT: Pipeline for predicting future species distributions
## B: Mapping and ensembling of predictions for future horizons
## Date: 25-09-2021
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

# Target model algorithms
models<-mod_algo

# Target RCP for future predictions
rcps<-proj_scenarios

# Target period for future predictions
pers<-proj_periods

# SBATCH array
array<-expand.grid(species=species, scenarios=rcps)
ispi_name <- array[arrayID,"species"]
rcp<-array[arrayID,"scenarios"]

for (per in pers){

cat(paste('Ready for mapping and ensembling', rcp, per, 'future GLO predictions obtained for', ispi_name, '...\n', sep=" "))

### =========================================================================
### C- Save prediction raster
### =========================================================================
for(i in 1:length(mod_algo)){
# C.1 Load raw prediction data
model_name<-mod_algo[i]
pred_path<-paste0(scr_path,"/outputs/",project,"/d13_preds-fut/glo/",rcp,"/",per)
full_pred_path<-paste0(paste(pred_path, ispi_name, model_name, sep="/"))
pred_file<-list.files(full_pred_path, pattern=".rds", full.names=TRUE)
pred<-readRDS(pred_file)

# C.2 Predict 
map_i<-nsdm.map(template=pred$template,
                nona_ix=pred$nona_ix, # indices for non-na cells (stored in clim_df_loc)
                species_name=ispi_name,
				model_name=model_name,
				level="glo",
				rcp_name=rcp,
				period_name=per,
                pred=pred$ndata_bck) 

# C.3 Save
nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="rds", save_path=paste0(scr_path,"/outputs/",project,"/d14_maps-fut/glo/",rcp,"/",per))
cat(paste0(model_name,' predictions saved \n'))
}

### =========================================================================
### D- Ensemble predictions
### =========================================================================
ensemble_glo<-nsdm.ensemble(model_names= mod_algo, # models for ensembling
                           species_name=ispi_name,
						   level="glo",
						   rcp_name=rcp,
						   period_name=per,
                           map_path=paste0(scr_path,"/outputs/",project,"/d14_maps-fut/glo/",rcp,"/",per), # path where prediction rasters are stored
                           score_path=paste0(scr_path,"/outputs/",project,"/d3_evals/glo"), # path where model evaluation tables are stored
                           weighting=do_weighting, # use weights when ensembling
                           weight_metric=weight_metric, # evaluation metric for weighting/discarding
                           discthre=disc_thre) # threshold under which to discard a model

nsdm.savemap(maps=ensemble_glo$ensemble, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d15_ensembles-fut/glo/",rcp,"/",per))
nsdm.savemap(maps=ensemble_glo$ensemble_cv, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d16_ensembles-cv-fut/glo/",rcp,"/",per))

}

cat(paste0('Ensemble predictions saved \n'))
cat(paste0('Finished!\n'))