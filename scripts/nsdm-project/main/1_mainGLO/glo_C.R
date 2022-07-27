#############################################################################
## 1_mainGLO
## C: ensembling and mapping
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/1_mainGLO","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","tmp",gsub("/main/1_mainGLO","",getwd())),"/settings/nsdm-settings.RData"))

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
ispi_name<-species[arrayID]

cat(paste0('Ready for mapping and ensembling GLO predictions obtained for ', ispi_name, ' ...\n'))

### =========================================================================
### C- Save prediction raster
### =========================================================================
for(i in 1:length(mod_algo)){
# Load raw prediction data
model_name<-mod_algo[i]
pred_path<-paste0(scr_path,"/outputs/",project,"/d6_preds/glo")
full_pred_path<-paste0(paste(pred_path, ispi_name, model_name, sep="/"))
pred_file<-list.files(full_pred_path, pattern=".rds", full.names=TRUE)
pred<-readRDS(pred_file)

# Predict 
map_i<-nsdm.map(template=pred$template,
                nona_ix=pred$nona_ix,
                species_name=ispi_name, model_name=model_name, level="glo",
                pred=pred$ndata_bck) 

# Save
nsdm.savemap(maps=map_i, species_name=ispi_name, model_name=model_name, format="rds", save_path=paste0(scr_path,"/outputs/",project,"/d7_maps/glo"))
cat(paste0(model_name,' predictions saved \n'))
}

### =========================================================================
### D- Ensemble predictions
### =========================================================================
ensemble_glo<-nsdm.ensemble(model_names= mod_algo,
                           species_name=ispi_name,
						   level="glo",
                           map_path=paste0(scr_path,"/outputs/",project,"/d7_maps/glo"), 
                           score_path=paste0(scr_path,"/outputs/",project,"/d3_evals/glo"), 
                           weighting=do_weighting, 
                           weight_metric=weight_metric, 
                           discthre=disc_thre) 

nsdm.savemap(maps=ensemble_glo$ensemble, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d8_ensembles/glo"))
nsdm.savemap(maps=ensemble_glo$ensemble_cv, species_name=ispi_name, model_name=NULL, save_path=paste0(scr_path,"/outputs/",project,"/d9_ensembles-cv/glo"))

cat(paste0('Ensemble predictions saved \n'))
cat(paste0('Finished!\n'))