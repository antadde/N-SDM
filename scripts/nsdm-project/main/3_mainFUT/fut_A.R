#############################################################################
## 3_mainFUT
## A: global-level projections
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
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

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
array<-expand.grid(model=models, species=species, scenarios=rcps)
ispi_name <- array[arrayID,"species"]
model_name <- array[arrayID,"model"]
rcp<-array[arrayID,"scenarios"]

for (per in pers){

cat(paste('Ready for', rcp, per, toupper(model_name),  'future GLO predictions for', ispi_name, '...\n', sep=" "))

### =========================================================================
### C- Load GLO (bioclimatic) model
### =========================================================================
# C.1.1 Load GLO (bioclimatic) data
d0_datasets<-nsdm.loadthis(species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/glo"))

# A.1.2 Load GLO (bioclimatic) model
prmod<-nsdm.loadthis(model_name=model_name, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))

# Specific loading strategy for lgb.booster			  
if("lgb.Booster" %in% class(prmod)){
prmod<-readRDS.lgb.Booster(paste0(scr_path,"/outputs/",project,"/d2_models/glo/",ispi_name,"/gbm/",ispi_name,"_gbm.rds"))
prmod2<-prmod
prmod<-nsdm.loadthis(model_name="glm", species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))
prmod@fits[[1]][[1]]<-prmod2	  
}

# A.2 List covariates	  
cov<-unlist(strsplit(prmod@meta$env_vars,", "))
			  
### =========================================================================
### B. Load future bioclimatic layers
### =========================================================================
# B.1 List available layers
lr_loc<-readRDS(paste0(w_path,"outputs/",project,"/settings/covariates-list.rds"))$lr_loc
lr_fut<-lr_loc[intersect(grep("/future/", lr_loc), grep(paste(per,rcp,sep="/"), lr_loc))]
lr_fut<-grep(paste0(cov,".rds", collapse="|"), lr_fut, value=T)
   
# B.2 Load target layers
clim_stk_fut<-nsdm.fastraster(files=lr_fut, nsplits=ncores)
names(clim_stk_fut)<-gsub(".*_", "", names(clim_stk_fut))
clim_stk_fut<-clim_stk_fut[[which(names(clim_stk_fut) %in% cov)]]
lr_fut<-grep(paste0(cov,".rds", collapse="|"), lr_fut, value=T)

### =========================================================================
### C- Spatial predictions
### =========================================================================
## C.1 Prepare covariate data for predictions
clim_df_loc<-nsdm.retrieve4pred(covstk=clim_stk_fut, # subset for selected covariates
                               scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")]) # scaling parameters to be reapplied

## C.2 Clean workspace to free some memory before predicting
template<-clim_stk_fut[[1]]
rm(d0_datasets, clim_stk_fut)
gc()

## C.3 Predict
ndata_bck<-nsdm.predict(models=prmod,
                        nwdata=clim_df_loc$covdf, # covariate data for predictions
                        nsplits=ncores)

## C.4 Save
nsdm.savethis(object=list(ndata_bck=ndata_bck, template=template, nona_ix=clim_df_loc$covdf_ix),
              model_name=model_name, species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d13_preds-fut/glo/",rcp,"/",per))
}

cat(paste0('Predictions calculated and saved \n'))
cat(paste0('Finished!\n'))