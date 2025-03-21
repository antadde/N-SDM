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
load(paste0(gsub("scripts","tmp",gsub("/main/3_mainFUT","",getwd())),"/settings/nsdm-settings.RData"))

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
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

# Target model algorithms
models<-mod_algo

# Target scenario for future predictions
scenars<-proj_scenarios

# Target period for future predictions
pers<-proj_periods

# SBATCH array
# array<-expand.grid(model=models, species=species, scenarios=scenars)
array<-expand.grid(species=species, scenarios=scenars)
ispi_name <- array[arrayID,"species"]
# model_name <- array[arrayID,"model"]
scenar<-array[arrayID,"scenarios"]


for(model_name in models){

for (per in pers){

cat(paste('Ready for', scenar, per, toupper(model_name),  'future GLO predictions for', ispi_name, '...\n', sep=" "))

### =========================================================================
### C- Load GLO model
### =========================================================================
# C.1.1 Load GLO data
d0_datasets<-nsdm.loadthis(species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/glo"))

# A.1.2 Load GLO model
prmod<-nsdm.loadthis(model_name=model_name, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))

# Specific loading strategy for lgb.booster			  
if("lgb.Booster" %in% class(prmod)){
prmod<-lgb.load(paste0(scr_path,"/outputs/",project,"/d2_models/glo/",ispi_name,"/gbm/",ispi_name,"_gbm.rds"))
prmod2<-prmod
prmod<-nsdm.loadthis(model_name="glm", species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))
prmod@fits[[1]][[1]]<-prmod2	  
}

# A.2 List covariates	  
cov<-unlist(strsplit(prmod@meta$env_vars,", "))
			  
### =========================================================================
### B. Load future layers
### =========================================================================
# Retrieve list of candidate covariates and covinfo table
lr<-readRDS(paste0(w_path,"tmp/",project,"/settings/covariates-list.rds"))
cov_info<-lr$cov_info
cov_info$ID<-paste(cov_info$cada, cov_info$variable, cov_info$attribute, cov_info$focal, sep="_")
cov_ID<-cov_info$ID[match(cov, gsub(".rds", "", basename(cov_info$file)))]
if(n_levels==1) cov_info_pres<-cov_info[cov_info$ID %in% cov_ID & cov_info$period!="future",]
if(n_levels==2) cov_info_pres<-cov_info[cov_info$ID %in% cov_ID & cov_info$period!="future" & cov_info$level=="reg",]
lr_pres_ID<-cov_info_pres$ID

# B.1 List available future layers for cov_ID
cov_info_fut<-cov_info[cov_info$ID %in% cov_ID & cov_info$period=="future" & cov_info$scenario==scenar & cov_info$start_year==gsub("_.*", "", per) & cov_info$end_year==gsub(".*_", "", per),]
lr_fut<-cov_info_fut$file
lr_fut_ID<-cov_info_fut$ID
lr_fut_names<-gsub(".rds", "", basename(cov_info_pres$file))[match(lr_fut_ID, lr_pres_ID)]

# B.2 Load target future layers if available
stk_fut<-stack()
if(length(lr_fut)>0){
stk_fut<-nsdm.fastraster(files=lr_fut, nsplits=ncores)
names(stk_fut)<-lr_fut_names
}

# B.3 Complete with (past) "static" layers if needed
remainders<-setdiff(lr_pres_ID, lr_fut_ID)
if(length(remainders)>0){
lr_remain<-cov_info_pres[cov_info_pres$ID %in% remainders,]$file
lr_remain_names<-gsub(".rds", "", basename(lr_remain))
stk_remain<-nsdm.fastraster(files=lr_remain, nsplits=ncores)
names(stk_remain)<-lr_remain_names
stk_fut<-stack(stk_fut, stk_remain)
}

### =========================================================================
### C- Spatial predictions
### =========================================================================
## C.1 Prepare covariate data for predictions
clim_df_reg<-nsdm.retrieve4pred(covstk=stk_fut,
                               scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")])

## C.2 Clean workspace to free some memory before predicting
template<-stk_fut[[1]]
rm(d0_datasets, stk_fut)
gc()

## C.3 Predict
ndata_bck<-nsdm.predict(models=prmod,
                        nwdata=clim_df_reg$covdf,
                        nsplits=ncores)

## C.4 Save
nsdm.savethis(object=list(ndata_bck=ndata_bck, template=template, nona_ix=clim_df_reg$covdf_ix),
              model_name=model_name, species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d13_preds-fut/glo/",scenar,"/",per))
}
}

cat(paste0('Predictions calculated and saved \n'))
cat(paste0('Finished!\n'))