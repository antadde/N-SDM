#############################################################################
## 3_mainFUT: Pipeline for predicting future species distributions
## C: Compute local (habitat) projections for future horizons
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
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# SBATCH param
args<-eval(parse(text=args))
arrayID<-eval(parse(text=arrayID))

# Target species
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))

# Target model algorithms
models<-mod_algo

# Scale-nesting methods for combining GLO and LOC predictions
nesting_methods<-nesting_methods

# Target RCP for future predictions
rcps<-proj_scenarios

# Target period for future predictions
pers<-proj_periods

# SBATCH array
array<-expand.grid(nesting=nesting_methods, model=models, species=species, scenarios=rcps)
ispi_name <- array[arrayID,"species"]
model_name <- array[arrayID,"model"]
nesting_method <- array[arrayID,"nesting"]
rcp<-array[arrayID,"scenarios"]


for (per in pers){

cat(paste('Ready for', rcp, per, toupper(model_name),  'future LOC predictions of', ispi_name, ' using the ',nesting_method,' method for scale-nesting ...\n'))

### =========================================================================
### C- Load LOC model
### =========================================================================
# C.1.1 Load LOC data
d0_datasets<-nsdm.loadthis(species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/loc"))

# C.1.2 Load LOC model
prmod<-nsdm.loadthis(model_name=model_name, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/loc/",nesting_method))

# Specific loading strategy for lgb.booster			  
if("lgb.Booster" %in% class(prmod)){
prmod<-readRDS.lgb.Booster(paste0(scr_path,"/outputs/",project,"/d2_models/loc/",nesting_method,"/",ispi_name,"/gbm/",ispi_name,"_gbm.rds"))
prmod2<-prmod
prmod<-nsdm.loadthis(model_name="glm", species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d2_models/loc/",nesting_method))
prmod@fits[[1]][[1]]<-prmod2	  
}

# C.2 List covariates	  
cov<-unlist(strsplit(prmod@meta$env_vars,", "))
			  
### =========================================================================
### D- Load future layers
### =========================================================================
# D.1 List available layers
lr_loc<-readRDS(paste0(w_path,"outputs/",project,"/settings/covariates-list.rds"))$lr_loc
lr_fut<-lr_loc[intersect(grep("/future/", lr_loc), grep(paste(per,rcp,sep="/"), lr_loc))]
lr_fut<-grep(paste0(cov,".rds", collapse="|"), lr_fut, value=T)
		   
# D.2 Load target future layers if available
hab_stk_fut<-stack()
if(length(lr_fut)>0){
hab_stk_fut<-nsdm.fastraster(files=lr_fut, nsplits=ncores)
names(hab_stk_fut)<-gsub(".*_", "", names(hab_stk_fut)) # Rename to match with glo names
hab_stk_fut<-hab_stk_fut[[which(names(hab_stk_fut) %in% cov_glo)]]}

# D.3 Load future GLO output if needed
if("mainGLO" %in% cov){
glo_out<-list.files(paste0(scr_path,"/outputs/",project,"/d15_ensembles-fut/glo/",rcp,"/",per,"/",ispi_name), pattern=".rds", full.names = TRUE)
hab_stk_fut<-stack(hab_stk_fut, readRDS(glo_out))
names(hab_stk_fut)[nlayers(hab_stk_fut)]<-"mainGLO"
}

# D.4 Complete with (past) "static" layers if needed
remainders<-setdiff(cov, names(hab_stk_fut))
if(length(remainders)>0){
lr_remain<-grep(paste0(remainders,".rds", collapse="|"), lr_loc, value=T)
hab_stk_fut<-stack(hab_stk_fut, nsdm.fastraster(files=lr_remain, nsplits=ncores))
}
 
### =========================================================================
### E- Spatial predictions
### =========================================================================
## E.1 Prepare covariate data for predictions
hab_df_loc<-nsdm.retrieve4pred(covstk=hab_stk_fut,
                               observational=grep(paste0(cov_obser, collapse="|"), names(hab_stk_fut), value=T),# Flatten observational covariates
							   obsval=cov_observ_val,
							   mask=mask_pred, # mask to be applied on predictions
                               scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")]) # scaling parameters to be reapplied

## E.2 Clean workspace to free some memory before predicting
template<-hab_stk_fut[[1]]
rm(d0_datasets, hab_stk_fut)
gc()

## E.3 Predict
ndata_bck<-nsdm.predict(models=prmod,
                        nwdata=hab_df_loc$covdf, # covariate data for predictions
                        nsplits=ncores)

## E.4 Save
nsdm.savethis(object=list(ndata_bck=ndata_bck, template=template, nona_ix=hab_df_loc$covdf_ix),
              model_name=model_name, species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d13_preds-fut/loc/",nesting_method,"/",rcp,"/",per))
}

cat(paste0('Predictions calculated and saved \n'))
cat(paste0('Finished!\n'))