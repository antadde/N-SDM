#############################################################################
## 1_mainGLO: Pipeline for fitting global-level (bioclimatic) species distribution models
## B: Modelling and predictions
## Date: 25-09-2021
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/1_mainGLO","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","outputs",gsub("/main/1_mainGLO","",getwd())),"/settings/nsdm-settings.RData"))

# Set permissions for new files
Sys.umask(mode="000")

# Set your working directory
setwd(w_path)

# Set lib path
.libPaths(lib_path)

# Load required packages
invisible(lapply(c(packs_data, packs_modl, packs_eval), require, character.only = TRUE))

# Source custom functions
invisible(sapply(list.files(f_path, pattern = ".R", full.names = TRUE, recursive=TRUE), source))

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

# SBATCH array
array<-expand.grid(model=models, species=species)
ispi_name <- array[arrayID,"species"]
model_name <- array[arrayID,"model"]

cat(paste0('Ready for global-level ', toupper(model_name), ' modelling of ',ispi_name, '...\n'))

### =========================================================================
### C- Load glo_A outputs
### =========================================================================
## Load modelling sets
d0_datasets<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/glo/",ispi_name,"/",ispi_name,".rds"))

## Load selected covariates
d1_covsels<-readRDS(paste0(scr_path,"/outputs/",project,"/d1_covsels/glo/",ispi_name,"/",ispi_name,".rds"))

cat(paste0('glo_A outputs loaded \n'))

### =========================================================================
### D- Define model parameters 
### =========================================================================
# Define all possible combinations of model parameters
modinp<-nsdm.setparam(model_name=model_name, covariate_names=names(d1_covsels$pseu.abs_i@env_vars),
                      param_grid=param_grid,
                      tmp_path=paste0(scr_path,"/tmp/",project),					  
                      weights=d0_datasets$weights,
                      ncov.esm=ncov_esm, comb.esm=comb_esm, # ESM settings
                      nthreads=ncores)

cat(paste0('Modelling parameters defined \n'))

### =========================================================================
### E- Fit, evaluate and save models
### =========================================================================
eval_list<-list()

suppressWarnings(for(i in 1:length(modinp)){
modinp_i<-modinp[i]

## E.1 Fit model
cat(paste0('\n fitting ',modinp_i[[1]]@tag, '...\n'))
mod<-NULL
ptm <- proc.time()
mod<-try(nsdm.flex2(x=d1_covsels$pseu.abs_i,
                   replicatetype=replicate_type,
                   reps=reps, 
                   mod_args=modinp_i,
				   ncores=ncores,
				   level="glo",
				   tmp_path=paste0(scr_path,"/tmp/",project)),TRUE)
timer<-c(proc.time() - ptm)

if(all(lapply(mod@fits, class) == "try-error")) next

## E.2 Evaluate model
evals<-NULL
try(evals<-nsdm.evaluate2(mod,crit=eval_crit, ncores=ncores,level="glo",
				          tmp_path=paste0(scr_path,"/tmp/",project)),TRUE)
  if(class(evals)!="NULL"){
    smev<-nsdm.summary(evals)
    t<-rbind(t.user=timer[1],t.elapsed=timer[3])
    smev<-rbind(smev,t)
    print(smev)
    eval_list[[i]]<-smev}

# E.3 Save model
  if(class(mod)!="NULL"){
    nsdm.savethis(object=list(model=mod, parameters=modinp_i),
                  species_name=ispi_name, model_name=model_name,
                  tag=modinp_i[[1]]@tag,
                  save_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))}	
})

## E.4 Save evaluation table
eval_list<-eval_list[lengths(eval_list) != 0] # filter for possible model failures 
nsdm.savethis(object=eval_list,
              model_name=model_name, species_name=ispi_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d3_evals/glo"))

cat(paste0('\n\nModels fitted and evaluated \n'))

### =========================================================================
### F- Refit top model for prediction
### =========================================================================
## F.1.1 Identify best model
if(model_name!="esm"){
  if(length(eval_list)>1){
    smev<-do.call(cbind,eval_list)
    ord<-sort(smev[best_met,],decreasing=T) # Target metric for identifying "best" model
    modinp_top<-modinp[names(ord[1])]
  } else {
    modinp_top<-modinp}}
## F.1.2 ... or discard models with Score < (for esm)
if(model_name=="esm"){
  smev<-do.call(cbind,eval_list)
  ord<-sort(smev[best_met,],decreasing=T) # Target metric for identifying "best" model
  ord<-ord[ord>best_thre_esm]
  modinp_top<-modinp[names(ord)]
}

## F.2a Refit model(s) using full dataset
if(model_name!="esm"){
suppressWarnings(prmod<-nsdm.flex2(x=d1_covsels$pseu.abs_i,
                                  replicatetype="none",
                                  reps=1,
								  ncores=1,
								  level="glo",
                                  mod_args=modinp_top,
								  tmp_path=paste0(scr_path,"/tmp/",project)))
}

## F.2b Refit model(s) using full dataset for ESMs
if(model_name=="esm"){
suppressWarnings(prmod<-nsdm.flex(x=d1_covsels$pseu.abs_i,
                                  replicatetype="none",
                                  reps=1,
                                  mod_args=modinp_top))
}
								  
nsdm.savethis(object=prmod,
              species_name=ispi_name, model_name=model_name, tag=model_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))
			  
if(model_name=="gbm") saveRDS.lgb.Booster(prmod@fits[[1]][[1]], paste0(scr_path,"/outputs/",project,"/d2_models/glo/",ispi_name,"/gbm/",ispi_name,"_",model_name,".rds"))
	
								  
cat(paste0('Top model ',names(modinp_top), ' refitted on full dataset for predictions \n'))

### =========================================================================
### G- Compute covariate importance and response curves
### =========================================================================
## G.1 Covariate importance
print(imp<-nsdm.varimp(prmod))
nsdm.savethis(object=imp,model_name=model_name, species_name=ispi_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d4_varimps/glo"))

## G.2 Response curves
Data<-d1_covsels$pseu.abs_i@env_vars
respcurves<-nsdm.respcurve(prmod, Data=Data,factor=100,
                           scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")],
                           model_name=model_name, species_name=ispi_name,
			               plotting=TRUE, ncores=ncores, save_path=paste0(scr_path,"/outputs/",project,"/plots/respcurves/glo")) 
						   
nsdm.savethis(object=respcurves, model_name=model_name, species_name=ispi_name,
              compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d5_respcurves/glo"))

cat(paste0('\n\nVariable importance scores and response curves computed \n'))

### =========================================================================
### H- Spatial predictions
### =========================================================================
## H.1 Prepare covariate data for predictions
clim_df_loc<-nsdm.retrieve4pred(covstk=d1_covsels$covstk, # subset for selected covariates
                               scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")]) # scaling parameters to be reapplied

## H.2 Clean workspace to free some memory before predicting
template<-d1_covsels$covstk[[1]]
rm(d0_datasets, d1_covsels, respcurves, imp, eval_list)
gc()

## H.3 Predict
ndata_bck<-nsdm.predict(models=prmod,
                        nwdata=clim_df_loc$covdf, # covariate data for predictions
                        nsplits=ncores)
						
nsdm.savethis(object=list(ndata_bck=ndata_bck, template=template, nona_ix=clim_df_loc$covdf_ix),
              model_name=model_name, species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d6_preds/glo"))

cat(paste0('Predictions calculated and saved \n'))

cat(paste0('Finished!\n'))