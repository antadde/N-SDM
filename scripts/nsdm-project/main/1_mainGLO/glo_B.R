#############################################################################
## 1_mainGLO
## B: modelling and predictions
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
# Number of cores to be used during parallel operations
ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# SBATCH param
args<-eval(parse(text=args))
arrayID<-eval(parse(text=arrayID))

# Target species
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

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
## E.1 Fit model
mod<-try(nsdm.flex3(x=d1_covsels$pseu.abs_i,
                   replicatetype=replicate_type,
                   reps=reps, 
                   mod_args=modinp,
				   ncores=ncores,
				   level="glo",
				   timer=TRUE,
				   tmp_path=paste0(scr_path,"/tmp/",project)),TRUE)

## E.2 Evaluate model
### Assessment metrics
try(evals<-nsdm.evaluate3(mod$mod,crit=eval_crit, ncores=ncores, level="glo",
				          tmp_path=paste0(scr_path,"/tmp/",project)),TRUE)

### Computation time		  				  
t<-list()
for(k in 1:length(mod$time)){t_k<-mod$time[[k]]; t[[k]]<-t_k}
t<-do.call(cbind,lapply(1:length(mod$time), function(z){rowMeans(as.matrix(t[[z]]))}))
colnames(t)<-names(mod$mod@fits)

### Bind
  if(class(evals)!="NULL"){
    smev<-nsdm.summary(evals)
	smev<-rbind(smev, t)
    print(smev)}

# E.3 Save model
  if(class(mod)!="NULL"){
  suppressWarnings(nsdm.savethis(object=list(model=mod, parameters=modinp),
                  species_name=ispi_name, model_name=model_name,
                  tag=paste0(model_name,"_tune"),
                  save_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo")))}

## E.4 Save evaluation table
suppressWarnings(nsdm.savethis(object=smev,
              model_name=model_name, species_name=ispi_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d3_evals/glo")))

cat(paste0('\n\nModels fitted and evaluated \n'))

### =========================================================================
### F- Refit top model for prediction
### =========================================================================
## F.1.1 Identify best model
if(model_name!="esm"){
  if(ncol(smev)>1){
    ord<-sort(smev[best_met,],decreasing=T)
    modinp_top<-modinp[names(ord[1])]
  } else {
    modinp_top<-modinp}}
## F.1.2 ... or discard models with best_thre_esm < (for esm)
if(model_name=="esm"){
  ord<-sort(smev[best_met,],decreasing=T)
  ord<-ord[ord>best_thre_esm]
  modinp_top<-modinp[names(ord)]
}

## F.2a Refit model(s) using full dataset
suppressWarnings(prmod<-nsdm.flex3(x=d1_covsels$pseu.abs_i,
                                  replicatetype="none",
                                  reps=1,
								  ncores=1,
								  level="glo",
                                  mod_args=modinp_top,
								  tmp_path=paste0(scr_path,"/tmp/",project)))

suppressWarnings(nsdm.savethis(object=prmod,
                 species_name=ispi_name, model_name=model_name, tag=model_name,
			     compression=TRUE,
                 save_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo")))
			  
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
respcurves<-nsdm.respcurve(prmod, Data=Data,
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
if(length(cov_observ)>0){
cov_obs<-grep(paste0(cov_observ, collapse="|"), names(d1_covsels$covstk), value=T)
} else {
cov_obs<-NULL
}

stk_df<-nsdm.retrieve4pred(covstk=d1_covsels$covstk, 
                           observational=cov_obs,
						   obsval=cov_observ_val,
						   mask=mask_pred, 
                           scaleparam=attributes(d0_datasets$env_vars)[c("scaled:center","scaled:scale")]) 

## H.2 Clean workspace to free some memory before predicting
template<-d1_covsels$covstk[[1]]
rm(d0_datasets, d1_covsels, respcurves, imp, eval_list)
gc()

## H.3 Predict
ndata_bck<-nsdm.predict(models=prmod,
                        nwdata=stk_df$covdf,
                        nsplits=ncores)
						
nsdm.savethis(object=list(ndata_bck=ndata_bck, template=template, nona_ix=stk_df$covdf_ix),
              model_name=model_name, species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d6_preds/glo"))

cat(paste0('Predictions calculated and saved \n'))

cat(paste0('Finished!\n'))