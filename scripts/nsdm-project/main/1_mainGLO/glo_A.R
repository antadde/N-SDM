#############################################################################
## 1_mainGLO
## A: covariate extraction and covariate selection
## Date: 20-05-2022
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

### =========================================================================
### C- Covariate data
### =========================================================================
# Retrieve glo and loc lists of candidate covariates and covinfo table
lr<-readRDS(paste0(w_path,"outputs/",project,"/settings/covariates-list.rds"))
lr_glo<-lr$lr_glo[intersect(grep("/present/", lr$lr_glo), grep("/bioclim/", lr$lr_glo))]
lr_loc<-lr$lr_loc[intersect(grep("/present/", lr$lr_loc), grep("/bioclim/", lr$lr_loc))]
cov_info<-lr$cov_info

# Retrieve glo and loc reference rasters
rsts_ref<-readRDS(paste0(w_path,"outputs/",project,"/settings/ref-rasters.rds"))

### =========================================================================
### D- Species data
### =========================================================================
# Target species arrayID
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))
ispi_name<-species[arrayID]

# Load species data
sp_dat<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/base/",ispi_name,"/",ispi_name,".rds"))

cat(paste0('Ready for modelling dataset preparation and covariate selection for ', ispi_name, '...\n'))

### =========================================================================
### E- Covariate extraction
### =========================================================================
# E.1 GLO
pseu.abs_i_glo<-nsdm.bigextract(cov=gsub(".rds", ".fst", lr_glo),
                                data=sp_dat$pseu.abs_i_glo,
							    rst_ref=rsts_ref$rst_glo,
							    cov_info=cov_info,
							    t_match=tmatch_glo,
							    nsplits=ncores)

pseu.abs_i_glo_copy<-pseu.abs_i_glo
                         
# E.2 LOC
pseu.abs_i_loc<-nsdm.bigextract(cov=gsub(".rds", ".fst", lr_loc),
                               data=sp_dat$pseu.abs_i_loc,
							   rst_ref=rsts_ref$rst_loc,
							   cov_info=cov_info,
							   t_match=FALSE,
							   p_int=pint_glo,
							   nsplits=ncores)
				   							 
# E.3 Combine GLO and LOC data
pa<-c(pseu.abs_i_glo@pa, pseu.abs_i_loc@pa)
pseu.abs_i_glo@pa<-pa
env_vars<-scale(rbind(pseu.abs_i_glo@env_vars, pseu.abs_i_loc@env_vars)) # keep scaling parameters to backtransform predictions later
pseu.abs_i_glo@env_vars<-data.frame(env_vars)
pseu.abs_i_glo@xy<-matrix() # empty xy matrix for safety (mix of 2 coord systems)

# E.4 Define weights
wi<-which(pseu.abs_i_glo@pa==1)
wt<-rep(1,length(pseu.abs_i_glo@pa))
wt[wi]<-round((length(pseu.abs_i_glo@pa)-length(wi))/length(wi))
if(unique(wt[wi]==0)) wt[wi]<-1

# E.5 Save modelling set
nsdm.savethis(object=list(group=sp_dat$group,
                          pseu.abs_i_glo_copy=pseu.abs_i_glo_copy,
                          pseu.abs_i_glo=pseu.abs_i_glo,
                          pseu.abs_i_loc=pseu.abs_i_loc,
						  weights=wt,
						  env_vars=env_vars),
              species_name=ispi_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/glo"))
			  
cat(paste0('Modelling dataset prepared \n'))

### =========================================================================
### F- Covariate selection
### =========================================================================
counter <- 0
while(TRUE){
counter <- sum(counter, 1)

# F.1 Step 1: Filtering for colinearity
cat('covariate selection S1: Filtering...\n')
cov.filter_i<-try(nsdm.filtersel(pa=pseu.abs_i_glo@pa, # pa vector
                             covdata=pseu.abs_i_glo@env_vars, # data.frame of environmental covariates extracted at pa
                             weights=wt, # weight vector
                             datasets=rep("bioclim", length(pseu.abs_i_glo@env_vars)),
                             method=sel_met, # univariate ranking method to be used
                             corcut=cor_cut), silent=TRUE) # correlation cutoff for colinearity

# F.2 Step 2: Model-specific embedding
cat('covariate selection S2: Embedding...\n')
cov.embed_i<-try(nsdm.embedsel(pa=pseu.abs_i_glo@pa,
                           covdata=cov.filter_i, # filtered covariate set from Step 1
                           weights=wt,
                           nthreads=ncores), silent=TRUE)

# F.3 Step 3: Overall ranking
cat('covariate selection S3: Ranking...\n')
cov.rk_i<-try(nsdm.covselrk(embed=cov.embed_i, # embedded results (S2)
                        species_name=ispi_name), silent=TRUE)

cat('covariate selection S4: Final subsetting...\n')
# F.4 Step 4: Final covariate subset
clim_stk_loc<-nsdm.fastraster(grep(paste0("_",min(pint_glo),"_",max(pint_glo),"_"), lr_loc, value=T), ncores)
names(clim_stk_loc)<-gsub(".*_", "", names(clim_stk_loc))
cov.sub_i<-try(nsdm.covsub(covdata=pseu.abs_i_glo@env_vars,
            rasterdata=clim_stk_loc,
            ranks=cov.rk_i, # ranking from S3
            thre=ceiling(log2(table(pseu.abs_i_glo@pa)['1'])), # target number of covariates (log2 nobs)
			max.thre=max_thre), silent=TRUE) # max number of possible covariates in model

if(counter > 5 | !is(cov.sub_i, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}
			
# update with cov.sub outputs
clim_stk_loc<-cov.sub_i$rasterdata
pseu.abs_i_glo@env_vars<-cov.sub_i$covdata

# F.5 Save covariate selection results
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i_glo, filter=names(cov.filter_i), embed=cov.embed_i, ranking=cov.rk_i,covstk=clim_stk_loc),
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/glo"))

cat(paste0('Covariate selection done \n'))

cat(paste0('Finished \n'))