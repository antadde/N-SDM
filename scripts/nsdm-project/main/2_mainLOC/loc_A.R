#############################################################################
## 2_mainLOC: Pipeline for fitting local-level species distribution models
## A: Modelling dataset preparation and covariate selection
## Date: 25-09-2021
## Author: Antoine Adde
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/2_mainLOC","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","outputs",gsub("/main/2_mainLOC","",getwd())),"/settings/nsdm-settings.RData"))

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
### C- species data
### =========================================================================
# Target species
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))
ispi_name<-species[arrayID]

# Load species data
sp_dat<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/glo/",ispi_name,"/",ispi_name,".rds"))
pseu.abs_i<-sp_dat$pseu.abs_i_loc
group_name<-gsub("\\s*\\([^\\)]+\\)","",as.character(sp_dat$group))

cat(paste0('Ready for modelling dataset preparation and covariate selection for ', ispi_name, '...\n'))

### =========================================================================
### D- Covariate data
### =========================================================================
# D.1 Retrieve loc list of candidate covariates and covinfo table
lr<-readRDS(paste0(w_path,"outputs/",project,"/settings/covariates-list.rds"))
lr_loc<-lr$lr_loc[!grepl("/future/", lr$lr_loc) & !grepl("/bioclim/", lr$lr_loc)]
cov_info<-data.frame(lr$cov_info[match(lr_loc, lr$cov_info$file),])

# D.2 Subset with list of expert-filtered candidate covariates, if available
expert_tab<-try(read_excel(expert_table, .name_repair = "minimal"), silent=T)
if(exists("expert_tab")){
colnames(expert_tab)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(expert_tab)))
dup_ix<-which(duplicated(colnames(expert_tab)))
if(length(dup_ix)>0) expert_tab<-expert_tab[,-which(duplicated(colnames(expert_tab)))]
expert_tab<-expert_tab[which(expert_tab[,group_name]=="1"),]
cov_info<-merge(cov_info, expert_tab)
lr_loc<-cov_info$file
}

# Retrieve glo and loc reference rasters
rsts_ref<-readRDS(paste0(w_path,"outputs/",project,"/settings/ref-rasters.rds"))

cat(paste0('Covariate data listed \n'))

### =========================================================================
### E- Prepare modelling dataset
### =========================================================================
# E.1 Retrieve predictions from mainGLO
glo_out<-list.files(paste0(scr_path,"/outputs/",project,"/d8_ensembles/glo/",ispi_name), pattern=".rds", full.names = TRUE)

# E.2 Extract new habitat covariates
pseu.abs_i<-nsdm.bigextract(cov=c(gsub(".rds", ".fst", lr_loc), glo_out),
                            data=pseu.abs_i,
							rst_ref=rsts_ref$rst_loc,
							cov_info=cov_info,
							t_match=tmatch_loc,
							tmatch_scheme=tmatch_scheme_loc,
							p_int=pint_loc,
							ex_pint=FALSE,
							nsplits=ncores)

# Update cov_info table
cov_info<-na.omit(cov_info[match(colnames(pseu.abs_i@env_vars), gsub(".rds","",basename(cov_info$file))),])

# Scale env_vars							
colnames(pseu.abs_i@env_vars)[ncol(pseu.abs_i@env_vars)]<-"mainGLO"
env_vars<-scale(pseu.abs_i@env_vars) # keep scaling parameters for backtransforming predictions later
pseu.abs_i@env_vars<-data.frame(env_vars)

# E.3 Define weights
wi<-which(pseu.abs_i@pa==1)
wt<-rep(1,length(pseu.abs_i@pa))
wt[wi]<-round((length(pseu.abs_i@pa)-length(wi))/length(wi))
if(unique(wt[wi]==0)) wt[wi]<-1

# E.4 Save modelling set
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i,
                          weights=wt,
						  cov_info=cov_info,
						  env_vars=env_vars,
						  group=sp_dat$group),
                          species_name=ispi_name,
			  			  compression=TRUE,
                          save_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/loc"))

cat(paste0('Modelling dataset prepared \n'))

### =========================================================================
### F- Covariate selection with mainGLO forced
### =========================================================================
covstk_res<-list()
covdata_res<-list()

if("covariate" %in% nesting_methods){
counter <- 0
while(TRUE){
counter <- sum(counter, 1)

# F.1 Step 1: Filtering for colinearity
cat('covariate selection with mainGLO forced S1: Filtering...\n')
cov.filter_i<-try(nsdm.filtersel(pa=pseu.abs_i@pa, # pa vector
                             covdata=pseu.abs_i@env_vars, # data.frame of environmental covariates extracted at pa
                             weights=wt, # weight vector
                             datasets=c(cov_info$cada, "mainGLO"),  
							 variables=c(gsub("_NA", "", paste(cov_info$variable, cov_info$attribute, sep="_")), "mainGLO"), 
                             focals=unique(cov_info[which(cov_info$focal!="NA"),"cada"]), # datasets with focal window selection
                             method=sel_met, # univariate ranking method to be used
                             corcut=cor_cut,
							 force=c("mainGLO")), silent=TRUE) # correlation cutoff for colinearity

# F.2 Step 2: Model-specific embedding
cat('covariate selection with mainGLO forced S2: Embedding...\n')
cov.embed_i<-try(nsdm.embedsel(pa=pseu.abs_i@pa,
                           covdata=cov.filter_i, # filtered covariate set from Step 1
                           weights=wt,
						   force=c("mainGLO"),
                           nthreads=ncores), silent=TRUE)
						   
# F.3 Step 3: Overall ranking
cat('covariate selection with mainGLO forced S3: Ranking...\n')
cov.rk_i<-try(nsdm.covselrk(embed=cov.embed_i, # embedded results (S2)
                        species_name=ispi_name), silent=TRUE)
					
cat('covariate selection with mainGLO forced S4: Final subsetting...\n')
# F.4 Step 4: Final covariate subset
hab_stk_loc<-try(nsdm.fastraster(files=na.omit(cov_info$file[match(cov.rk_i$var, gsub(".rds","",basename(cov_info$file)))])[1:max_thre], nsplits=ncores), silent=TRUE)
cov.sub_i_cov <-try(nsdm.covsub(covdata=pseu.abs_i@env_vars,
            rasterdata=hab_stk_loc,
            ranks=cov.rk_i[-(which(cov.rk_i$var=="mainGLO")),], # ranking from S3
            thre=ceiling(log2(table(pseu.abs_i@pa)['1']))-1, # target number of covariates
			max.thre=max_thre-1, # max number of possible covariates in model
			glo.out=readRDS(glo_out), # add mainGLO outputs to covariate subset
			glo.xy=pseu.abs_i@xy, #xy coordinates where to extract glo.out
			), silent=TRUE)
			
covstk_res[["cov"]]<-cov.sub_i_cov$rasterdata
covdata_res[["cov"]]<-cov.sub_i_cov$covdata

if(counter > 5 | !is(cov.sub_i_cov, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}}

### =========================================================================
### G- Covariate selection without main GLO
### =========================================================================
if("multiply" %in% nesting_methods){
counter <- 0
while(TRUE){
counter <- sum(counter, 1)
# Remove mainGLO
pseu.abs_i@env_vars<-subset(pseu.abs_i@env_vars, select=-c(mainGLO))

# G.1 Step 1: Filtering for colinearity
cat('covariate selection without mainGLO S1: Filtering...\n')
cov.filter_i<-try(nsdm.filtersel(pa=pseu.abs_i@pa, # pa vector
                             covdata=pseu.abs_i@env_vars, # data.frame of environmental covariates extracted at pa
                             weights=wt, # weight vector
                             datasets=cov_info$cada,  
							 variables=gsub("_NA", "", paste(cov_info$variable, cov_info$attribute, sep="_")), 
                             focals=unique(cov_info[which(cov_info$focal!="NA"),"category"]), # categories with focal window selection
                             method=sel_met, # univariate ranking method to be used
                             corcut=cor_cut), silent=TRUE) # correlation cutoff for colinearity

# G.2 Step 2: Model-specific embedding
cat('covariate selection without mainGLO S2: Embedding...\n')
cov.embed_i<-try(nsdm.embedsel(pa=pseu.abs_i@pa,
                           covdata=cov.filter_i, # filtered covariate set from Step 1
                           weights=wt,
                           nthreads=ncores), silent=TRUE)

# G.3 Step 3: Overall ranking
cat('covariate selection without mainGLO S3: Ranking...\n')
cov.rk_i<-try(nsdm.covselrk(embed=cov.embed_i, # embedded results (S2)
                        species_name=ispi_name), silent=TRUE)
		
cat('covariate selection without mainGLO S4: Final subsetting...\n')
# G.4 Step 4: Final covariate subset
hab_stk_loc<-try(nsdm.fastraster(files=na.omit(cov_info$file[match(cov.rk_i$var, gsub(".rds","",basename(cov_info$file)))])[1:max_thre], nsplits=ncores), silent=TRUE)
cov.sub_i_mul <-try(nsdm.covsub(covdata=pseu.abs_i@env_vars,
            rasterdata=hab_stk_loc,
            ranks=cov.rk_i, # ranking from S3
            thre=ceiling(log2(table(pseu.abs_i@pa)['1'])), # target number of covariates
			max.thre=max_thre # max number of possible covariates in model
			), silent=TRUE) # buffer extraction radius for glo.out
			
covstk_res[["mul"]]<-cov.sub_i_mul$rasterdata
covdata_res[["mul"]]<-cov.sub_i_mul$covdata

if(counter > 5 | !is(cov.sub_i_mul, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}}

# H.5 Save covariate selection results
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i,
              covstk=covstk_res,
			  covdata=covdata_res,
			  env_vars=env_vars),
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/loc"))

cat(paste0('Covariate selection for ',ispi_name,' done \n'))
cat(paste0('Finished \n'))