#############################################################################
## 2_mainREG
## A: covariate extraction and covariate selection
## Date: 20-05-2022
## Author: Antoine Adde
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/2_mainREG","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","tmp",gsub("/main/2_mainREG","",getwd())),"/settings/nsdm-settings.RData"))

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
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))
ispi_name<-species[arrayID]

# Load species data
sp_dat<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/glo/",ispi_name,"/",ispi_name,".rds"))
pseu.abs_i<-sp_dat$pseu.abs_i_reg
group_name<-gsub("\\s*\\([^\\)]+\\)","",as.character(sp_dat$group))

cat(paste0('Ready for modelling dataset preparation and covariate selection for ', ispi_name, '...\n'))

### =========================================================================
### D- Covariate data
### =========================================================================
# D.1 Retrieve reg list of candidate covariates and covinfo table
lr<-readRDS(paste0(w_path,"tmp/",project,"/settings/covariates-list.rds"))
cov_info<-lr$cov_info
lr_reg<-lr$lr_reg[intersect(grep("/future/", lr$lr_reg, invert=TRUE), grep(paste0("/",unique(cov_info$category[cov_info$level=="glo"]),"/"), lr$lr_reg, invert=T))]
cov_info<-data.frame(lr$cov_info[match(lr_reg, lr$cov_info$file),])

# D.2 Subset with list of expert-filtered candidate covariates, if available
expert_tab<-try(read_excel(expert_table, .name_repair = "minimal"), silent=T)
if(!"try-error" %in% class(expert_tab)){
colnames(expert_tab)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(expert_tab)))
dup_ix<-which(duplicated(colnames(expert_tab)))
if(length(dup_ix)>0) expert_tab<-expert_tab[,-which(duplicated(colnames(expert_tab)))]
expert_tab<-expert_tab[which(expert_tab[,group_name]=="1"),]
cov_info<-merge(cov_info, expert_tab)
lr_reg<-cov_info$file
}

# Retrieve glo and reg reference rasters
rsts_ref<-readRDS(paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))

cat(paste0('Covariate data listed \n'))

### =========================================================================
### E- Prepare modelling dataset
### =========================================================================
# E.1 Retrieve predictions from mainGLO
glo_out_f<-list.files(paste0(scr_path,"/outputs/",project,"/d8_ensembles/glo/",ispi_name), pattern=".rds", full.names = TRUE)

# E.2 Extract new habitat covariates
pseu.abs_i<-nsdm.bigextract(cov=c(gsub(".rds", ".fst", lr_reg), glo_out_f),
                            data=pseu.abs_i,
							rst_ref=rsts_ref$rst_reg,
							cov_info=cov_info,
							t_match=tmatch_reg,
							tmatch_scheme=tmatch_scheme_reg,
                                                                                    nzvt=16,
							ex_pint=FALSE,
							nsplits=ncores)

# Update cov_info table
cov_info<-na.omit(cov_info[match(colnames(pseu.abs_i@env_vars), gsub(".rds","",basename(cov_info$file))),])

# Scale env_vars							
colnames(pseu.abs_i@env_vars)[ncol(pseu.abs_i@env_vars)]<-"mainGLO"
env_vars<-scale(pseu.abs_i@env_vars)
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
                          save_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/reg"))

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

# F.1 Step 1: Filtering for collinearity
cat('covariate selection with mainGLO forced S1: Filtering...\n')
cov.filter_i<-try(covsel.filter(
  pseu.abs_i@env_vars,
  pseu.abs_i@pa,
  variables = c(cov_info$variable, "mainGLO"),
  categories = c(cov_info$cada, "mainGLO"),
  weights = wt,
  force = c("mainGLO"),
  corcut = cor_cut), silent=TRUE)

# F.2 Step 2: Model-specific embedding
cat('covariate selection with mainGLO forced S2: Embedding...\n')
cov.embed_i<-try(covsel.embed(
  cov.filter_i,
  pseu.abs_i@pa,
  weights = wt,
  algorithms = c("glm", "gam", "rf"),
  force = c("mainGLO"),
  ncov = if(!"esm" %in% mod_algo){ceiling(log2(table(pseu.abs_i@pa)['1']))-1} else {ncov_esm},
  maxncov  = max_thre,
  nthreads = ncores), silent=TRUE)
						   
if(counter > 5 | !is(cov.embed_i, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}
# F.3 Finalize
hab_stk_reg<-try(nsdm.fastraster(files=na.omit(cov_info$file[match(cov.embed_i$ranks_2$covariate, gsub(".rds","",basename(cov_info$file)))]), nsplits=ncores), silent=TRUE)
glo_out<-readRDS(glo_out_f); names(glo_out)<-"mainGLO"
covstk_res[["cov"]]<-stack(hab_stk_reg, glo_out)
covdata_res[["cov"]]<-cov.embed_i$covdata
}

### =========================================================================
### G- Covariate selection without main GLO
### =========================================================================
if("multiply" %in% nesting_methods){
counter <- 0
while(TRUE){
counter <- sum(counter, 1)

# Remove mainGLO
pseu.abs_i@env_vars<-subset(pseu.abs_i@env_vars, select=-c(mainGLO))

# G.1 Step 1: Filtering for collinearity
cat('covariate selection S1: Filtering...\n')
cov.filter_i<-try(covsel.filter(
  pseu.abs_i@env_vars,
  pseu.abs_i@pa,
  variables = cov_info$variable,
  categories = cov_info$cada,
  weights = wt,
  corcut = cor_cut), silent=TRUE)

# G.2 Step 2: Model-specific embedding
cat('covariate selection S2: Embedding...\n')
cov.embed_i<-try(covsel.embed(
  cov.filter_i,
  pseu.abs_i@pa,
  weights = wt,
  algorithms = c("glm", "gam", "rf"),
  ncov = if(!"esm" %in% mod_algo){ceiling(log2(table(pseu.abs_i@pa)['1']))-1} else {ncov_esm},
  maxncov  = max_thre,
  force=NULL,
  nthreads = ncores), silent=TRUE)
						   
if(counter > 5 | !is(cov.embed_i, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}
# G.3 Finalize
covstk_res[["mul"]]<-try(nsdm.fastraster(files=na.omit(cov_info$file[match(cov.embed_i$ranks_2$covariate, gsub(".rds","",basename(cov_info$file)))]), nsplits=ncores), silent=TRUE)
covdata_res[["mul"]]<-try(cov.embed_i$covdata, silent=TRUE)
}

# Save
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i,
              covstk=covstk_res,
			  covdata=covdata_res,
			  env_vars=env_vars),
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/reg"))

cat(paste0('Covariate selection for ',ispi_name,' done \n'))
cat(paste0('Finished \n'))