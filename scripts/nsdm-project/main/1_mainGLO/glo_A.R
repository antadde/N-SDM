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

### =========================================================================
### C- Species data
### =========================================================================
# Target species arrayID
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))
ispi_name<-species[arrayID]

# Load species data
sp_dat<-readRDS(paste0(scr_path,"/outputs/",project,"/d0_datasets/base/",ispi_name,"/",ispi_name,".rds"))
group_name<-sp_dat$group

cat(paste0('Ready for modelling dataset preparation and covariate selection for ', ispi_name, '...\n'))

### =========================================================================
### D- Covariate data
### =========================================================================
# Retrieve lists of candidate covariates and covinfo table
lr<-readRDS(paste0(w_path,"tmp/",project,"/settings/covariates-list.rds"))
cov_info<-lr$cov_info
# Refine glo set
lr_glo<-lr$lr_glo[grep("/future/", lr$lr_glo, invert=TRUE)]
cov_info_glo<-data.frame(lr$cov_info[match(lr_glo, lr$cov_info$file),])
# Refine reg set
if(n_levels>1){
lr_reg<-lr$lr_reg[intersect(grep("/future/", lr$lr_reg, invert=TRUE), grep(paste0("/",unique(cov_info_glo$category),"/"), lr$lr_reg))]
cov_info_reg<-data.frame(lr$cov_info[match(lr_reg, lr$cov_info$file),])
}


# Subset with list of expert-filtered candidate covariates, if available
expert_tab<-try(read_excel(expert_table, .name_repair = "minimal"), silent=T)
if(!"try-error" %in% class(expert_tab)){
colnames(expert_tab)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(expert_tab)))
dup_ix<-which(duplicated(colnames(expert_tab)))
if(length(dup_ix)>0) expert_tab<-expert_tab[,-which(duplicated(colnames(expert_tab)))]
expert_tab<-data.frame(expert_tab[which(expert_tab[,group_name]=="1"),])
## glo
cov_info_glo<-merge(cov_info_glo, expert_tab)
lr_glo<-cov_info_glo$file
if(n_levels>1){
## reg
cov_info_reg<-merge(cov_info_reg, expert_tab)
lr_reg<-cov_info_reg$file
}
}

# Retrieve glo and reg reference rasters
rsts_ref<-readRDS(paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))

### =========================================================================
### E- Covariate extraction
### =========================================================================
# E.1 GLO
pseu.abs_i_glo<-nsdm.bigextract(cov=gsub(".rds", ".fst", lr_glo),
                                data=sp_dat$pseu.abs_i_glo,
							    rst_ref=rsts_ref$rst_glo,
							    cov_info=cov_info_glo,
							    t_match=tmatch_glo,
							    nsplits=ncores)

pseu.abs_i_glo_copy<-pseu.abs_i_glo
                         
# E.2 REG
if(n_levels>1){
pseu.abs_i_reg<-nsdm.bigextract(cov=gsub(".rds", ".fst", lr_reg),
                               data=sp_dat$pseu.abs_i_reg,
							   rst_ref=rsts_ref$rst_reg,
							   cov_info=cov_info_reg,
							   t_match=tmatch_reg,
							   nsplits=ncores)
							   
# E.3 Combine GLO and REG data
pseu.abs_i_glo@pa<-c(pseu.abs_i_glo@pa, pseu.abs_i_reg@pa)
pseu.abs_i_glo@years<-c(pseu.abs_i_glo@years, pseu.abs_i_reg@years)
pseu.abs_i_glo@xy<-rbind(pseu.abs_i_glo@xy, pseu.abs_i_reg@xy)
m_ord<-match(cov_info_reg$variable, cov_info_glo$variable)
names(pseu.abs_i_glo@env_vars)[m_ord]<-names(pseu.abs_i_reg@env_vars) 
env_vars<-scale(rbind(pseu.abs_i_glo@env_vars, pseu.abs_i_reg@env_vars))
pseu.abs_i_glo@env_vars<-data.frame(env_vars)							   
} else {
env_vars<-scale(pseu.abs_i_glo@env_vars)
pseu.abs_i_glo@env_vars<-data.frame(env_vars)
# E.4 Update cov_info table
cov_info_glo<-na.omit(cov_info_glo[match(colnames(pseu.abs_i_glo@env_vars), gsub(".rds","",basename(cov_info_glo$file))),])					   
}

# E.5 Define weights
wi<-which(pseu.abs_i_glo@pa==1)
wt<-rep(1,length(pseu.abs_i_glo@pa))
wt[wi]<-round((length(pseu.abs_i_glo@pa)-length(wi))/length(wi))
if(unique(wt[wi]==0)) wt[wi]<-1

# E.6 Save modelling set
if(n_levels>1){
l<-list(group=sp_dat$group,
        pseu.abs_i_glo_copy=pseu.abs_i_glo_copy,
        pseu.abs_i_glo=pseu.abs_i_glo,
        pseu.abs_i_reg=pseu.abs_i_reg,
		weights=wt,
		env_vars=env_vars)
} else {
l<-list(group=sp_dat$group,
        pseu.abs_i_glo_copy=pseu.abs_i_glo_copy,
        pseu.abs_i_glo=pseu.abs_i_glo,
		weights=wt,
		env_vars=env_vars)
}

nsdm.savethis(l,
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
cov.filter_i<-try(covsel.filter(
  pseu.abs_i_glo@env_vars,
  pseu.abs_i_glo@pa,
  variables = cov_info_glo$variable,
  categories = cov_info_glo$cada,
  weights = wt,
  corcut = cor_cut), silent=TRUE)

# F.2 Step 2: Model-specific embedding
cat('covariate selection S2: Embedding...\n')
cov.embed_i<-try(covsel.embed(
  cov.filter_i,
  pseu.abs_i_glo@pa,
  weights = wt,
  ncov = if(!"esm" %in% mod_algo){ceiling(log2(table(pseu.abs_i_glo@pa)['1']))-1} else {ncov_esm},
  maxncov  = max_thre,
  nthreads = ncores), silent=TRUE)
  
if(counter > 5 | !is(cov.embed_i, 'try-error')) break
cat(paste0("an error occured; tentative number ", counter, " for covariate selection...\n"))
}

# F.3 Finalize
if(n_levels>1){
stk<-try(nsdm.fastraster(files=na.omit(cov_info_reg$file[match(cov.embed_i$ranks_2$covariate, gsub(".rds","",basename(cov_info_reg$file)))]), nsplits=ncores), silent=TRUE)
} else {
stk<-try(nsdm.fastraster(files=na.omit(cov_info_glo$file[match(cov.embed_i$ranks_2$covariate, gsub(".rds","",basename(cov_info_glo$file)))]), nsplits=ncores), silent=TRUE)}


# F.4 Save
pseu.abs_i_glo@env_vars<-cov.embed_i$covdata
nsdm.savethis(object=list(pseu.abs_i=pseu.abs_i_glo, filter=names(cov.filter_i), embed=cov.embed_i, covstk=stk),
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/glo"))

cat(paste0('Covariate selection done \n'))

cat(paste0('Finished \n'))