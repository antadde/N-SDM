#############################################################################
## 0_mainPRE
## A: read N-SDM settings; load input data; occurrences disaggregation; covariate formatting
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Main N-SDM settings
### =========================================================================
# Project name
project<-gsub("/main","",gsub(".*scripts/","",getwd()))

# Set permissions for new files
Sys.umask(mode="000")

# Load and retrieve main settings
settings<-read.csv2("./settings/settings.csv")
parameters<-settings$parameter
values<-settings$value
for(i in 1:length(parameters)){
if(grepl(pattern="NULL", values[i])){
do.call("<-",list(parameters[i], NULL))
} else {
if(grepl(pattern=paste0(c("c\\('|paste", ":20"), collapse="|"), values[i])){
do.call("<-",list(parameters[i], eval(parse(text=values[i]))))
} else {
 if(is.na(as.numeric(values[i]))) {
do.call("<-",list(parameters[i], values[i]))
 } else {
do.call("<-",list(parameters[i], as.numeric(values[i])))
}
}
}
}

# Additional settings
ssl_id<-readLines(paste0(w_path,"tmp/",project,"/settings/tmp/ssl_id.txt"))
cov_path<-paste0(w_path,"data/",project,"/covariates/")
spe_glo<-list.files(paste0(w_path,"data/",project,"/species/glo"), full.names=T, pattern=".rds")
if(n_levels>1) spe_reg<-list.files(paste0(w_path,"data/",project,"/species/reg"), full.names=T, pattern=".rds")
param_grid<-paste0(w_path,"scripts/",project,"/main/settings/", param_grid)
if(length(expert_table)>0) expert_table<-paste0(w_path,"scripts/",project,"/main/settings/", expert_table)
if(length(forced_species)>0) forced_species<-paste0(w_path,"scripts/",project,"/main/settings/", forced_species)

# Check and refine masks
if(n_levels>1) mask_reg<-paste0(w_path,"data/",project,"/masks/", mask_reg)
if(length(mask_pred)>0) mask_pred<-paste0(w_path,"data/",project,"/masks/", mask_pred)

# Save settings
rm(settings, parameters, values, i)
save.image(paste0(w_path,"tmp/",project,"/settings/nsdm-settings.RData"))

print(paste0("N-SDM settings defined"))

### =========================================================================
### B- Prepare covariate data
### =========================================================================
# Set lib path
.libPaths(lib_path)

# Load nsdm package
require(nsdm)

ncores<-as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

# Generate covariate info table
nsdm.covinfo(cov_path=cov_path, save_path=paste0(w_path,"scripts/",project,"/main/settings/"),
             time_cov=cov_time, focal_cov=cov_focal)
			 
cov_info<-read_excel(paste0(w_path,"scripts/", project, "/main/settings/predictors-available.xlsx"))
			 
# Rewrite glo rasters masked with mask_reg if not already done
lr_glo<-cov_info$file[cov_info$level=="glo"]
if(n_levels>1){
if(length(grep("msk", lr_glo)) != length(lr_glo)/2){
msk_reg<-readRDS(mask_reg)
for(m in lr_glo){
r<-readRDS(m)
r_msk<-raster::mask(r, msk_reg, inverse=TRUE)
saveRDS(r_msk, gsub(".rds", "_msk.rds", m, fixed=T))
nsdm.covinfo(cov_path=cov_path, save_path=paste0(w_path,"scripts/",project,"/main/settings/"),
             time_cov=cov_time, focal_cov=cov_focal)			 
cov_info<-read_excel(paste0(w_path,"scripts/", project, "/main/settings/predictors-available.xlsx"))
}
}
}

# list availabe covariates
if(n_levels>1){
lr_glo<-cov_info$file[cov_info$level=="glo" & cov_info$attribute=="msk"]
} else {
lr_glo<-cov_info$file[cov_info$level=="glo"]
}
if(n_levels>1){lr_reg<-cov_info$file[cov_info$level=="reg"]
lr<-c(lr_glo, lr_reg)
} else {
lr<-lr_glo}

# Create fst versions 
threads_fst(nr_of_threads = 1)

# if not already done: check if fst version exists
k<-sapply(lr, function(f){file.exists(gsub(".rds", ".fst", f, fixed=T))})
lr<-lr[!k]

# otherwise do it now
if(length(lr) > 0){
pp<-mclapply(lr, function(c){
r<-readRDS(c)
r_df<-as.data.frame(r)
write.fst(r_df, gsub(".rds", ".fst", c, fixed=T), 75)
names(r)
}, mc.cores=ncores)}

# Create regional and global reference rasters
## for glo if bioclim layer available use it; else first one
if(length(grep("bio1", lr_glo, value=T))>0){
rst_glo<-readRDS(grep("bio1", lr_glo, value=T)[1])
}else{
rst_glo<-readRDS(lr_glo[1])}

## for reg if bioclim layer available use it; else first one
if(n_levels>1){
if(length(grep("bio1", lr_reg, value=T))>0){
rst_reg<-readRDS(grep("bio1", lr_reg, value=T)[1])
}else{
rst_reg<-readRDS(lr_reg[1])}
}

# Save covariate data settings
## reference rasters
if(n_levels>1){
l<-list(rst_reg=rst_reg, rst_glo=rst_glo)
} else {
l<-list(rst_glo=rst_glo)}

saveRDS(l,  
        paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))
		
## covariates list and info
if(n_levels>1){
l2<-list(lr_reg=lr_reg, lr_glo=lr_glo, cov_info=cov_info)
} else {
l2<-list(lr_glo=lr_glo, cov_info=cov_info)}

saveRDS(l2, 
        paste0(w_path,"tmp/",project,"/settings/covariates-list.rds"))

print(paste0("Covariate settings defined"))

### =========================================================================
### C- Refine the list of species to be modelled
### =========================================================================
# Check if .fst versions of regional and global species data exist or create them
threads_fst(nr_of_threads = ncores)

if(n_levels>1){
if(!file.exists(gsub(".rds", ".fst", spe_reg))){
spe_reg_pts<-readRDS(spe_reg)
if(class(spe_reg_pts)=="SpatialPointsDataFrame"){spe_reg_pts_dat<-spe_reg_pts@data
} else {
spe_reg_pts_dat<-spe_reg_pts}
write_fst(spe_reg_pts_dat, gsub(".rds", ".fst", spe_reg), 0)
}
}

if(!file.exists(gsub(".rds", ".fst", spe_glo))){
spe_glo_pts<-readRDS(spe_glo)
if(class(spe_glo_pts)=="SpatialPointsDataFrame"){spe_glo_pts_dat<-spe_glo_pts@data
} else {
spe_glo_pts_dat<-spe_glo_pts}
write_fst(spe_glo_pts_dat, gsub(".rds", ".fst", spe_glo), 0)
}

# Load regional and global species fst data
if(n_levels>1){
spe_reg_fst <-read_fst(gsub(".rds", ".fst", spe_reg))
if("canonicalname" %in% names(spe_reg_fst)) spe_reg_fst$species<-spe_reg_fst$canonicalname
spe_reg_names<-unique(spe_reg_fst$species)
print(paste0("Initial number of species in REG species dataset is: ",length(spe_reg_names)))
}

spe_glo_fst<-read_fst(gsub(".rds", ".fst", spe_glo))
if("canonicalname" %in% names(spe_glo_fst)) spe_glo_fst$species<-spe_glo_fst$canonicalname
spe_glo_names<-unique(spe_glo_fst$species)
print(paste0("Initial number of species in GLO species dataset is: ",length(spe_glo_names)))

# Intersect global and regional lists of species (matching)
if(n_levels>1){
species<-sort(intersect(spe_glo_names, spe_reg_names))
} else {
species<-spe_glo_names}

print(paste0("Number of remaining species after intersecting GLO and REG species dataset is: ",length(species)))

# If a vector of species to be forced is available use it
if(length(forced_species)>0) forced_species<-read.csv2(forced_species, header=FALSE)[,1]
if(is.character(forced_species)){
species<-forced_species
}

# discard species to match n_spe (number of species to be modelled) if needed
if(n_spe<length(species)){
species<-species[1:n_spe]
}

print(paste0("Total number of species considered for this N-SDM run is: ",length(species)))

### =========================================================================
### D- Spatiotemporal disaggregation of species data
### =========================================================================
# Subset regional species data
if(n_levels>1){

spe_reg_pts<-spe_reg_fst[spe_reg_fst$species %in% species,]

print(paste0("Ready for spatiotemporal disaggregation of regional data for ",length(species)," species"))

# Disaggregate regional species data
if(tmatch_reg==FALSE) spe_reg_pts$year<-NULL # no temporal matching
if(disag_reg==TRUE){
spe_reg_dis<-nsdm.disaggregate(pres=spe_reg_pts, rst=rst_reg, thindist=thin_dist_reg, thinyear=thin_time_reg, max_uncertain=max_uncertain_reg, min_occ=min_occ_reg, ncores=ncores)
} else { # if no disaggregation is requested
spe_reg_dis<-spe_reg_pts
}

# Update species list
species<-unique(spe_reg_dis$species)
}

# Subset global species data
spe_glo_pts<-spe_glo_fst[spe_glo_fst$species %in% species,]

print(paste0("Ready for spatiotemporal disaggregation of global data for ",length(species)," species"))

# disaggregate global species data
if(tmatch_glo==FALSE) spe_glo_pts$year<-NULL # no temporal matching
if(disag_glo==TRUE){
spe_glo_dis<-nsdm.disaggregate(pres=spe_glo_pts, rst=rst_glo, thindist=thin_dist_glo, thinyear=thin_time_glo, max_uncertain=max_uncertain_glo, min_occ=min_occ_glo, ncores=ncores)
} else { # if no disaggregation is requested
spe_glo_dis<-spe_glo_pts
}

# Save species data settings
## datasets
if(n_levels>1){
l3<-list(spe_reg_dis=spe_reg_dis, spe_glo_dis=spe_glo_dis)
} else {
l3<-list(spe_glo_dis=spe_glo_dis)}

saveRDS(l3,
        paste0(w_path,"tmp/",project,"/settings/species-occurences.rds"))
		
## species list
saveRDS(species,  
        paste0(w_path,"tmp/",project,"/settings/species-list.rds"))
		
## number of cluster runs needed
spe_runs<-length(splitIndices(length(species), length(species)/n_mx_spe))
writeLines(as.character(spe_runs),        paste0(w_path,"tmp/",project,"/settings/tmp/spe_runs.txt"))
writeLines(as.character(length(species)), paste0(w_path,"tmp/",project,"/settings/tmp/n_spe.txt"))

print(paste0("Species data and settings defined"))

print("Finished!")