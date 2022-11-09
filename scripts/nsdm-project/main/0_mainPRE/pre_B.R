#############################################################################
## 0_mainPRE
## B: pseudo-absences; modelling object generation
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################

### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/0_mainPRE","",gsub(".*scripts/","",getwd()))

# Load N-SDM settings
load(paste0(gsub("scripts","tmp",gsub("/main/0_mainPRE","",getwd())),"/settings/nsdm-settings.RData"))

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

# Retrieve reference rasters
rsts_ref<-readRDS(paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))

### =========================================================================
### C- species data
### =========================================================================
# Load species data
sp_dat<-readRDS(paste0(w_path,"tmp/",project,"/settings/species-occurences.rds"))
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

# Loop over species
for(ispi_name in species){

cat(paste0('Ready for GLO and REG modelling dataset preparation for ', ispi_name, '...\n'))

### =========================================================================
### D- Sample background absences
### =========================================================================
# D.1 GLO
sp_dat_glo<-sp_dat$spe_glo_dis[sp_dat$spe_glo_dis$species==ispi_name,]
pseu.abs_i_glo<-nsdm.pseuabs(n=n_back,
							 rst_ref=rsts_ref$rst_glo,
							 type=pa_po_glo,
                             pres=sp_dat_glo,
                             taxon=ispi_name)

cat(paste0('GLO dataset prepared (n_occ=',length(which(pseu.abs_i_glo@pa==1)),')...\n'))
                            
# D.2 REG
if(n_levels>1){
sp_dat_reg<-sp_dat$spe_reg_dis[sp_dat$spe_reg_dis$species==ispi_name,]
pseu.abs_i_reg<-nsdm.pseuabs(n=n_back,
							 rst_ref=rsts_ref$rst_reg,
							 type=pa_po_reg,
                             pres=sp_dat_reg,
                             taxon=ispi_name)
						   
cat(paste0('REG dataset prepared (n_occ=',length(which(pseu.abs_i_reg@pa==1)),')...\n'))

if(length(which(pseu.abs_i_reg@pa==1))==0) next}
if(length(which(pseu.abs_i_glo@pa==1))==0) next

# D.3 Save spatial points as shapefiles
dir<-paste0(scr_path,"/outputs/",project,"/d0_datasets/shp/",ispi_name)
dir.create(dir, recursive=T)
if(n_levels<2) rsts_ref$rst_reg<-rsts_ref$rst_glo				   
suppressWarnings(glo_pts<-spTransform(SpatialPointsDataFrame(coords=pseu.abs_i_glo@xy[which(pseu.abs_i_glo@pa==1),],
                               data=data.frame(pa=rep(1,length(which(pseu.abs_i_glo@pa==1)))),
                               proj4string=crs(rsts_ref$rst_glo)), crs(rsts_ref$rst_reg)))
writeOGR(glo_pts, dir, paste0(ispi_name, "_glo"), driver = "ESRI Shapefile", overwrite=T)

if(n_levels>1){
suppressWarnings(reg_pts<-SpatialPointsDataFrame(coords=pseu.abs_i_reg@xy[which(pseu.abs_i_reg@pa==1),],
                              data=data.frame(pa=rep(1,length(which(pseu.abs_i_reg@pa==1)))),
                               proj4string=crs(rsts_ref$rst_reg)))
writeOGR(reg_pts, dir, paste0(ispi_name, "_reg"), driver = "ESRI Shapefile", overwrite=T)}

### =========================================================================
### E- Save
### =========================================================================
if(n_levels>1){
l<-list(group=unique(sp_dat_reg$class),
        pseu.abs_i_glo=pseu.abs_i_glo,
        pseu.abs_i_reg=pseu.abs_i_reg)
} else {
l<-list(group=unique(sp_dat_glo$class),
        pseu.abs_i_glo=pseu.abs_i_glo)
}
		  						  
nsdm.savethis(object=l,
              species_name=ispi_name,
              save_path=paste0(scr_path,"/outputs/",project,"/d0_datasets/base"))
}

print("Finished!")