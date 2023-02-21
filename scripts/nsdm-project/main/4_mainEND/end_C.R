#############################################################################
## 4_mainEND
## Pre-fill ODMAP protocol
## Date: 21-02-2023
## Author: Antoine Adde 
#############################################################################
### =========================================================================
### A- Preparation
### =========================================================================
project<-gsub("/main/4_mainEND","",gsub(".*scripts/","",getwd()))

# Load nsdm settings
load(paste0(gsub("scripts","tmp",gsub("/main/4_mainEND","",getwd())),"/settings/nsdm-settings.RData"))

# Set permissions for new files
Sys.umask(mode="000")

# Set your working directory
setwd(w_path)

# Set lib path
.libPaths(lib_path)

# Load nsdm package
require(nsdm)

# Load ODMAP template
ODMAP<-read.csv(paste0(w_path, "/data/",project,"/odmap/ODMAP.csv"), sep=";")
ODMAP_short<-ODMAP[,c(1,5)]

# Load template raster
r<-readRDS(paste0(w_path,"tmp/",project,"/settings/ref-rasters.rds"))[[1]]

# Load covinfo table
lr<-readRDS(paste0(w_path,"tmp/",project,"/settings/covariates-list.rds"))
cov_info<-lr$cov_info

# species list
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

# Pre-fill with available information
ODMAP_short[9,2]<-as.character(extent(r)) # Spatial Extent (Lon / Lat)	9
ODMAP_short[10,2]<-res(r)[1]# Spatial resolution	10
ODMAP_short[16,2]<-paste(unique(cov_info$category), collapse=", ") # Climatic, topographic, edaphic, habitat, etc.	16
ODMAP_short[19,2]<-paste(unique(mod_algo), collapse=", ") # Model algorithms	19
ODMAP_short[21,2]<-paste0("weighting=", do_weighting) # Is model averaging/ensemble modelling used?	21
ODMAP_short[26,2]<-paste(unique(species), collapse=", ") # Taxon names	26
if(disag_reg==TRUE | disag_glo==TRUE){
ODMAP_short[33,2]<-"spatial and/or temporal disaggregation"  # Details on scaling, if applicable: e.g., rasterisation of polygon maps, spatial and temporal thinning, measures to address spatial uncertainties	33
} else {
ODMAP_short[33,2]<-"no spatial and/or temporal disaggregation"  # Details on scaling, if applicable: e.g., rasterisation of polygon maps, spatial and temporal thinning, measures to address spatial uncertainties	33
}
ODMAP_short[36,2]<-paste(n_back, back_strat, "background points", sep=" ") # Details on background data derivation, if applicable: e.g., spatial and temporal extent, spatial and temporal buffer, bias correction (e.g. target group sampling)	36
ODMAP_short[38,2]<-replicate_type # Selection of validation data (withheld from model fitting, used for estimating prediction error for model selection, model averaging or ensemble): e.g., cross-validation method	38
ODMAP_short[44,2]<-proj4string(r) # Map projection (coordinate reference system)	44
ODMAP_short[50,2]<-as.character(extent(r)) # Spatial extent	50
ODMAP_short[51,2]<-res(r)[1] # Spatial resolution	51
if(length(proj_periods) != 0){
ODMAP_short[52,2]<-paste(c("current", proj_periods), collapse=", ") # temporal extent/time period	52
} else {
ODMAP_short[52,2]<-"current" # temporal extent/time period	52
}
if(length(proj_scenarios) != 0) ODMAP_short[54,2]<-paste(unique(proj_scenarios), collapse=", ")# Models and scenarios used	54
if(do_weighting==TRUE){
ODMAP_short[65,2]<-paste(weight_metric, "metric derived weights", sep=" ")# Method for model averaging: e.g. derivation of weights 	65
} else {
ODMAP_short[65,2]<-"no weighting done"# Method for model averaging: e.g. derivation of weights 	65
}
if(do_weighting==TRUE){
ODMAP_short[66,2]<-"weighted average" # Ensemble method	66
} else {
ODMAP_short[66,2]<-"unweighted average" # Ensemble method	66
}
ODMAP_short[71,2]<-"performance statistics on training data: AUC, AUC_S, RMSE, CBI, maxSensitivity, maxSpecificity, maxAccuracy, maxPPV, maxNPV, maxJaccard, maxTSS, maxKappa, maxSEDI" #Performance statistics estimated on training data	71
ODMAP_short[72,2]<-"performance statistics on validation data: same" #Performance  statistics estimated on validation data (from data partitioning)	72

# Update and save ODMAP
ODMAP[,5]<-ODMAP_short[,2]
write.table(ODMAP, paste0(w_path,"/tmp/",project,"/ODMAP/ODMAP.csv"), quote=F, sep=";", row.names = FALSE)
