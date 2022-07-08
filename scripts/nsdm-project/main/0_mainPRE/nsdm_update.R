#############################################################################
## 0_mainPRE
## Update N-SDM routine
## Date: 20-05-2022
## Author: Antoine Adde 
#############################################################################
project<-gsub("/main","",gsub(".*scripts/","",getwd()))

# Load N-SDM settings
load(paste0(gsub("scripts","tmp",gsub("/main","",getwd())),"/settings/nsdm-settings.RData"))

# Load species list
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/species-list.rds"))

# Split species into cluster run chunks to match n_mx_spe
library(parallel)
splits<-splitIndices(length(species), length(species)/n_mx_spe)
run_id<-read.table(paste0(w_path,"tmp/",project,"/settings/tmp/run_id.txt"), )$V1
species<-species[splits[[run_id]]]
writeLines(as.character(length(species)), paste0(w_path,"tmp/",project,"/settings/tmp/n_spe.txt"))

# Save species list for specific run
saveRDS(species, paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))

# Update N-SDM settings
save.image(paste0(w_path,"tmp/",project,"/settings/nsdm-settings.RData"))