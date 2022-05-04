#############################################################################
## 0_mainPRE
## Update routine
## Date: 25-09-2021
## Author: Antoine Adde 
#############################################################################
project<-gsub("/main","",gsub(".*scripts/","",getwd()))

# Load nsdm.ine settings
load(paste0(gsub("scripts","outputs",gsub("/main","",getwd())),"/settings/nsdm-settings.RData"))

# Load species list
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/species-list.rds"))

# Split species into cluster run chunks to match n_mx_spe
library(parallel)
splits<-splitIndices(length(species), length(species)/n_mx_spe)
run_id<-read.table(paste0(w_path,"outputs/",project,"/settings/tmp/run_id.txt"), )$V1
species<-species[splits[[run_id]]]
writeLines(as.character(length(species)), paste0(w_path,"outputs/",project,"/settings/tmp/n_spe.txt"))

# Save species list for specific run
saveRDS(species, paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))

# Update N-SDM settings
save.image(paste0(w_path,"outputs/",project,"/settings/nsdm-settings.RData"))