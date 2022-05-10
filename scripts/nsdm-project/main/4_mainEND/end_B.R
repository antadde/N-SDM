#############################################################################
## Pipeline sacct output analysis
## Date: 25-09-2021
## Author: Antoine Adde 
#############################################################################
project<-gsub("/main/4_mainEND","",gsub(".*scripts/","",getwd()))

# Load nsdm.ine settings
load(paste0(gsub("scripts","outputs",gsub("/main/4_mainEND","",getwd())),"/settings/nsdm-settings.RData"))

# Set permissions for new files
Sys.umask(mode="000")

# Set your working directory
setwd(w_path)

# Set lib path
.libPaths(lib_path)

# Load nsdm package
require(nsdm)

# Target species
species<-readRDS(paste0(w_path,"outputs/",project,"/settings/tmp/species-list-run.rds"))

# Retrieve and preprocess sacct outputs
sacct<-read.csv2(list.files(paste0("outputs/",project,"/sacct/"), pattern=paste0(ssl_id,"_",run_id,"_sacct"), recursive=T, full.names=T))
sacct<-sacct[-grep("extern", unlist(sacct), fixed=T),]
sacct<-read.table(text=sacct, skip = 1,
          fill=TRUE, col.names=c("JobID", "JobName", "Elapsed", "NCPUs", "TotalCPU", "CPUTime", "ReqMem", "MaxRSS", "MaxDiskRead", "MaxDiskWrite", "State", "ExitCode"))
ix2rep<-grep(names(sort(table(sacct$JobName), decreasing=T)[1]), sacct$JobName)
ix4rep<-ix2rep-1
sacct$JobName[ix2rep]<-sacct$JobName[ix4rep]
sacct$JobName<-toupper(sacct$JobName)
sacct<-sacct[-ix4rep,]
jobs<-c(paste0("PRE_",c("A","B")), paste0("GLO_",c("A", "B", "C")), paste0("LOC_", c("A", "B", "C")), paste0("FUT_", c("A", "B", "C", "D")))

sacct<-sacct[-which(!sacct$JobName %in% jobs),]
sacct$species<-NA

sacct<-sacct[head(which(match(sacct$JobName, "GLO_A")==1),1):nrow(sacct),]

j_na<-names(table(sacct$JobName))
v_na<-table(sacct$JobName)

n_spe<-length(species)

for(i in 1:length(j_na)){
j_na_i<-j_na[i]
if(round(v_na[i]/n_spe) == v_na[i]/n_spe){
if(grepl("FUT_A", j_na_i) | grepl("FUT_C", j_na_i)){
sacct$species[which(match(sacct$JobName, j_na_i)==1)]<-rep(species, each=length(mod_algo), times=length(periods)*length(scenarios))
}else if(grepl("FUT_B", j_na_i) | grepl("FUT_D", j_na_i)){
sacct$species[which(match(sacct$JobName, j_na_i)==1)]<-rep(species, each=1, times=length(periods)*length(scenarios))
}else{
sacct$species[which(match(sacct$JobName, j_na_i)==1)]<-rep(species, each=v_na[i]/n_spe)
}
} else {
sacct$species[which(match(sacct$JobName, j_na_i)==1)]<-"ALL"
}
}

spth<-paste0(scr_path,"/outputs/",project,"/d18_sacct/",paste0(ssl_id,"_run_",run_id))
dir.create(spth, recursive=T)
fwrite(sacct, paste0(spth,"/sacct_table_",ssl_id,"_run_",run_id,".csv"))

sacct$JobName<- factor(sacct$JobName, levels = jobs)
sacct$Elapsed<-60 * hours(times(sacct$Elapsed)) + minutes(times(sacct$Elapsed))

# Synthetize sacct info
n_spe<-length(species)
times_sm<-aggregate(times(sacct$Elapsed), list(sacct$JobName, sacct$species), sum)
times_mx<-aggregate(times(sacct$Elapsed), list(sacct$JobName, sacct$species), max)
times_sm$t_sum<-times_sm$x
times_mx$t_mx<-times_mx$x

convb <- function(x){
  ptn <- "(\\d*(.\\d+)*)(.*)"
  num  <- as.numeric(sub(ptn, "\\1", x))
  unit <- sub(ptn, "\\3", x)
  unit[unit==""] <- "1"

  mult <- c("1"=1, "K"=1000, "M"=1000^2, "G"=1000^3)
  num * unname(mult[unit])
}

mem_mn<-aggregate(convb(sacct$MaxRSS), list(sacct$JobName, sacct$species), mean)
mem_mx<-aggregate(convb(sacct$MaxRSS), list(sacct$JobName, sacct$species), max)
mem_mn$m_mn<-mem_mn$x/1000^3
mem_mx$m_mx<-mem_mx$x/1000^3

nsdm.savethis(object=list(mem_mn, mem_mx, times_sm, times_mx), species=paste0(ssl_id,"_run_",run_id),
              save_path=paste0(scr_path,"/outputs/",project,"/d18_sacct"))

# Plots
pals<-c(
get_palette(palette = "Greys", length(grep("PRE", unique(sacct$JobName)))),
get_palette(palette = "Blues", length(grep("GLO", unique(sacct$JobName)))),
get_palette(palette = "Greens", length(grep("LOC", unique(sacct$JobName)))),
get_palette(palette = "Reds", length(grep("FUT", unique(sacct$JobName)))))

spth<-paste0(scr_path,"/outputs/",project,"/plots/sacct/",ssl_id,"_run_",run_id)
dir.create(spth, recursive=T)
wd=10
ht=5

pdf(file = paste0(spth,"/m_mn.pdf"),   
    width = wd, 
    height = ht) 
print(ggboxplot(mem_mn, x="Group.1", y="m_mn",  fill="Group.1", ylab="max RSS_mn (gb)", xlab="step", title = "Max RSS (mean)", palette = pals, legend="none"))
dev.off()

pdf(file = paste0(spth,"/m_mx.pdf"),   
    width = wd, 
    height = ht)
print(ggboxplot(mem_mx, x="Group.1", y="m_mx", fill="Group.1", ylab="max RSS_mx (gb)", xlab="step", title = "Max RSS (max)", palette = pals, legend="none"))
dev.off()

pdf(file = paste0(spth,"/t_mx.pdf"),  
    width = wd, 
    height = ht) 
print(ggboxplot(times_mx, x="Group.1", y="t_mx",  fill="Group.1", ylab="t_max (min)", xlab="step", title = "Time (max)", palette = pals, legend="none"))
dev.off()

pdf(file = paste0(spth,"/t_sum.pdf"),
    width = wd,
    height = ht) 
print(ggboxplot(times_sm, x="Group.1", y="t_sum",  fill="Group.1", ylab="t_sum (min)", xlab="step", title = "Time (sum)", palette = pals, legend="none"))
dev.off()
