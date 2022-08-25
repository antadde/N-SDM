#############################################################################
## 4_mainEND
## A: final evaluation ensembles of nested-ensembles
## Date: 20-05-2022
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

### =========================================================================
### B- Definitions
### =========================================================================
# SBATCH param
args<-eval(parse(text=args))
arrayID<-eval(parse(text=arrayID))

# SBATCH array
species<-readRDS(paste0(w_path,"tmp/",project,"/settings/tmp/species-list-run.rds"))
ispi_name <- species[arrayID]

# Scale-nesting methods for combining GLO and LOC predictions
nesting_methods<-nesting_methods

# Target model algorithms
models<-mod_algo

cat(paste0('Ready for evaluating ensemble predictions obtained for ', ispi_name, ' ...\n'))

### =========================================================================
### C- Ensemble predictions
### =========================================================================

### =========================================================================
### GLO
### =========================================================================
################# Retrieve predictions
level<-"glo"
## Loop on target algorithms
pred_all<-list()
for(model in models){
# Retrieve evaluation table and identify best model
eval_list<-nsdm.loadthis(model_name=model, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d3_evals/glo"))

if(model!="esm"){
  if(ncol(eval_list)>1){
    ord<-sort(eval_list[best_met,],decreasing=T)
    modinp_top<-names(ord[1])
  } else {
    modinp_top<-colnames(eval_list)}}
if(model=="esm"){
  ord<-sort(eval_list[best_met,],decreasing=T)
  ord<-ord[ord>best_thre_esm]
  modinp_top<-names(ord)
}

# Load model(s)
mod_m<-nsdm.loadthis(species_name=ispi_name, model_name=model,
               tag=paste0(model,"_tune"),
               read_path=paste0(scr_path,"/outputs/",project,"/d2_models/glo"))$model

# Loop over models (one in regular cases; more for ESMs)
pop<-list()
for(n in 1:length(modinp_top)){
modinp_top_n<-modinp_top[n]

# Retrieve test and training data 
testa_GLO<-lapply(mod_m$m@tesdat,function(x){
y<-x[,-which(colnames(x)=="Presence"),drop=FALSE]
})

 papa_GLO<-lapply(mod_m$m@tesdat,function(x){
 y<-x[,"Presence"]})
  
# Predict
outerloop<-length(mod_m$m@tesdat)
tmp_path_gbm<-paste0(scr_path,"/tmp/",project,"/gbm")
pred<-list()
for(k in 1:outerloop){
if("lgb.Booster" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
mod_m$m@fits[[modinp_top_n]][[k]]<-readRDS.lgb.Booster(paste0(tmp_path_gbm, "/", ispi_name,"_rep",k,"_mod",gsub(".*-","",modinp_top_n),"_",level,".rds"))
testa_GLO[[k]]<- testa_GLO[[k]][,-which(colnames(testa_GLO[[k]]) %in% c("X","Y"))]}
if("try-error" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
pred_i<-rep(NA, nrow(testa_GLO[[k]]))
} else {
pred_i<-nsdm.prd(mod_m$m@fits[[modinp_top_n]][[k]], testa_GLO[[k]])}
pred[[k]]<-pred_i
}
pop_n<-do.call(cbind, pred)
pop[[n]]<-pop_n
}

if(length(pop)>1){
pop<-simplify2array(pop)
pop<-apply(pop, c(1,2), mean)
}else{
pop<-pop[[1]]}

# List results
pred_all[[model]]<-pop
}

GLO_preds<-simplify2array(pred_all)

if(n_levels>1){
### =========================================================================
### LOC MULTIPLY
### =========================================================================
if("multiply" %in% nesting_methods){
################# Retrieve predictions
level<-"loc_multiply"
## Loop on target algorithms
pred_all<-list()
for(model in models){
# Retrieve evaluation table and identify best model
eval_list<-nsdm.loadthis(model_name=model, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d3_evals/loc/multiply"))

if(model!="esm"){
  if(ncol(eval_list)>1){
    ord<-sort(eval_list[best_met,],decreasing=T)
    modinp_top<-names(ord[1])
  } else {
    modinp_top<-colnames(eval_list)}}
if(model=="esm"){
  ord<-sort(eval_list[best_met,],decreasing=T)
  ord<-ord[ord>best_thre_esm]
  modinp_top<-names(ord)
}

# Load best model
mod_m<-nsdm.loadthis(species_name=ispi_name, model_name=model,
               tag=paste0(model,"_tune"),
               read_path=paste0(scr_path,"/outputs/",project,"/d2_models/loc/multiply"))$model

# Loop over models (one in regular cases; more for ESMs)
pop<-list()
for(n in 1:length(modinp_top)){
modinp_top_n<-modinp_top[n]

# Retrieve test and training data 
testa_LOC_multiply<-lapply(mod_m$m@tesdat,function(x){
y<-x[,-which(colnames(x)=="Presence"),drop=FALSE]
})

 papa_LOC_multiply<-lapply(mod_m$m@tesdat,function(x){
 y<-x[,"Presence"]})
  
# Predict
outerloop<-length(mod_m$m@tesdat)
tmp_path_gbm<-paste0(scr_path,"/tmp/",project,"/gbm")
pred<-list()
glo_prob<-list(); glo_out<-readRDS(list.files(paste0(scr_path,"outputs/",project,"/d8_ensembles/glo/", ispi_name), pattern=".rds", full.names=T))
for(k in 1:outerloop){
if("lgb.Booster" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
mod_m$m@fits[[modinp_top_n]][[k]]<-readRDS.lgb.Booster(paste0(tmp_path_gbm, "/", ispi_name,"_rep",k,"_mod",gsub(".*-","",modinp_top_n),"_",level,".rds"))
testa_LOC_multiply[[k]]<- testa_LOC_multiply[[k]][,-which(colnames(testa_LOC_multiply[[k]]) %in% c("X","Y"))]}
if("try-error" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
pred_i<-rep(NA, nrow(testa_LOC_multiply[[k]]))
} else {
pred_i<-nsdm.prd(mod_m$m@fits[[modinp_top_n]][[k]], testa_LOC_multiply[[k]])}
pred[[k]]<-pred_i
glo_prob[[k]]<-raster::extract(glo_out, data.frame(mod_m$m@tesdat[[k]]$X,mod_m$m@tesdat[[k]]$Y))/100
}
pop_n<-do.call(cbind, pred)
pop[[n]]<-pop_n
}

if(length(pop)>1){
pop<-simplify2array(pop)
pop<-apply(pop, c(1,2), mean)
}else{
pop<-pop[[1]]}

# List results
pred_all[[model]]<-pop
}

LOC_multiply_preds<-simplify2array(pred_all)
}

### =========================================================================
### LOC COVARIATE
### =========================================================================
if("covariate" %in% nesting_methods){
# Load covariate matrix to rescale GLO predictions
mat<-nsdm.loadthis(species_name=ispi_name, read_path=paste0(scr_path,"/outputs/",project,"/d1_covsels/loc"))$env_vars

################# Retrieve predictions
level<-paste0("loc_covariate")
## Loop on target algorithms
pred_all<-list()
for(model in models){
# Retrieve evaluation table and identify best model
eval_list<-nsdm.loadthis(model_name=model, species_name=ispi_name,
              read_path=paste0(scr_path,"/outputs/",project,"/d3_evals/loc/covariate"))

if(model!="esm"){
  if(ncol(eval_list)>1){
    ord<-sort(eval_list[best_met,],decreasing=T)
    modinp_top<-names(ord[1])
  } else {
    modinp_top<-colnames(eval_list)}}
if(model=="esm"){
  ord<-sort(eval_list[best_met,],decreasing=T)
  ord<-ord[ord>best_thre_esm]
  modinp_top<-names(ord)
}

# Load best model
mod_m<-nsdm.loadthis(species_name=ispi_name, model_name=model,
               tag=paste0(model,"_tune"),
               read_path=paste0(scr_path,"/outputs/",project,"/d2_models/loc/covariate"))$model

# Loop over models (one in regular cases; more for ESMs)
pop<-list()
for(n in 1:length(modinp_top)){
modinp_top_n<-modinp_top[n]

# Retrieve test and training data 
testa_LOC_covariate<-lapply(mod_m$m@tesdat,function(x){
y<-x[,-which(colnames(x)=="Presence"),drop=FALSE]
})

 papa_LOC_covariate<-lapply(mod_m$m@tesdat,function(x){
 y<-x[,"Presence"]})
   
# Predict
outerloop<-length(mod_m$m@tesdat)
tmp_path_gbm<-paste0(scr_path,"/tmp/",project,"/gbm")
pred<-list()
for(k in 1:outerloop){
if("lgb.Booster" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
mod_m$m@fits[[modinp_top_n]][[k]]<-readRDS.lgb.Booster(paste0(tmp_path_gbm, "/", ispi_name,"_rep",k,"_mod",gsub(".*-","",modinp_top_n),"_",level,".rds"))
testa_LOC_covariate[[k]]<- testa_LOC_covariate[[k]][,-which(colnames(testa_LOC_covariate[[k]]) %in% c("X","Y"))]}
if("try-error" %in% class(mod_m$m@fits[[modinp_top_n]][[k]])){
pred_i<-rep(NA, nrow(testa_LOC_covariate[[k]]))
} else {
pred_i<-nsdm.prd(mod_m$m@fits[[modinp_top_n]][[k]], testa_LOC_covariate[[k]])}
pred[[k]]<-pred_i
}
pop_n<-do.call(cbind, pred)
pop[[n]]<-pop_n
}

if(length(pop)>1){
pop<-simplify2array(pop)
pop<-apply(pop, c(1,2), mean)
}else{
pop<-pop[[1]]}

# List results
pred_all[[model]]<-pop
}
LOC_covariate_preds<-simplify2array(pred_all)
}
}

### =========================================================================
### D- Evaluate ensemble predictions
### =========================================================================
scores_ensemble<-list()
# GLO-level ensemble
target<-GLO_preds
scores<-list()
for (z in 1:outerloop){
score<-nsdm.ceval(f=rowMeans(as.data.frame(target[,z,]), na.rm=T),
                   pa=papa_GLO[[z]],
                   tesdat=testa_GLO[[z]],
                   crit=eval_crit)
scores[[z]]<-score}
scores_ensemble[["GLO"]]<-scores	

if(n_levels>1){
# LOC-level ensemble without nesting (only possible in multiple mode)
if("multiply" %in% nesting_methods){
target<-LOC_multiply_preds
scores<-list()
for (z in 1:outerloop){
score<-nsdm.ceval(f=rowMeans(as.data.frame(target[,z,]),na.rm=T),
                   pa=papa_LOC_multiply[[z]],
                   tesdat=testa_LOC_multiply[[z]],
                   crit=eval_crit)
scores[[z]]<-score}
scores_ensemble[["LOC"]]<-scores
}	

# Covariate nested ensemble
if("covariate" %in% nesting_methods){
target<-LOC_covariate_preds
scores<-list()
for (z in 1:outerloop){
score<-nsdm.ceval(f=rowMeans(as.data.frame(target[,z,]),na.rm=T),
                   pa=papa_LOC_covariate[[z]],
                   tesdat=testa_LOC_covariate[[z]],
                   crit=eval_crit)
scores[[z]]<-score}
scores_ensemble[["COV"]]<-scores
}

# Multiply nested ensemble
if("multiply" %in% nesting_methods){
target<-LOC_multiply_preds
glo_prob2<-lapply(glo_prob, function(eux){
eux_i<-which(eux<0)
eux[eux_i]<-0
return(eux)
})
scores<-list()
for (z in 1:outerloop){
score<-nsdm.ceval(f=sqrt(rowMeans(as.data.frame(target[,z,]), na.rm=T)*glo_prob2[[z]]),
                   pa=papa_LOC_multiply[[z]],
                   tesdat=testa_LOC_multiply[[z]],
                   crit=eval_crit)
scores[[z]]<-score}
scores_ensemble[["MUL"]]<-scores	
}
}

#### Save scores
scores_array<-lapply(scores_ensemble, simplify2array)
scores_array<-lapply(scores_array, rowMeans)
print(scores_array)

nsdm.savethis(object=scores_array,
              species_name=ispi_name,
			  compression=TRUE,
              save_path=paste0(scr_path,"/outputs/",project,"/d12_evals-final"))
			  
cat(paste0('Nested ensembles evaluated \n'))
cat(paste0('Finished!\n'))
