#' covsel.embed
#'
#' Covariate selection with model-specific embedding (Step-2)
#'
#' @param covdata data.frame containing covariate data (continuous values) extracted at 'pa' locations
#' @param pa numeric vector of species presences (1) and absences (0)
#' @param weights numeric vector containing the weights for each value in 'pa' (of length 'pa')
#' @param force optional character vector indicating the name(s) of the covariate(s) to be forced in the final set
#' @param algorithms character vector indicating the name(s) of the algorithms(s) to be used for the embedding procedure (options: 'glm', 'gam', 'rf')
#' @param ncov target number of covariates to include in the final set
#' @param maxncov maximum possible number of covariates in the final set
#' @param nthreads number of cores to be used during parallel operations
#'
#' @return A list with three objects: (i) a data.frame with the covariates selected after the regularization/penalization and ranking procedures (covdata), (ii) a data.frame with the individual ranks of all covariates for each target algorithm (ranks_1), (iii) a data.frame with the final average ranks of selected covariates (ranks_2)
#' @author Antoine Adde (antoine.adde@unil.ch)
#' @examples
#' library(covsel)
#' covdata<-data_covfilter
#' dim(covdata)
#' covdata_embed<-covsel.embed(covdata, pa=data_covsel$pa, algorithms=c('glm','gam','rf'))
#' dim(covdata_embed$covdata)
#' @export

covsel.embed <- function(covdata, pa, weights=NULL, force=NULL, algorithms=c('glm','gam','rf'), ncov=ceiling(log2(length(which(pa==1)))), maxncov=12, nthreads=detectCores()/2){

# pre-settings
ranks_1<-data.frame()
if(!is.numeric(weights)) weights<-rep(1, length(pa))

###
###
# Embedding
###
###
if('glm' %in% algorithms){
###
# GLM (elastic-net)
###
# embedded covariate selection
form<-as.formula(paste0("as.factor(pa) ~ " ,paste(paste0("poly(",names(covdata),",2)"),collapse=" + "),"-1"))
x <- model.matrix(form, covdata)
mdl.glm <- suppressWarnings(cv.glmnet(x, as.factor(pa), alpha=0.5, weights=weights, family = "binomial", type.measure = "deviance", parallel = TRUE))
# Extract results
glm.beta<-as.data.frame(as.matrix(coef(mdl.glm, s=mdl.glm$lambda.1se)))
glm.beta<-data.frame(covariate = row.names(glm.beta), coef=as.numeric(abs(glm.beta[,1])))[which(glm.beta != 0),][-1,]
if(nrow(glm.beta)<1){
glm.beta<-as.data.frame(as.matrix(coef(mdl.glm, s=mdl.glm$lambda.min)))
glm.beta<-data.frame(covariate = row.names(glm.beta), coef=as.numeric(abs(glm.beta[,1])))[which(glm.beta != 0),][-1,]
}
if(nrow(glm.beta)<1){
print("No covariate selected after elastic-net regularization, skipping to next algorithm")
} else {
glm.beta<-data.frame(glm.beta[order(glm.beta$coef, decreasing = TRUE),], model="glm")
glm.beta$covariate<-stri_sub(glm.beta$covariate,6,-6)
glm.beta<-data.frame(setDT(glm.beta)[, .SD[which.max(coef)], by=covariate])
glm.beta$rank<-1:nrow(glm.beta)
ranks_1<-rbind(ranks_1, glm.beta[,c("covariate","rank","model")])
}
}

if('gam' %in% algorithms){
###
# GAM (null-space penalization)
###
# put aside forced covariates with < 10 unique points (required for default mgcv settings)
if(is.character(force)){
pointless10<-integer(1); names(pointless10)<-"pointless10"
df_force<-data.frame(covdata[,force]); names(df_force)<-force
pointless10<-which(apply(df_force, 2, function(x) length(unique(x)))<10)
# embedded covariate selection
if(length(pointless10)>0){
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata)[names(covdata) != names(pointless10)],",bs='cr')"),collapse=" + ")))
} else {
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata),",bs='cr')"),collapse=" + ")))
}
} else {
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata),",bs='cr')"),collapse=" + ")))
}
mdl.gam <- suppressWarnings(mgcv::bam(form, data=cbind(covdata, as.factor(pa)), weights=weights, family="binomial", method="fREML", select=TRUE, discrete=TRUE, control=list(nthreads=nthreads)))
t<-try(summary(mdl.gam), TRUE)
if(class(t)=="try-error"){
if(is.character(force)){
if(length(pointless10)>0){
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata)[names(covdata) != names(pointless10)],",bs='ts')"),collapse=" + ")))
} else {
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata),",bs='ts')"),collapse=" + ")))
}
}
mdl.gam <- suppressWarnings(mgcv::bam(form, data=cbind(covdata, as.factor(pa)), weights=weights, family="binomial", method="fREML", select=TRUE, discrete=TRUE, control=list(nthreads=nthreads)))
}
# Extract results
gam.beta<-data.frame(covariate=names(mdl.gam$model)[! names(mdl.gam$model) %in% c('(weights)', 'pa')], summary(mdl.gam)$s.table, row.names = NULL)
gam.beta<-gam.beta[gam.beta$p.value<0.9,]
if(nrow(gam.beta)<1){print("No covariate selected after GAM (null-space penalization), skipping to next algorithm")
} else {
gam.beta<-data.frame(gam.beta[order(abs(gam.beta$Chi.sq), decreasing = TRUE),], rank = 1:nrow(gam.beta), model="gam")
ranks_1<-rbind(ranks_1, gam.beta[,c("covariate","rank","model")])
}
}

if('rf' %in% algorithms){
###
# RF (guided regularized random forest)
###
# embedded covariate selection
rf <- RRF(covdata,as.factor(pa), flagReg = 0)
impRF <- rf$importance[,"MeanDecreaseGini"]
imp <- impRF/(max(impRF))
gamma <- 0.5
coefReg <- (1-gamma)+gamma*imp
mdl.rf <- RRF(covdata, as.factor(pa), classwt=c("0"=min(weights), "1"=max(weights)), coefReg=coefReg, flagReg=1)
# Extract results
rf.beta<-data.frame(covariate = row.names(mdl.rf$importance), mdl.rf$importance,  row.names=NULL)
rf.beta<-rf.beta[which(rf.beta$MeanDecreaseGini > 0),]
if(nrow(rf.beta)<1){print("No covariate selected after RF (guided regularized random forest)")
} else {
rf.beta<-data.frame(rf.beta[order(rf.beta$MeanDecreaseGini, decreasing = TRUE),], rank = 1:nrow(rf.beta), model="rf")
ranks_1<-rbind(ranks_1, rf.beta[,c("covariate","rank","model")])
}
}

###
###
# Ranking
###
###
if(nrow(ranks_1)<1){print("No covariate selected after the embedding procedure ...")
return(NULL)
} else {
# Rank covariates selected commonly by the algorithms (intersect)
intersect.tmp<-ranks_1[ranks_1$covariate %in% names(which(table(ranks_1$covariate) == length(unique(ranks_1$model)))),]
intersect.tmp<-aggregate(intersect.tmp[,c("rank")], list(intersect.tmp$covariate), sum); colnames(intersect.tmp)<-c("covariate","rank")
intersect.sel<-data.frame(intersect.tmp[order(intersect.tmp$rank, decreasing = FALSE),], rank.f = 1:nrow(intersect.tmp))
# Rank and add other covariates (union), if needed
union.tmp<-ranks_1[ranks_1$covariate%in%names(which(table(ranks_1$covariate) < length(unique(ranks_1$model)))),]
if(nrow(union.tmp)>0){
union.tmp<-aggregate(union.tmp[,c("rank")], list(union.tmp$covariate), sum); colnames(union.tmp)<-c("covariate","rank")
union.sel.tmp<-data.frame(union.tmp[order(union.tmp$rank, decreasing = FALSE),], rank.f = (max(intersect.sel$rank.f+1)):(max(intersect.sel$rank.f)+nrow(union.tmp)))
ranks_2<-rbind(intersect.sel,union.sel.tmp)
} else {
ranks_2<-intersect.sel
}

###
###
# Sub-setting
###
###
# subset ranking
if(ncov>maxncov) ncov<-maxncov
if(ncov>nrow(ranks_2)) ncov<-nrow(ranks_2)
ranks_2<-ranks_2[1:ncov,]
# forced covariate(s)
if(is.character(force)){
tf<-force[which(!(force%in%ranks_2$covariate))]
if(length(tf>1)){
toforce<-data.frame(covariate=tf, rank="forced", rank.f="forced")
ranks_2[c(nrow(ranks_2)-nrow(toforce)+1):c(nrow(ranks_2)),]<-toforce
}
}
ranks_2<-ranks_2[,c("covariate", "rank.f")]
# subset covariate set
covdata<-covdata[sub('.*\\.', '', unlist(ranks_2["covariate"]))]

# Return results
return(list(covdata=covdata, ranks_1=ranks_1, ranks_2=ranks_2))
}
}
