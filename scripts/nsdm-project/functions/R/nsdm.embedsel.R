#' nsdm.embedsel
#'
#' Embedded covariate selection
#'
#' @param covdata A data.frame containing covariate data
#' @param pa A vector of presences (1) and absences (0)
#' @param weights Vector of weights of length pa
#' @param force A character vector indicating the name(s) of covariates to be forced in final set (e.g.: mainEU)
#' @param nthreads Number (numeric) of cores to be used during parallel operations
#'
#' @return A ranked list of covariates selected after the embedded procedure
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.embedsel <- function(covdata, pa, weights, force=NULL, nthreads=detectCores()/2){
###
# GLM with elastic-net
###
# Do embeddded covariate selection
form<-as.formula(paste0("as.factor(pa) ~ " ,paste(paste0("poly(",names(covdata),",2)-1"),collapse=" + ")))
x <- model.matrix(form, covdata)
mdl.glm <- cv.glmnet(x, as.factor(pa), alpha=0.5, weights=weights, family = "binomial", type.measure = "deviance", parallel = TRUE)
# Extract results
glm.beta<-as.data.frame(as.matrix(coef(mdl.glm, s=mdl.glm$lambda.1se)))
if(plyr::empty(glm.beta)) glm.beta<-as.data.frame(as.matrix(coef(mdl.glm, s=mdl.glm$lambda.min)))
glm.beta<-data.frame(var = row.names(glm.beta)[which(glm.beta != 0)], coef=abs(glm.beta$s1))[-1,]
glm.beta<-data.frame(glm.beta[order(glm.beta$coef, decreasing = TRUE),], model="glm")
glm.beta$var<-stri_sub(glm.beta$var,6,-6)
glm.beta<-data.frame(setDT(glm.beta)[, .SD[which.max(coef)], by=var])
glm.beta$rank<-1:nrow(glm.beta)

print("1/3")

###
# GAM with null-space penalization
###
# Do embeddded covariate selection
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata),",bs='cr')"),collapse=" + ")))
mdl.gam <- mgcv::bam(form, data=cbind(covdata, as.factor(pa)), weights=weights, family="binomial", method="fREML", select=TRUE, discrete=TRUE, control=list(nthreads=nthreads))
t<-try(summary(mdl.gam), TRUE)
if(class(t)=="try-error"){
form<-as.formula(paste0("pa ~ " ,paste(paste0("s(",names(covdata),",bs='ts')"),collapse=" + ")))
mdl.gam <- mgcv::bam(form, data=cbind(covdata, as.factor(pa)), weights=weights, family="binomial", method="fREML", select=TRUE, discrete=TRUE, control=list(nthreads=nthreads))
}

# Extract results
gam.beta<-data.frame(var=names(mdl.gam$model)[-c(1,length(names(mdl.gam$model)))], summary(mdl.gam)$s.table, row.names = NULL)
gam.beta<-gam.beta[gam.beta$p.value<0.9,]
gam.beta<-data.frame(gam.beta[order(abs(gam.beta$Chi.sq), decreasing = TRUE),], rank = 1:nrow(gam.beta), model="gam")

print("2/3")

###
# Guided regularized RF
###
# Do embeddded covariate selection
rf <- RRF(covdata,as.factor(pa), flagReg = 0)
impRF <- rf$importance[,"MeanDecreaseGini"]
imp <- impRF/(max(impRF))#normalize the importance score
gamma <- 0.5
coefReg <- (1-gamma)+gamma*imp #weighted average
mdl.rf <- RRF(covdata,as.factor(pa), classwt=c("0"=min(weights), "1"=max(weights)), coefReg=coefReg, flagReg=1)
# Extract results
rf.beta<-data.frame(var = row.names(mdl.rf$importance), mdl.rf$importance,  row.names=NULL)
rf.beta<-rf.beta[which(rf.beta$MeanDecreaseGini > 0),]
rf.beta<-data.frame(rf.beta[order(rf.beta$MeanDecreaseGini, decreasing = TRUE),], rank = 1:nrow(rf.beta), model="rf")

print("3/3")

###
# Combine and return
###
res<-rbind(glm.beta[,c("var","rank","model")], gam.beta[,c("var","rank","model")], rf.beta[,c("var","rank","model")])

# If covariate to force not included, add it
if(is.character(force)){
if(length(res[which(res$var==force),"model"])<3){
models<-c("glm","rf","gam")
missing<-setdiff(models, res[which(res$var==force),"model"])
res<-rbind(data.frame(res), data.frame(var=force, rank=1, model=missing))
}}

# Return results
return(res)
}