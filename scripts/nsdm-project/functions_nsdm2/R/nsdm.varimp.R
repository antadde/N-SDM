#' nsdm.varimp
#'
#' Covariate importance computation
#'
#' @param models A nsdm.fit object containing fitted model(s)
#'
#' @return A list containing the values of the relative importance for each individual covariate and model
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.varimp <- function(models){

imp_scaled_f<-list() # list where results will be stored
  
for(m in 1:length(models@fits)){ # Loop over models (1 in regular cases but much more in esm settings)

# Retrieve model fit
model<-models@fits[[m]][[1]]	  
  
# In case of esm model uses its name its index for naming and subset Data
if(model_name=="esm"){
  model_nameu<-names(models@fits)[m]
}else{
  model_nameu<-model_name}

# Compute importance  
if(c("lgb.Booster") %in% class(model)){
imp<-data.frame(lgb.importance(model))
imp_scaled <- data.frame(Variable=imp$Feature, Importance=scale(imp$Gain,FALSE,max(imp$Gain)))
}

if(c("glm") %in% class(model)){
imp<-varImp(model, scale=FALSE)
imp$variable<-gsub("[\\(\\,]", "", regmatches(rownames(imp), gregexpr("\\(.*?\\,", rownames(imp))))
imp<-aggregate(imp[,"Overall"], list(imp$variable), mean)
imp_scaled <- data.frame(Variable=imp$Group.1, Importance=scale(imp$x,FALSE,max(imp$x)))
}

if(c("gam") %in% class(model)){
imp<-summary(model)$s.table[,"Chi.sq"]
imp_scaled <- data.frame(Variable=names(imp), Importance=scale(imp,FALSE,max(imp)))
imp_scaled$Variable<-gsub("[\\(\\)]", "", regmatches(imp_scaled$Variable, gregexpr("\\(.*?\\)", imp_scaled$Variable)))
}

if(c("randomForest") %in% class(model)){
imp<-randomForest::importance(model)
imp_scaled <- data.frame(Variable=rownames(imp), Importance=scale(imp,FALSE,max(imp)))
}

if(c("ranger") %in% class(model)){
  imp <- model$variable.importance
  imp_scaled <- data.frame(
    Variable = names(imp),
    Importance = scale(imp, center = FALSE, scale = max(imp))
  )
}

if(c("maxnet") %in% class(model)){
variables<-names(model$levels)
imp<-data.frame()
for(i in 1:length(variables)){
ix<-grep(variables[i],names(model$betas))
score<-sum(abs(model$betas[ix]))
imp<-rbind(imp, data.frame(score=score, variable=variables[i]))
}
imp_scaled <- data.frame(Variable=imp$variable, Importance=scale(imp$score,FALSE,max(imp$score)))
}

row.names(imp_scaled)<-NULL
names(imp_scaled)<-c("Covariate", "Importance")
imp_scaled<-imp_scaled[order(imp_scaled$Importance, decreasing=T),]

imp_scaled_f[[m]]<-imp_scaled
names(imp_scaled_f)[m]<-model_nameu
}
if(model_name=="esm")return(imp_scaled_f)
if(model_name!="esm")return(imp_scaled_f[[1]])}