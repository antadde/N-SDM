#' nsdm.predict
#'
#' Compute predictions (parallel) for each target modelling algorithm
#'
#' @param models A nsdm.flex object containing fitted models
#' @param nwdata A data.frame object containing covariate data for predictions
#' @param nsplits Number of splits for parallel operations
#'
#' @return A vector of predicted values
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.predict <- function(models, nwdata, nsplits=NULL){

ndata_bck_f<-list() # list where results will be stored

# Split indices for parallel predict
splits_ix<-parallel::splitIndices(nrow(nwdata), nsplits)

for(m in 1:length(models@fits)){ # Loop over models (1 in regular cases but much more in esm settings)
  
# Retrieve model fit
model<-models@fits[[m]][[1]]

# Reorder covariates to match model order (needed for gbm booster ..)
cov<-unlist(strsplit(models@meta$env_vars,", "))
nwdata <- nwdata[,cov]

# In case of esm model uses its name its index for naming
if(model_name=="esm"){
model_nameu<-names(models@fits)[m]
}else{
  model_nameu<-model_name}

## GLM
if(c("glm") %in% class(model) & !("gam") %in% class(model)){
ndata<-parallel::mclapply(splits_ix, function(x){round(data.frame(fit=predict(model, nwdata[x,], se.fit = FALSE, type = "response")),2)}, mc.cores = nsplits)
ndata_bck<-rbind.fill(ndata)
}

## GAM
if(c("gam") %in% class(model)){
ndata<-parallel::mclapply(splits_ix, function(x){data.frame(fit=round(predict(model, nwdata[x,], se.fit = FALSE, gc.level=1, type = "response", block.size=round(nrow(nwdata[x,])/3)),2))}, mc.cores = nsplits)
ndata_bck<-rbind.fill(ndata)
}

## Maxent
if(c("maxnet") %in% class(model)){
ndata<-parallel::mclapply(splits_ix, function(x){data.frame(fit=round(predict(model, nwdata[x,], type="cloglog", clamp=F),2))}, mc.cores = nsplits)
ndata_bck<-rbind.fill(ndata)
}

## RF
if(c("randomForest") %in% class(model)){
ndata<-parallel::mclapply(splits_ix, function(x){data.frame(fit=predict(model, nwdata[x,], type="prob")[,2])}, mc.cores = nsplits)
ndata_bck<-rbind.fill(ndata)
}

## RF with ranger
if ("ranger" %in% class(model)) {
  # ndata <- parallel::mclapply(
    # splits_ix,
    # function(x) {
      # preds <- predict(model, data = nwdata[x, ], type = "response", verbose = FALSE)$predictions
      # if (is.matrix(preds)) {
        # data.frame(fit = preds[, 2])  # probability of class 1
      # } else {
        # data.frame(fit = preds)  # if only 2-class case returns vector
      # }
    # },
    # mc.cores = nsplits
  # )
  # ndata_bck <- plyr::rbind.fill(ndata)
future::plan("multicore", workers = nsplits)  # or "multicore" on Linux

ndata <- future.apply::future_lapply(
  splits_ix,
  function(x) {
    preds <- predict(model, data = nwdata[x, ], type = "response")$predictions
    if (is.matrix(preds)) {
      data.frame(fit = preds[, 2])
    } else {
      data.frame(fit = preds)
    }
  }
)
ndata_bck <- data.table::rbindlist(ndata, fill = TRUE)
}

## GBM
if(c("lgb.Booster") %in% class(model)){
ndata_bck<-data.frame(fit=predict(model, as.matrix(nwdata), num_threads=nsplits)) # already makes use of all available cores
}

ndata_bck_f[[m]]<-ndata_bck
names(ndata_bck_f)[m]<-model_nameu

}
if(model_name=="esm")return(ndata_bck_f)
if(model_name!="esm")return(ndata_bck_f[[1]])
}
