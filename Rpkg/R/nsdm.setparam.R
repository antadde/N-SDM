#' nsdm.setparam
#'
#' Define model parameter settings for various modeling algorithms.
#'
#' @param model_name A character string specifying the modeling algorithm to use (e.g., "glm", "gam", "rf", "gbm", "max", "esm").
#' @param param_grid A character string indicating the complete file path where the `.xlsx` document containing the parameter settings is stored.
#' @param covariate_names A character vector specifying the names of covariates to be included in the models.
#' @param tmp_path A character string specifying the directory where temporary outputs should be stored.
#' @param ncov.esm (Optional) A numeric value specifying the total number of covariates evaluated when using the ESM (Ensemble of Small Models) framework.
#' @param comb.esm (Optional) A numeric value indicating the number of covariates combined in each individual small model under the ESM framework.
#' @param weights A numeric vector specifying the weights to be applied in the models. If a single value is provided, no weighting is applied.
#' @param nthreads A numeric value specifying the number of cores to be used for parallel computations.
#'
#' @return A list containing elements of class `'multi.input'`, each representing a model configuration to be fitted.
#' 
#' @author Antoine Adde (\email{antoine.adde@eawag.ch})
#' @export

nsdm.setparam <- function(model_name, param_grid, covariate_names, tmp_path, ncov.esm=NULL, comb.esm=NULL, weights=1, nthreads){
# Load data
grid_list <- fread(param_grid)  # Ensure 'param_grid' is correctly assigned

# Initialize a list to store subdataframes
subdataframes <- list()

# Get unique algorithms
unique_algos <- unique(grid_list$Algorithm)

# Loop over each unique algorithm
for (algo in unique_algos) {
  # Subset data for the current algorithm
  df <- grid_list[Algorithm == algo]
  
  # Store original parameter names
  param_names <- df$Parameter  
  
  # Remove Algorithm and Parameter columns before transposing
  df_t <- df[, -c(1,2)]  

  # Transpose the dataframe
  df_t <- t(df_t)

  # Convert to a data frame
  df_t <- as.data.frame(df_t, stringsAsFactors = FALSE)

  # Assign correct column names
  colnames(df_t) <- param_names

  # Reset row names
  rownames(df_t) <- NULL  

  # Store in list with algorithm name as key
  subdataframes[[algo]] <- df_t
}

subdataframes <- lapply(subdataframes, function(df) {
  df[] <- lapply(df, function(col) {
    # Remove leading/trailing spaces
    col <- trimws(col)
    
    # Check if column is character and all non-NA values are numeric-like
    if (is.character(col) && all(grepl("^[-]?[0-9]*[.,]?[0-9]+$", col[!is.na(col)]))) {
      as.numeric(gsub(",", ".", col))  # Convert numbers with comma decimals to proper numeric
    } else {
      col  # Keep as character
    }
  })
  return(df)
})

# Weighting?
if(length(weights)==1){
  weighting=FALSE
} else {
  weighting=TRUE}
  
# Set model settings
modinp<-list() # list where all model settings will be stored

## GLM
if(model_name=="glm"){
params_glm<-na.omit(data.frame(subdataframes[["glm"]]))

multis_glm<-list()
for(p in 1:nrow(params_glm)){
param_glm<-params_glm[p,]
form_glm<-as.formula(paste("Presence~", paste(paste0("poly(", covariate_names,",", as.numeric(param_glm), ", raw=FALSE)"), collapse="+")))
multi_glm<-nsdm.multi("glm", list(formula=form_glm, family="binomial"), tag=paste0("glm-",p), weight=weighting)
modinp<-append(modinp, multi_glm)
}
}

## GAM
if(model_name=="gam"){
mgcv_gam<<-mgcv::bam
mgcv_gam_fx<<-mgcv::bam

## with penalization
if(nrow(subdataframes[["gam.auto"]])>1){
  params_gam_auto<-na.omit(expand.grid(data.frame(subdataframes[["gam.auto"]])))
} else {
params_gam_auto<-na.omit(data.frame(subdataframes[["gam.auto"]]))
}

multis_gam_auto<-list()
for(p in 1:nrow(params_gam_auto)){
param_gam<-params_gam_auto[p,]
form_gam<-as.formula(paste0("Presence~ " , paste(paste0("s(", covariate_names, ",bs='", param_gam$reg.spline,"')"), collapse=" + ")))
multi_gam<-nsdm.multi("mgcv_gam", list(formula=form_gam, family="binomial", method=as.character(param_gam$method), select=FALSE, control=list(nthreads=nthreads)), tag=paste0("gam-",p), weight=weighting)
multis_gam_auto<-append(multis_gam_auto, multi_gam)
}

## pure regression splines without penalization
params_gam_fx<-na.omit(data.frame(subdataframes[["gam.fx"]]))

multis_gam_fx<-list()
for(p in 1:nrow(params_gam_fx)){
param_gam<-params_gam_fx[p,]
form_gam<-as.formula(paste0("Presence~ " ,paste(paste0("s(",covariate_names,", k=", as.numeric(param_gam),", fx=TRUE)"),collapse=" + ")))
multi_gam<-nsdm.multi("mgcv_gam", list(formula=form_gam, family="binomial", method="REML", select=FALSE),tag=paste0("gam-",p+nrow(params_gam_auto)), weight=weighting)
multis_gam_fx<-append(multis_gam_fx, multi_gam)
}
modinp<-append(multis_gam_auto,multis_gam_fx)
}

## Maxent; Note for Maxent default weights: "upweighting of background points" (https://par.nsf.gov/servlets/purl/10079053)
if(model_name=="max"){
if(nrow(subdataframes[["max"]])>1){
  params_maxent<-na.omit(expand.grid(data.frame(subdataframes[["max"]])))
} else {
  params_maxent<-na.omit(data.frame(subdataframes[["max"]]))
}

for(p in 1:nrow(params_maxent)){
param_maxent<-params_maxent[p,]
multi_maxent<-nsdm.multi("maxnet",list(classes = as.character(param_maxent$classes), addsamplestobackground=F, regmult = as.numeric(param_maxent$regmult)), tag=paste0("max-",p))
modinp<-append(modinp, multi_maxent)
}
}

# ## RF
if(model_name=="rf"){
if(nrow(subdataframes[["rf"]])>1){
  params_rf<-na.omit(expand.grid(data.frame(subdataframes[["rf"]])))
} else {
  params_rf<-na.omit(data.frame(subdataframes[["rf"]]))
}

for(p in 1:nrow(params_rf)){
param_rf<-params_rf[p,]
form_rf<-as.formula(paste("Presence~", paste(covariate_names, collapse=" + ")))
multi_rf<-nsdm.multi("randomForest",list(formula=form_rf,ntree=as.numeric(param_rf$num.trees), mtry=floor(sqrt(length(covariate_names))), nodesize = as.numeric(param_rf$min.node.size), classwt=c("0"=min(weights), "1"=max(weights))), tag=paste0("rf-",p), weight=weighting)
modinp<-append(modinp, multi_rf)
}
}

## (light)GBM
if(model_name=="gbm"){
if(nrow(subdataframes[["gbm"]])>1){
  params_gbm<-na.omit(expand.grid(data.frame(subdataframes[["gbm"]])))
} else {
  params_gbm<-na.omit(data.frame(subdataframes[["gbm"]]))
}
tmp_path_gbm<-paste0(tmp_path, "/gbm")
dir.create(tmp_path_gbm, recursive=TRUE)

for(p in 1:nrow(params_gbm)){
param_gbm<-params_gbm[p,]		   
multi_gbm <- nsdm.multi(
  "lgb.train",
  list(
    params = list(
      num_leaves = 2^as.numeric(param_gbm$max_depth) - 1,
      min_data_in_leaf = as.numeric(param_gbm$min_data_in_leaf),
      max_depth = as.numeric(param_gbm$max_depth),
      objective = "binary",
      learning_rate = as.numeric(param_gbm$learning_rate)
    ),
    nrounds = as.numeric(param_gbm$num_iterations),
    verbose = -10
  ),
  weight = weighting,
  tag = paste0("gbm-", p)
)

modinp<-append(modinp, multi_gbm)
}
}

## ESM
if(model_name=="esm"){
covs<-covariate_names[1:ncov.esm]
combs<-combn(covs, comb.esm)

### Generate list of all possible formulas
formulas<-list()
for (j in 1:ncol(combs)){
  formula <- as.formula(paste0("Presence ~", paste0("poly(", combs[,j],",2)", collapse = "+")))
  formulas <- c(formulas, formula)
}

### Set esm settings
for(p in 1:length(formulas)){
  form_glm<-formulas[[p]]
  multi_esm<-nsdm.multi("glm", list(formula=form_glm, family="binomial"), tag=paste0("esm-",p), weight=weighting)
  modinp<-append(modinp, multi_esm)
}
}
names(modinp)<- unlist(lapply(modinp, function(x){x@tag}))
return(modinp)}