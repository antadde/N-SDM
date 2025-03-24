#' nsdm.respcurve
#'
#' Covariate response curves computation
#'
#' @param models A nsdm.fit object containing fitted model(s)
#' @param Data A data.frame object containing covariate values used for model fitting
#' @param scaleparam Attribute part of a scaled matrix containing the scaling parameters to be reapplied
#' @param factor When applicable, the value by which original covariates have been multiplied (e.g. by 10, 100, etc.) for storage purposes
#' @param model_name A character string indicating the name of the target modelling algorithm
#' @param species_name A character string indicating the name of the target taxon
#' @param plotting Logical. Plot response curves?
#' @param save_path A character string indicating the path where the plots will be saved
#' @param ncores Number of cores to be used during parallel operations
#'
#' @return A list containing the values of the response curves for each individual covariate and model
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.respcurve <- function(models, Data, scaleparam=NULL, factor=1, model_name, species_name, plotting=TRUE, save_path=NULL, ncores=1){

#########
## MAX ##
#########
if(class(models@fits[[1]][[1]])[1]=="maxnet"){
temp<-list() # list where results will be stored
NbVar <- ncol(Data)
mod<-models@fits[[1]][[1]]
type <- "cloglog"
mm <- mod$samplemeans
temp_m <- array(0, dim=c(100, 2, NbVar), dimnames=list(NULL, c("Var", "Pred"), colnames(Data)))
for(v in 1:NbVar){
max <- mod$varmax[v]
min <- mod$varmin[v]
levels <- unlist(mod$levels[v])
nr <- if (is.null(levels)) 100 else length(levels)
m <- data.frame(matrix(mm, nr, length(mm), byrow = T))
colnames(m) <- names(mm)
m[,v] <- if (!is.null(levels)) levels else 
      seq(min - 0.1*(max-min), max+0.1*(max-min), length=100)
preds <- suppressWarnings(predict(mod, m, type = type))
temp_m[,1,v]<-m[,v]
temp_m[,2,v]<-preds
}

# Rescale using scale parameters
	if(!is.null(scaleparam)){
    for(i in 1:NbVar){
      temp_m[,1,i]<-temp_m[,1,i]*scaleparam$"scaled:scale"[names(Data)[i]] + scaleparam$"scaled:center"[names(Data)[i]]
	}}
	
# Plot
	if(plotting==TRUE){
    save_this_path<-paste(save_path, species_name, model_name, sep="/")
    suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
    pdf(paste0(save_this_path,"/", paste(species_name,model_name,sep="_"),".pdf"))
  
    ## plotting window
    W.width <- ceiling(sqrt(NbVar))
    W.height <- ceiling(NbVar/W.width)
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
	  par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n", axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves:", species_name, toupper(model_name), sep=" "),col="#4c57eb")
    par(mar = c(3,3,1.5,1))
	
	## do plot
	for(i in 1:dim(temp_m)[3]){
	temp_i<-temp_m[,,i]
	temp_i<-temp_i[order(temp_i[,1]),]
	plot(temp_i[,1]/factor, temp_i[,2], ylim=c(0,1), xlab="", ylab="", type="l",  cex.axis=0.8, cex.main=0.8, main=names(Data)[i])
	}
	dev.off()
	temp[[1]]<-temp_m
    names(temp)[1]<-model_name
	}
	}

#########
## GBM ##
#########
if(class(models@fits[[1]][[1]])[1]=="lgb.Booster"){
temp<-list() # list where results will be stored

# Retrieve model fit
model<-models@fits[[1]][[1]]

# Set number of covariates
NbVar <- ncol(Data)

# Calculate responses
pts<-100 
parplot<-function(x){
ix<<-x
nsdm.partialPlot.gbm(model, as.matrix(Data), xname = colnames(Data)[ix], n.pt=pts)}

temp_m_l<-mclapply(1:length(names(Data)), parplot, mc.cores=ncores)

temp_m <- array(0, dim=c(pts, 2, NbVar), dimnames=list(NULL, c("Var", "Pred"), colnames(Data)))

for(i in 1:NbVar){
X1<-temp_m_l[i]
temp_m[,1,i]<-X1[[1]]$x
temp_m[,2,i]<-X1[[1]]$y
}

# Rescale using scale parameters
	if(!is.null(scaleparam)){
    for(i in 1:NbVar){
      temp_m[,1,i]<-temp_m[,1,i]*scaleparam$"scaled:scale"[names(Data)[i]] + scaleparam$"scaled:center"[names(Data)[i]]
	}}
	
# Plot
	if(plotting==TRUE){
    save_this_path<-paste(save_path, species_name, model_name, sep="/")
    suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
    pdf(paste0(save_this_path,"/", paste(species_name,model_name,sep="_"),".pdf"))
  
    ## plotting window
    W.width <- ceiling(sqrt(NbVar))
    W.height <- ceiling(NbVar/W.width)
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
	  par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n", axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves:", species_name, toupper(model_name), sep=" "),col="#4c57eb")
    par(mar = c(3,3,1.5,1))
	
	## do plot
	for(i in 1:dim(temp_m)[3]){
	temp_i<-temp_m[,,i]
	temp_i<-temp_i[order(temp_i[,1]),]
	plot(temp_i[,1]/factor, temp_i[,2], ylim=c(0,1), xlab="", ylab="", type="l",  cex.axis=0.8, cex.main=0.8, main=names(Data)[i])
	}
	dev.off()
	temp[[1]]<-temp_m
    names(temp)[1]<-model_name
	}
	}

########
## RF ##
########
if(class(models@fits[[1]][[1]])[2]=="randomForest"){
temp<-list() # list where results will be stored

# Retrieve model fit
model<-models@fits[[1]][[1]]

# Set number of covariates
NbVar <- ncol(Data)

# Calculate responses 
parplot<-function(x){
ix<<-x
randomForest::partialPlot(model, Data, x.var = colnames(Data)[ix], n.pt=25, plot=FALSE, which.class="1")}

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

temp_m_l<-mclapply(1:length(names(Data)), parplot, mc.cores=ncores)

temp_m <- array(0, dim=c(25, 2, NbVar), dimnames=list(NULL, c("Var", "Pred"), colnames(Data)))

for(i in 1:NbVar){
X1<-temp_m_l[i]
temp_m[,1,i]<-X1[[1]]$x
temp_m[,2,i]<-logit2prob(X1[[1]]$y)
}

# Rescale using scale parameters
	if(!is.null(scaleparam)){
    for(i in 1:NbVar){
      temp_m[,1,i]<-temp_m[,1,i]*scaleparam$"scaled:scale"[names(Data)[i]] + scaleparam$"scaled:center"[names(Data)[i]]
	}}
	
# Plot
	if(plotting==TRUE){
    save_this_path<-paste(save_path, species_name, model_name, sep="/")
    suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
    pdf(paste0(save_this_path,"/", paste(species_name,model_name,sep="_"),".pdf"))
  
       ## plotting window
    W.width <- ceiling(sqrt(NbVar))
    W.height <- ceiling(NbVar/W.width)
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
	  par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n", axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves:", species_name, toupper(model_name), sep=" "),col="#4c57eb")
    par(mar = c(3,3,1.5,1))
	
	## do plot
	for(i in 1:dim(temp_m)[3]){
	temp_i<-temp_m[,,i]
	temp_i<-temp_i[order(temp_i[,1]),]
	plot(temp_i[,1]/factor, temp_i[,2], ylim=c(0,1), xlab="", ylab="", type="l",  cex.axis=0.8, cex.main=0.8, main=names(Data)[i])
	}
	dev.off()
	temp[[1]]<-temp_m
    names(temp)[1]<-model_name
	}
	}

########
## Ranger ##
########
if (inherits(models@fits[[1]][[1]], "ranger")) {
  temp <- list()
  model <- models@fits[[1]][[1]]
  NbVar <- ncol(Data)
  Data_sub <- Data[sample(nrow(Data), min(1000, nrow(Data))), ]

  
  # Custom partial dependence function
parplot_fast <- function(ix) {
  var_name <- colnames(Data_sub)[ix]
  grid_vals <- seq(min(Data_sub[[ix]]), max(Data_sub[[ix]]), length.out = 25)
  
  # Repeat Data_sub 25 times and modify one column
  repeated_data <- Data_sub[rep(1:nrow(Data_sub), times = 25), ]
  repeated_data[[var_name]] <- rep(grid_vals, each = nrow(Data_sub))
  
  preds <- predict(model, data = repeated_data)$predictions
  preds_matrix <- matrix(preds, ncol = 25)
  avg_preds <- colMeans(preds_matrix)
  
  list(x = grid_vals, y = avg_preds)
}

    # Apply in parallel
temp_m_l <- mclapply(1:NbVar, parplot_fast, mc.cores = ncores)
  
  # Store results
  temp_m <- array(0, dim = c(25, 2, NbVar), dimnames = list(NULL, c("Var", "Pred"), colnames(Data)))
  for (i in 1:NbVar) {
    temp_m[, 1, i] <- temp_m_l[[i]]$x
    temp_m[, 2, i] <- temp_m_l[[i]]$y
  }
  
  # Rescale using provided scale parameters
  if (!is.null(scaleparam)) {
    for (i in 1:NbVar) {
      temp_m[, 1, i] <- temp_m[, 1, i] * scaleparam$`scaled:scale`[names(Data)[i]] + 
        scaleparam$`scaled:center`[names(Data)[i]]
    }
  }
  
  # Plotting
  if (plotting == TRUE) {
    save_this_path <- file.path(save_path, species_name, model_name)
    suppressWarnings(dir.create(save_this_path, recursive = TRUE))
    pdf(file.path(save_this_path, paste0(species_name, "_", model_name, ".pdf")))
    
    # Layout
    W.width <- ceiling(sqrt(NbVar))
    W.height <- ceiling(NbVar / W.width)
    mat <- matrix(c(rep(1, W.width), 1:(W.height * W.width) + 1), ncol = W.width, byrow = TRUE)
    layout(mat, widths = rep(1, W.width), heights = c(0.3, rep(1, W.height)))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = "n", axes = FALSE)
    polygon(c(-2, -2, 2, 2), c(-2, 2, 2, -2), col = "#f5fcba", border = NA)
    text(0.5, 0.8, labels = paste("Response curves:", species_name, toupper(model_name)), 
         cex = 1.6, col = "#4c57eb", pos = 1)
    par(mar = c(3, 3, 1.5, 1))
    
    # Draw curves
    for (i in 1:dim(temp_m)[3]) {
      temp_i <- temp_m[,,i]
      temp_i <- temp_i[order(temp_i[,1]),]
      plot(temp_i[,1] / factor, temp_i[,2], ylim = c(0, 1), xlab = "", ylab = "", 
           type = "l", main = colnames(Data)[i], cex.axis = 0.8, cex.main = 0.8)
    }
    
    dev.off()
  }
  
  temp[[1]] <- temp_m
  names(temp)[1] <- model_name
}
}


##############
## GLM, GAM ##
##############

if(class(models@fits[[1]][[1]])[2]!="randomForest"  & class(models@fits[[1]][[1]])[1]!="lgb.Booster" & class(models@fits[[1]][[1]])[1]!="maxnet"){
temp<-list() # list where results will be stored

for(m in 1:length(models@fits)){ # Loop over models (1 in regular cases but much more in esm settings)

# Duplicate Data for loop reuse
Data2<-Data
  
# Retrieve model fit
model<-models@fits[[m]][[1]]		  
  
# In case of esm model uses its name its index for naming and subset Data
if(model_name=="esm"){
model_nameu<-names(models@fits)[m]
Data2<-Data2[,(gsub(".*[(]([^.]+)[,].*", "\\1", rownames(model$R)[c(2,5)]))]
}else{
model_nameu<-model_name}

# Set number of covariates
NbVar <- ncol(Data2)

# Calculate responses
temp_m <- array(0, dim=c(nrow(Data2), 2, NbVar), dimnames=list(NULL, c("Var", "Pred"), colnames(Data2)))
Xp  <- as.data.frame(matrix(NA, ncol=NbVar, nrow=nrow(Data2), dimnames=list(NULL, colnames(Data2))))
	
for(i in 1:NbVar){
Xp[,i] <- mean(Data2[,i])
}   

    for(i in 1:NbVar){
    xr <- range(Data2[,i])
    Xp1 <- Xp
    Xp1[,i] <- seq(xr[1], xr[2],  len=nrow(Data2))
    
	  if(class(model)[1]=="glm") Xf <- predict(model, as.data.frame(Xp1), type="response")
    if(c("gam") %in% class(model)) Xf <- predict(model, as.data.frame(Xp1), type="response")
    if(class(model)[1]=="maxnet") Xf <-  predict(model, as.data.frame(Xp1), type="cloglog") 

    temp_m[,1,i] <-Xp1[ ,i]; temp_m[,2,i] <- Xf
	}
	
# Rescale using scale parameters
	if(!is.null(scaleparam)){
    for(i in 1:NbVar){
      temp_m[,1,i]<-temp_m[,1,i]*scaleparam$"scaled:scale"[names(Data2)[i]] + scaleparam$"scaled:center"[names(Data2)[i]]
	}}
	
# Plot
	if(plotting==TRUE){
    if(model_name!="esm") save_this_path<-paste(save_path, species_name, model_nameu, sep="/")
	if(model_name=="esm") save_this_path<-paste(save_path, species_name, model_name, model_nameu, sep="/")
    suppressWarnings(dir.create(save_this_path,  recursive = TRUE))
    pdf(paste0(save_this_path,"/", paste(species_name,model_nameu,sep="_"),".pdf"))
  
    ## plotting window
    W.width <- ceiling(sqrt(NbVar))
    W.height <- ceiling(NbVar/W.width)
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n", axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves:", species_name, toupper(model_nameu), sep=" "),col="#4c57eb")
    par(mar = c(3,3,1.5,1))
	
	## do plot
	for(i in 1:dim(temp_m)[3]){
	temp_i<-temp_m[,,i]
	temp_i<-temp_i[order(temp_i[,1]),]
	plot(temp_i[,1]/factor, temp_i[,2], ylim=c(0,1), xlab="", ylab="", type="l",  cex.axis=0.8, cex.main=0.8, main=names(Data2)[i])
	}
	dev.off()
	temp[[m]]<-temp_m
	names(temp)[m]<-model_nameu
	}
}
}

if(model_name=="esm")return(temp)
if(model_name!="esm")return(temp[[1]])
}