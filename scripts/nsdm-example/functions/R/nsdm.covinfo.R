#' nsdm.covinfo
#'
#' Generate covariate information table based on covariate paths/names
#'
#' @param cov_path A character string for the root path where covariates are stored
#' @param save_path A character string for the path to which write the .xlsx covariate information table
#' @param dyna_cov A character vector containing the names of datasets equipped with temporally dynamic covariates
#' @param focal_cat A character vector containing the names of categories equipped with temporally dynamic covariates
#'
#' @return An xlsx covariate information table
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.covinfo<-function(
cov_path,
save_path,
time_cov, 
focal_cov){ 

### List all available layers and extract their information
setwd(cov_path)    
pred_list<-list.files(".", pattern='.rds', recursive = TRUE)
full_pred_list <- list.files(cov_path, pattern='.rds', recursive = TRUE, full.names = TRUE) 
pred_table <- c()

## Generic information (level, category, dataset)
pred_table$level <- unlist(lapply(strsplit(pred_list, split='/'),  `[[`, 1))
pred_table$category <- unlist(lapply(strsplit(pred_list, split='/'),  `[[`, 2))
pred_table$dataset <- unlist(lapply(strsplit(pred_list, split='/'),  `[[`, 3))
pred_table$cada<-paste(pred_table$category, pred_table$dataset, sep="_")

## Non-generic information
for(i in 1:length(pred_list)){ 
   # Dynamic covariates
    if(pred_table$cada[i] %in% time_cov){
    pred_table$period[i] <-  strsplit(pred_list[i], split='/')[[1]][4]
   # Do we have future layers for this covariate?
   if(pred_table$period[i]=="future"){  
     pred_table$scale[i] <- strsplit(pred_list[i], split='/')[[1]][5]
     pred_table$start_year[i] <- strsplit(strsplit(pred_list[i], split='/')[[1]][6], split='_')[[1]][1]
     pred_table$end_year[i] <- strsplit(strsplit(pred_list[i], split='/')[[1]][6], split='_')[[1]][2]
     pred_table$scenario[i] <- strsplit(pred_list[i], split='/')[[1]][7]
	 pred_table$variable[i] <- strsplit(pred_list[i], split='/')[[1]][8]
	 pred_table$attribute[i]<- str_match(pred_list[i], paste0(pred_table$variable[i],"_\\s*(.*?)\\s*.rds"))[2]
     } else {
     pred_table$scale[i] <- strsplit(pred_list[i], split='/')[[1]][5]
     pred_table$start_year[i] <- strsplit(strsplit(pred_list[i], split='/')[[1]][6], split='_')[[1]][1]
     pred_table$end_year[i] <- strsplit(strsplit(pred_list[i], split='/')[[1]][6], split='_')[[1]][2]
	 pred_table$scenario[i]<-NA
     pred_table$variable[i] <- strsplit(pred_list[i], split='/')[[1]][7]
	 pred_table$attribute[i]<- str_match(pred_list[i], paste0(pred_table$variable[i],"_\\s*(.*?)\\s*.rds"))[2]   
    } 
	} else {
    # Static covariates
    pred_table$period[i] <- 'present'
    pred_table$start_year[i] <- NA
    pred_table$end_year[i] <- NA
	pred_table$scenario[i]<-NA
    pred_table$scale[i] <- strsplit(pred_list[i], split='/')[[1]][4]
    pred_table$variable[i] <- strsplit(pred_list[i], split='/')[[1]][5]
	pred_table$attribute[i] <- str_match(pred_list[i],  paste0("_",pred_table$variable[i],"_\\s*(.*?)\\s*.rds"))[2]
  }
   # Add focals for equipped datasets
   if(pred_table$cada[i] %in% focal_cov){
    focal<-sub("^.+_", "", file_path_sans_ext(full_pred_list[i]))
    suppressWarnings(true_focal<-as.numeric(focal))
    if(focal=="cavg")true_focal<-"cavg"
    pred_table$focal[i]<-true_focal
  }else{
    pred_table$focal[i]<-NA
  }
    
  # Add file path
   pred_table$file[i]<-full_pred_list[i]
}

# Format and clean
pred_table<-data.frame(pred_table)
pred_table$attribute<-unlist(lapply(1:nrow(pred_table), function(i){gsub(paste0("_",pred_table$focal[i]), "", pred_table$attribute[i])})) # remove focals from attributes column
pred_table$attribute[which(unlist(lapply(1:nrow(pred_table), function(z){pred_table$focal[z] == pred_table$attribute[z]}))==TRUE)]<-NA # remove focals from attributes column
pred_table[is.na(pred_table)] = "NA"

# Save covariate table
row.names(pred_table)<-NULL
write_xlsx(pred_table, paste0(save_path,"predictors-available.xlsx"))
}