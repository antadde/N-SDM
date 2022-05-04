#' Feed the predict functions depending on model type
#'
#' @author Philipp Brun (philipp.brun@wsl.ch)

nsdm.prd=function(mod,tst){

  # Generate probabilistic precitions
  if("maxnet"%in%class(mod)){

    pred<-nsdm.df_or_rast(mod,
                     nwdat=tst,
                     type="logistic")
    
  } else if(any(c("glm","ranger")%in%class(mod))){

    pred<-nsdm.df_or_rast(mod=mod,
                     nwdat=tst,
                     type="response")

  } else if("gbm"%in%class(mod)){

    pred<-nsdm.df_or_rast(mod,
                     nwdat=tst,
                     n.trees=mod$n.trees,
                     type="response")

  } else if("lgb.Booster"%in%class(mod)){
    
    pred<-nsdm.df_or_rast(mod=mod,
                     nwdat=tst)
    
  } else if("randomForest"%in%class(mod)){

    pred<-nsdm.df_or_rast(mod,
                     nwdat=tst,
                     type="prob")
  }

  # Convert to numeric
  if(class(tst)=="data.frame"){
    pred<-as.numeric(pred)
  }

  return(pred)

}
