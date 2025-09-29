#' nsdm.prd
#'
#' Feed the predict functions depending on model type
#'
#' @author Antoine Adde (antoine.adde@eawag.ch)
#' @export

nsdm.prd <- function(mod, tst) {
  
  if ("maxnet" %in% class(mod)) {
    
    pred <- predict(mod,
                    newdata = tst,
                    type = "cloglog")
    
  } else if (any(c("glm") %in% class(mod))) {
    
    pred <- predict(mod,
                    newdata = tst,
                    type = "response")
    
  } else if ("lgb.Booster" %in% class(mod)) {
    
    pred <- predict(mod,
                    data = as.matrix(tst))
    
  } else if ("randomForest" %in% class(mod)) {
    
    pred <- predict(mod,
                    newdata = tst,
                    type = "prob")[, 2]  # Probability for class 1
      
  }

  pred <- round(as.numeric(pred), 2)
  
  return(pred)
}
