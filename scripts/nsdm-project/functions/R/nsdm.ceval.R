#' nsdm.ceval
#'
#' Model evaluation core function
#'
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)
#' @export

nsdm.ceval<-function(f,pa,tesdat,crit,tre=numeric()){

  # If there are any presences in the evaluation data
  if(any(pa==1)){

    if(crit=="maxTSS" && length(tre)==0){
      
      ordf=order(f,decreasing = T)
      
      rdf=f[ordf]
      opa=pa[ordf]
      cpa=cumsum(opa)
      
      # Reduce to sensible range
      nsns=which(cpa==max(cpa))
      nsns=nsns[-1]

      if(length(nsns)!=0){
           rdf=rdf[-nsns]
           opa=opa[-nsns]
           cpa=cpa[-nsns]
      }

      
      tsss=apply(cbind(1:length(cpa),cpa),1,function(x,y,z){
        out=tss(x[2],x[1]-x[2],z-x[2],y-x[1]-(z-x[2]))
        return(out)
      }, y=length(pa),z=sum(pa))
      
      tre=rdf[which.max(tsss)]

    } else if (crit=="" && length(tre) == 0){
      tre = mean(pa)
    }

    # Calculate threshold-dependent metrics
    bina=ifelse(f<tre,0,1)
    tb = table(factor(bina,levels=c("0","1")),factor(pa,levels=c("0","1")))
    tdep = unlist(nsdm.all.metrics(tb[2,2],tb[2,1],tb[1,2],tb[1,1]))

   	# Boyce
	  boyce_s_pe=try(ecospat.boyce(fit=f, obs=f[pa==1], PEplot = FALSE, rm.duplicate = TRUE, method = 'spearman'), TRUE)
	  boyce_p_pe=try(ecospat.boyce(fit=f, obs=f[pa==1], PEplot = FALSE, rm.duplicate = TRUE, method = 'pearson'), TRUE)
	  boyce_p_npe=try(ecospat.boyce(fit=f, obs=f[pa==1], PEplot = FALSE, rm.duplicate = FALSE, method = 'pearson'), TRUE)
	  boyce_s_npe=try(ecospat.boyce(fit=f, obs=f[pa==1], PEplot = FALSE, rm.duplicate = FALSE, method = 'spearman'), TRUE)

	  if(is.numeric(boyce_s_pe$cor)) { boyce_s_pe=boyce_s_pe$cor } else { boyce_s_pe=NA}
	  if(is.numeric(boyce_p_pe$cor)) {boyce_p_pe=boyce_p_pe$cor } else { boyce_p_pe=NA}
	  if(is.numeric(boyce_s_npe$cor)) { boyce_s_npe=boyce_s_npe$cor } else { boyce_s_npe=NA}
	  if(is.numeric(boyce_p_npe$cor)) { boyce_p_npe=boyce_p_npe$cor } else { boyce_p_npe=NA}
	  
    # AUC
    z=prediction(f,pa)
    auc=performance(z,measure="auc")@y.values[[1]]
    rmse=performance(z,measure="rmse")@y.values[[1]]
	
	# Somer's AUC
	aucS=(2*auc)-1 # rescale AUC to -1 +1

	# # Consensus score
	score_vec<-na.omit(c(aucS, boyce_s_pe, tdep["TSS"]))
	score<-sum(score_vec)/length(score_vec)

    # Return results
    weg=c(auc, aucS, rmse, boyce_s_pe, boyce_p_pe, boyce_s_npe, boyce_p_npe,  score, tre, tdep)
    names(weg)[1:9]=c("AUC", "AUC_S", "RMSE", "Boyce_S_PE", "Boyce_P_PE", "Boyce_S_NPE", "Boyce_P_NPE", "Score", "threshold")
    
    return(unlist(weg))
  }
}
