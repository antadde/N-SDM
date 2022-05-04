#' nsdm.bigextract
#'
#' Parallel loading and spatiotemporal extraction of large raster sets
#'
#' @param n Number of background points (pseudoabsences) to be generated
#' @param pres Spatial point object containing species occurence data
#' @param taxon A character string for the target species
#' @param type A character string indicating the type of species data ("po"=presence only; "pa"=presence absence)
#' @param rst_ref Reference raster used in subsequent analyses
#' @param set_max_npres_to_nabs Logical (TRUE FALSE) for setting the maximal number of presences to the number of background points
#'
#' @return A nsdm.pseudoabsences object
#' @author Antoine Adde (aadde@unil.ch) and Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.pseuabs<-function(n=10000,
                       pres,
                       taxon=character(),
					   type="po",
					   rst_ref,
                       set_max_npres_to_nabs=TRUE){

  ### ------------------------
  ### generate nsdm.pseudoabsences object and add meta info
  ### ------------------------
  out<-preva.meta(type="pseudoabsence")
  out@meta$type="random"
  out@meta$taxon=taxon
  call=match.call()
  out@call<-call

  ### ------------------------
  ### Prepare presences
  ### ------------------------
      if(type=="pa"){
	  set_max_npres_to_nabs=FALSE
	  pres=pres[pres$pa==1,]}
      # set_max_npres_to_nabs if requested
	  if(set_max_npres_to_nabs && length(pres)>n){
	  ## stratum indicator
	  d<-NULL
	  d <- data.frame(id=1:length(pres), coordinates(pres))
	  if(d[1,2]>10000){ d$group<- round(d[,2]/10000*d[,3]/10000) 
      } else {
	  d$group <- round(d[,2]*d[,3])}
	  ## sample selection
      outs <- nsdm.stratified(d, "group", (n+100)/nrow(d))
	  pres<-pres[outs$id,]
	  }

  ### ------------------------
  ### Prepare (pseudo)absences
  ### ------------------------
  if(type=="po"){
  # A- Sample pseudoabsences
  nona<-which(!is.na(raster::values(rst_ref)))
  rnd.pts<-sample(nona, size=n*1.5)
  crds_abs<-coordinates(rst_ref)[rnd.pts,]
  abs<-SpatialPoints(coords=crds_abs, proj4string =rst_ref@crs)
  # B- Subsample to n+200 (some points will be dropped after NA covariate cleaning)
   if(length(abs)>n+200) abs<-abs[sample(1:length(abs),n+200),]
  } else if(type=="pa"){
  abs<-pres[pres$pa==0,]}
    
  ### ------------------------
  ### Prepare output
  ### ------------------------
    out@pa<-c(rep(1,length(pres)),
              rep(0,length(abs)))
	out@years<-as.numeric(c(pres$year, rep(NA, length(abs))))
    out@env_vars=data.frame()
    out@xy=rbind(coordinates(pres),coordinates(abs))
  
  # return
  return(out)
}
