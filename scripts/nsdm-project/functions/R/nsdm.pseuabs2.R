#' nsdm.pseuabs2
#'
#' Prepare background (pseudoabsences) points
#'
#' @param n Number of background points to be generated
#' @param pres Spatial point object containing species occurrence data
#' @param taxon A character string for the target species
#' @param type A character string indicating the type of species data ("po"=presence only; "pa"=presence absence)
#' @param rst_ref Reference raster used in subsequent analyses
#' @param set_max_npres_to_nabs Logical (TRUE FALSE) for setting the maximal number of presences to the number of background points
#' @param samp_area A character string to the basename of the polygon shapefile used for sampling the background points
#'
#' @return A nsdm.pseudoabsences object
#' @author Antoine Adde (aadde@unil.ch) and Philipp Brun (philipp.brun@wsl.ch)
#' @export
nsdm.pseuabs2<-function(n=10000,
                       pres,
                       taxon=character(),
					             type="po",
					             rst_ref,
                       set_max_npres_to_nabs=TRUE,
					             samp_area=NULL){

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
    abs=pres[pres$pa==0,]
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
  # A0- Refine target sampling area if needed
  if(length(samp_area)>0){
  shp_samp_area<-readOGR(dsn = paste0(w_path,"data/",project,"/masks/"), layer = gsub(".shp", "", samp_area))
  rst_ref <- crop(rst_ref, extent(shp_samp_area))
  rst_ref <- mask(rst_ref, shp_samp_area)
  }
  # A1- Sample pseudoabsences
  nona<-which(!is.na(raster::values(rst_ref)))
  rnd.pts<-sample(nona, size=n*1.5)
  crds_abs<-coordinates(rst_ref)[rnd.pts,]
  abs<-SpatialPointsDataFrame(coords=crds_abs, data=data.frame(X=crds_abs[,1],Y=crds_abs[,2]), proj4string =rst_ref@crs)
  # B- Subsample to n+200 (some points will be dropped after NA covariate cleaning)
  if(length(abs)>n+200) abs<-abs[sample(1:length(abs),n+200),]}

  ### ------------------------
  ### Prepare output
  ### ------------------------
  out@pa<-c(rep(1,nrow(pres)),
            rep(0,nrow(abs)))
  out@years<-as.numeric(c(pres$year, rep(NA, nrow(abs))))
  out@env_vars=data.frame()
  out@xy=as.matrix(rbind(data.frame(X=pres$X, Y=pres$Y), data.frame(X=abs$X, Y=abs$Y)))
  
  # return
  return(out)
}
