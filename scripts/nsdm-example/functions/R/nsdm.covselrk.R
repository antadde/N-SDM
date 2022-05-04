#' nsdm.covselrk
#'
#' Ranking of covariates after selection procedure
#'
#' @param embed An nsdm.embedsel output
#' @param species_name A character string indicating the name of the target taxon
#'
#' @return A ranked list of candidate covariates after selection procedure
#' @author Antoine Adde (aadde@unil.ch)
#' @export

nsdm.covselrk <- function(embed, species_name){
embed$species<-species_name

## Rank covariates selected commonly by the three algorithms (intersect)
# Check first if such covariate exists
if(length(which(table(embed$var) == length(unique(embed$model))))>0){
# if so, do intersect ranking
intersect.tmp<-embed[embed$var %in% names(which(table(embed$var) == length(unique(embed$model)))),]
intersect.tmp<-aggregate(intersect.tmp[,c("rank")], list(intersect.tmp$var), mean); colnames(intersect.tmp)<-c("var","rank")
intersect.sel<-data.frame(intersect.tmp[order(intersect.tmp$rank, decreasing = FALSE),], rank.f = 1:nrow(intersect.tmp), species=species_name)
# else skip to union ranking
} else { intersect.sel<-data.frame(var=NA, rank=0, rank.f=0, species=species_name)}

## Rank and add other covariates (union), if needed
union.tmp<-embed[embed$var%in%names(which(table(embed$var) < length(unique(embed$model)))),]
# Check first if union ranking is needed
if(nrow(union.tmp)>0){
union.tmp<-aggregate(union.tmp[,c("rank")], list(union.tmp$var), mean); colnames(union.tmp)<-c("var","rank")
union.sel.tmp<-data.frame(union.tmp[order(union.tmp$rank, decreasing = FALSE),], rank.f = (max(intersect.sel$rank.f+1)):(max(intersect.sel$rank.f)+nrow(union.tmp)),species=species_name)
union.sel<-rbind(intersect.sel,union.sel.tmp)
# else skip to results
} else {
union.sel<-intersect.sel
}
union.sel$type<-"union"
union.sel<-union.sel[complete.cases(union.sel), ]
# Return results
return(union.sel)
}