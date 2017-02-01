#' Get starting values (Par) on working scale based on design matrix
#'
#'
#' @return A list of parameter values (Par) that can be used as starting values in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}
#'
#' @examples
#' # m is a moveHMMext object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' Par <- getPar(m)
#'
#' @export
getParDM<-function(data,nbStates,dist,
                 Par,
                 estAngleMean=NULL,
                 DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL,stateNames=NULL){
  
  if(nbStates<1) stop('nbStates must be >0')
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
  distnames<-names(dist)
  if(any(is.na(match(distnames,names(data))))) stop(paste0(distnames[is.na(match(distnames,names(data)))],collapse=", ")," not found in data")
  if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
  Par <- Par[distnames]
  
  zeroInflation <- vector('list',length(distnames))
  names(zeroInflation) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% zeroInflationdists){
      if(length(which(data[[i]]==0))>0) {
        zeroInflation[[i]]<-TRUE
      }
      else 
        zeroInflation[[i]]<-FALSE
    }
    else zeroInflation[[i]]<-FALSE
  }
  
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,zeroInflation,DM,userBounds,cons,workcons,stateNames)
  
  DMinputs<-getDM(data,inputs$DM,dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,inputs$cons,inputs$workcons,zeroInflation,FALSE)
  fullDM <- DMinputs$fullDM
  DMind <- DMinputs$DMind
  cons <- DMinputs$cons
  workcons <- DMinputs$workcons
  
  wpar <- Par
  for(i in distnames){
    if(!is.null(DM[[i]])){
      if(DMind[[i]]){
        if(length(wpar[[i]])!=nrow(fullDM[[i]])) stop('Par$',i,' should be of length ',nrow(fullDM[[i]]))
        bounds<-inputs$p$bounds[[i]]
        if(any(wpar[[i]]<=bounds[,1] | wpar[[i]]>=bounds[,2])) stop('Par$',i,' must be within parameter bounds')
        bndInd <- which(!duplicated(getboundInd(fullDM[[i]])))
        a<-bounds[bndInd,1]
        b<-bounds[bndInd,2]
        par <- wpar[[i]][bndInd]
        if(any(wpar[[i]]!=par[getboundInd(fullDM[[i]])])) stop('Par$',i,' values are not consistent with DM$',i)
        piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
        ind1<-which(piInd)
        ind2<-which(!piInd)
        
        p<-numeric(length(bndInd))
        if(length(ind1)) p[ind1] <- (solve(unique(fullDM[[i]])[ind1,ind1],tan(par[ind1]/2))-workcons[[i]][ind1])^(1/cons[[i]][ind1])
        
        ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
        ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
        ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
        
        if(length(ind21)) p[ind21]<-solve(unique(fullDM[[i]])[ind21,ind21],log(par[ind21]-a[ind21])-workcons[[i]][ind21])^(1/cons[[i]][ind21])
        if(length(ind22)) p[ind22]<-solve(unique(fullDM[[i]])[ind22,ind22],logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))-workcons[[i]][ind22])^(1/cons[[i]][ind22])
        if(length(ind23)) p[ind23]<-solve(unique(fullDM[[i]])[ind23,ind23],-log(-par[ind23]+b[ind23])-workcons[[i]][ind23])^(1/cons[[i]][ind23])
        
        if(any(!is.finite(p))) stop(i," working scale parameters are not finite. Check natural parameter values, bounds, and constraints.")
        wpar[[i]]<-p
      } else stop('sorry, design matrices with individual covariates are not supported by getParDM')
    }
  }
  wpar
}