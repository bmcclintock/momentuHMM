#' Get starting values (Par) from momentuHMM object returned by fitHMM
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return A list of parameter values that can be used as starting values (Par) in fitHMM
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' Par <- getPar(m)
#'
#' @export
getPar<-function(m){
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  cons <- m$conditions$cons
  logitcons <- m$conditions$logitcons
  bounds <- m$conditions$bounds
  
  parindex <- c(0,cumsum(unlist(lapply(DM,ncol)))[-length(DM)])
  names(parindex) <- distnames
  
  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,bounds,DM)
  wpar <- m$mod$estimate
  
  Par <- list()
  for(i in distnames){
    if(!is.numeric(bounds[[i]])){
      bounds[[i]] <- gsub(i,"",bounds[[i]],fixed=TRUE)
    }
    Par[[i]] <- c(w2nDM(wpar[parindex[[i]]+1:ncol(DM[[i]])],bounds[[i]],DM[[i]],cons[[i]],p$boundInd[[i]],logitcons[[i]])$Par)
  }
  Par
}