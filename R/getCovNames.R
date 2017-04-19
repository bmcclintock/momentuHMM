#' Get names of any covariates used in probability distribution parameters
#' 
#' @param m \code{\link{momentuHMM}} object
#' @param p list returned by \code{\link{parDef}}
#' @param distname Name of the data stream
#' 
#' @return A list of:
#' \item{DMterms}{Names of all covariates included in the design matrix for the data stream}
#' \item{DMpartems}{A list of the names of all covariates for each of the probability distribution parameters}
#' 
#' @importFrom stats terms
getCovNames<-function(m,p,distname){
  nbStates<-length(m$stateNames)
  DMparterms<-list()
  DMterms<-character()
  if(!m$conditions$DMind[[distname]]){
    if(!is.list(m$conditions$DM[[distname]])){
      for(j in 1:length(p$parNames[[distname]])){
        DMparterms[[p$parNames[[distname]][j]]] <- vector('list',nbStates)
        for(jj in 1:nbStates){
          DMparterms[[p$parNames[[distname]][j]]][[jj]]<-unique(m$conditions$DM[[distname]][(j-1)*nbStates+jj,][suppressWarnings(which(is.na(as.numeric(m$conditions$DM[[distname]][(j-1)*nbStates+jj,]))))])
        }
      }
      DMterms<-unique(m$conditions$DM[[distname]][suppressWarnings(which(is.na(as.numeric(m$conditions$DM[[distname]]))))])
    } else {
      m$conditions$DM[[distname]]<-m$conditions$DM[[distname]][p$parNames[[distname]]]
      for(j in 1:length(p$parNames[[distname]])){
        DMparterms[[p$parNames[[distname]][j]]] <- vector('list',nbStates)
        tmpparnames <- rownames(attr(terms(m$conditions$DM[[distname]][[p$parNames[[distname]][j]]]),"factors"))
        if(!is.null(tmpparnames)) {
          for(jj in 1:nbStates){
            DMparterms[[p$parNames[[distname]][j]]][[jj]]<-tmpparnames
          }
        }
        DMterms <- c(DMterms,unlist(DMparterms[[p$parNames[[distname]][j]]]))
      }
    }
    DMterms <- unique(DMterms)
    for(j in 1:length(p$parNames[[distname]])){
      for(jj in 1:nbStates){
        if(length(DMparterms[[p$parNames[[distname]][j]]][[jj]])){
          for(k in 1:length(DMparterms[[p$parNames[[distname]][j]]][[jj]])){
            DMparterms[[p$parNames[[distname]][j]]][[jj]][k]<-all.vars(as.formula(paste0("~",DMparterms[[p$parNames[[distname]][j]]][[jj]][k])))
          }
        }
        DMparterms[[p$parNames[[distname]][j]]][[jj]]<-unique(DMparterms[[p$parNames[[distname]][j]]][[jj]])
      }
    }
    for(j in 1:length(DMterms)){
      DMterms[j]<-all.vars(as.formula(paste0("~",DMterms[j])))
    }
    DMterms <- unique(DMterms)
  }
  list(DMterms=DMterms,DMparterms=DMparterms)
}