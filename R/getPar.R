#' Get starting values from momentuHMM, miHMM, or miSum object returned by fitHMM, MIfitHMM, or MIpool
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return A list of parameter values (Par, beta, delta) that can be used as starting values in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}
#'
#' @seealso \code{\link{getPar0}}, \code{\link{getParDM}}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' Par <- getPar(m)
#'
#' @export
getPar<-function(m){
  
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  m <- delta_bc(m)
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  
  parCount<- lapply(m$conditions$fullDM,ncol)
  for(i in distnames[unlist(m$conditions$circularAngleMean)]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount)))
  names(parindex) <- c(distnames,"beta")
  
  Par <- list()
  if(is.miSum(m)){
    m$mle<-lapply(m$Par$real[distnames],function(x) x$est)
    for(i in distnames){
      if(!is.null(DM[[i]]))
        m$mle[[i]]<-nw2w((m$Par$beta[[i]]$est-m$conditions$workcons[[i]])^(1/m$conditions$cons[[i]]),m$conditions$workBounds[[i]])
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        m$mle[[i]]<-m$mle[[i]][-1,]
      Par[[i]] <- c(t(unname(m$mle[[i]])))
    }
    beta <- unname(nw2w(m$Par$beta$beta$est,m$conditions$workBounds$beta))
    if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
      delta <- unname(m$Par$real$delta$est[1,])
    } else {
      delta <- unname(nw2w(m$Par$beta$delta$est,m$conditions$workBounds$delta))
    }
  } else {
    for(i in distnames){
      if(is.null(DM[[i]])){
        par <- unname(m$mle[[i]])
        if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]]) 
          par <- par[-1,]
        par <- c(t(par))
      } else par <- unname(m$mod$estimate[parindex[[i]]+1:parCount[[i]]])#unname(nw2w((m$CIbeta[[i]]$est-m$conditions$workcons[[i]])^(1/m$conditions$cons[[i]]),m$conditions$workBounds[[i]]))
      Par[[i]] <- par
    }
    beta <- unname(matrix(m$mod$estimate[parindex[["beta"]]+1:length(m$mle$beta)],nrow(m$mle$beta),ncol(m$mle$beta)))#unname(nw2w(m$mle$beta,m$conditions$workBounds$beta))
    if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
      delta <- unname(m$mle$delta[1,])
    } else {
      delta <- unname(matrix(m$mod$estimate[parindex[["beta"]]+length(m$mle$beta)+1:length(m$CIbeta$delta$est)],nrow(m$CIbeta$delta$est),ncol(m$CIbeta$delta$est)))#unname(nw2w(m$CIbeta$delta$est,m$conditions$workBounds$delta))
    }
  }
  list(Par=Par,beta=beta,delta=delta)
}