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
  
  m$conditions$optInd <- numeric() # extra hack needed for bc
  m <- delta_bc(m)
  m$conditions$optInd <- NULL # extra hack needed for bc
  
  if(!is.null(m$mod$hessian) & inherits(m$CIbeta,"error")){
    m$mod$hessian <- NULL
    m$CIbeta <- tryCatch(CIbeta(m),error=function(e) e)
  }
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  
  parCount<- lapply(m$conditions$fullDM,ncol)
  for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount)))
  names(parindex) <- c(distnames,"beta")
  
  if(is.null(m$conditions$formulaPi)) {
    formPi <- ~1
  } else formPi <- m$conditions$formulaPi
  
  if(is.null(m$conditions$formulaDelta)) {
    formDelta <- ~1
  } else formDelta <- m$conditions$formulaDelta
  
  Par <- list()
  beta <- delta <- NULL
  if(is.miSum(m)){
    m$mle<-lapply(m$Par$real[distnames],function(x) x$est)
    for(i in distnames){
      if(!is.null(DM[[i]]))
        m$mle[[i]]<-nw2w((m$Par$beta[[i]]$est-m$conditions$workcons[[i]])^(1/m$conditions$cons[[i]]),m$conditions$workBounds[[i]])
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        m$mle[[i]]<-m$mle[[i]][-1,]
      Par[[i]] <- c(t(unname(m$mle[[i]])))
    }
    if(nbStates>1){
      beta <- unname(nw2w(m$Par$beta$beta$est,m$conditions$workBounds$beta))
      if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
        delta <- unname(m$Par$real$delta$est[1:m$conditions$mixtures,])
      } else {
        delta <- unname(nw2w(m$Par$beta$delta$est,m$conditions$workBounds$delta))
      }
      if(!is.null(m$conditions$recharge)){
        beta <- list(beta=beta,g0=nw2w(m$Par$beta$g0$est,m$conditions$workBounds$g0),theta=nw2w(m$Par$beta$theta$est,m$conditions$workBounds$theta))
      }
      if(m$conditions$mixtures>1){
        if(!length(attr(terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
          pie <- unname(m$Par$real$pi$est)
        } else {
          pie <- unname(nw2w(m$Par$beta$pi$est,m$conditions$workBounds$pi))
        }
      }
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
    if(nbStates>1){
      beta <- unname(matrix(m$mod$estimate[parindex[["beta"]]+1:length(m$mle$beta)],nrow(m$mle$beta),ncol(m$mle$beta)))#unname(nw2w(m$mle$beta,m$conditions$workBounds$beta))
      if(m$conditions$stationary & length(attr(terms.formula(m$conditions$formula),"term.labels"))>0){
        delta <- unname(m$mle$delta)
      } else if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
        delta <- unname(m$mle$delta[seq(1,nrow(m$mle$delta),nrow(m$mle$delta)/m$conditions$mixtures),])
      } else {
        delta <- unname(matrix(m$mod$estimate[parindex[["beta"]]+length(m$mle$beta)+length(m$CIbeta$pi$est)+1:length(m$CIbeta$delta$est)],nrow(m$CIbeta$delta$est),ncol(m$CIbeta$delta$est)))#unname(nw2w(m$CIbeta$delta$est,m$conditions$workBounds$delta))
      }
      if(!is.null(m$conditions$recharge)){
        beta <- list(beta=beta,g0=m$mod$estimate[parindex[["beta"]]+length(m$mle$beta)+length(m$CIbeta$pi$est)+length(m$CIbeta$delta$est)+1:length(m$mle$g0)],theta=m$mod$estimate[parindex[["beta"]]+length(m$mle$beta)+length(m$CIbeta$pi$est)+length(m$CIbeta$delta$est)+length(m$mle$g0)+1:length(m$mle$theta)])
      }
      if(m$conditions$mixtures>1){
        if(!length(attr(terms.formula(formPi),"term.labels")) & is.null(m$conditions$formulaPi)){
          pie <- unname(m$mle$pi[1,])
        } else {
          pie <- unname(matrix(m$mod$estimate[parindex[["beta"]]+length(m$mle$beta)+1:length(m$CIbeta$pi$est)],nrow(m$CIbeta$pi$est),ncol(m$CIbeta$pi$est)))#unname(nw2w(m$CIbeta$pi$est,m$conditions$workBounds$pi))
        }
      }
    }
  }
  if(m$conditions$mixtures>1){
    if(!is.list(beta)) beta <- list(beta=beta,pi=pie)
    else beta$pi <- pie
  }
  
  out <- list(Par=Par,beta=beta,delta=delta)
  
  # map hierarchical beta and delta
  if(is.momentuHierHMM(m)){
    if(is.list(beta)){
      Pi <- beta$pi
    } else {
      Pi <- NULL
    }
    inputHierHMM <- formatHierHMM(m$data,m$conditions$hierStates,m$conditions$hierDist,hierBeta=NULL,hierDelta=NULL,m$conditions$hierFormula,m$conditions$hierFormulaDelta,m$conditions$mixtures)
    hier <- mapHier(beta,Pi,delta,m$conditions$hierBeta,m$conditions$hierDelta,inputHierHMM$hFixPar,inputHierHMM$hBetaCons,inputHierHMM$hDeltaCons,m$conditions$hierStates,inputHierHMM$newformula,m$conditions$formulaDelta,inputHierHMM$data,m$conditions$mixtures,inputHierHMM$recharge)
    beta <- hier$hierBeta
    delta <- hier$hierDelta
    out <- list(Par=Par, hierBeta=beta, hierDelta=delta)
  }
  return(out)
}