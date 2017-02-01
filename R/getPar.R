#' Get starting values from momentuHMM or momentuHMMMI object returned by fitHMM or MI_summary
#'
#' @param m A \code{momentuHMM} or \code{momentuHMMMI} object.
#'
#' @return A list of parameter values (Par, beta, delta) that can be used as starting values in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' Par <- getPar(m)
#'
#' @export
getPar<-function(m){
  
  if(!is.momentuHMM(m) & !is.momentuHMMMI(m))
    stop("'m' must be a momentuHMM or momentuHMMMI object (as output by fitHMM or MI_summary)")
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  
  if(is.momentuHMMMI(m)){
    m$mle<-lapply(m$Par[distnames],function(x) x$est)
    beta <- unname(m$Par$beta$beta$est)
    delta <- unname(m$Par$real$delta$est)
  } else {
    beta <- unname(m$mle$beta)
    delta <- unname(m$mle$delta)
  }
  
  Par <- list()
  for(i in distnames){
    if(is.null(DM[[i]])){
      par <- unname(m$mle[[i]])
      if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]]) 
        par <- par[-1,]
      par <- c(t(par))
    } else par <- unname(m$CI_beta[[i]]$est)
    Par[[i]] <- par
  }
  list(Par=Par,beta=beta,delta=delta)
}