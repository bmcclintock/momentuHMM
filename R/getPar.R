#' Get starting values from momentuHMM or miSum object returned by fitHMM or MI_summary
#'
#' @param m A \code{momentuHMM} or \code{miSum} object.
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
  
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MI_summary)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  
  Par <- list()
  if(is.miSum(m)){
    m$mle<-lapply(m$Par$real[distnames],function(x) x$est)
    for(i in distnames){
      if(!is.null(DM[[i]]))
        m$mle[[i]]<-m$Par$beta[[i]]$est
      else if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]])
        m$mle[[i]]<-m$mle[[i]][-1,]
      Par[[i]] <- c(t(unname(m$mle[[i]])))
    }
    beta <- unname(m$Par$beta$beta$est)
    delta <- unname(m$Par$real$delta$est)
  } else {
    for(i in distnames){
      if(is.null(DM[[i]])){
        par <- unname(m$mle[[i]])
        if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]]) 
          par <- par[-1,]
        par <- c(t(par))
      } else par <- unname(m$CI_beta[[i]]$est)
      Par[[i]] <- par
    }
    beta <- unname(m$mle$beta)
    delta <- unname(m$mle$delta)
  }
  
  list(Par=Par,beta=beta,delta=delta)
}