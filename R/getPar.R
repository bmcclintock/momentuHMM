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
  
  if(!is.momentuHMM(m) & !is.momentuHMMMI(m))
    stop("'m' must be a momentuHMM or momentuHMMMI object (as output by fitHMM or MI_summary)")
  
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  
  if(is.momentuHMMMI(m)){
    m$mle<-lapply(m$Par[distnames],function(x) x$est)
  }
  
  Par <- list()
  for(i in distnames){
    if(is.null(DM[[i]])){
      par <- m$mle[[i]]
      if(dist[[i]] %in% angledists & !m$conditions$estAngleMean[[i]]) 
        par <- par[-1,]
      par <- c(t(par))
    } else par <- m$CI_beta[[i]]$est
    Par[[i]] <- par
  }
  Par
}