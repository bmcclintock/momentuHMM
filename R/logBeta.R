
#' Backward log-probabilities
#'
#' Used in \code{\link{stateProbs}}.
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return A list of length \code{model$conditions$mixtures} where each element is a matrix of backward log-probabilities for each mixture.
#'
#' @examples
#' \dontrun{
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' lb <- momentuHMM:::logBeta(m)
#' }

logBeta <- function(m)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  m <- delta_bc(m)
  
  if(is.miSum(m)){
    beta<-m$Par$beta$beta$est
    delta<-m$Par$real$delta$est
    g0<-m$Par$beta$g0$est
    theta<-m$Par$beta$theta$est
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
    g0 <- m$mle$g0
    theta <- m$mle$theta
  }
  
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  nbObs <- nrow(m$data)
  lbeta <- matrix(NA,nbObs,nbStates)
  
  # identify covariates
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,m$data,par=m$mle)
  covs <- reForm$covs
  
  probs <- allProbs(m)
  
  mixtures <- m$conditions$mixtures
  
  trMat <- lbeta <- list()
  lbeta[1:mixtures] <- list(matrix(NA,nbObs,nbStates))
  
  for(mix in 1:mixtures){
    if(nbStates>1)
      trMat[[mix]] <- trMatrix_rcpp(nbStates,beta[(mix-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],as.matrix(covs),m$conditions$betaRef)
    else
      trMat[[mix]] <- array(1,dim=c(1,1,nbObs))
  }
  
  aInd <- NULL
  #aInd2 <- NULL
  for(i in 1:nbAnimals){
    aInd <- c(aInd,max(which(m$data$ID==unique(m$data$ID)[i])))
    #aInd2 <- c(aInd2,which(m$data$ID==unique(m$data$ID)[i])[1])
  }
  
  for(mix in 1:mixtures){
    for(i in nbObs:1) {
      
      if(any(i==aInd)){
        foo <- rep(1,nbStates)
        lscale <- 0
      } else {
        gamma <- trMat[[mix]][,,i+1]
        #if(any(i==aInd2)){
        #  gamma <- delta %*% gamma
        #}
        foo <- gamma%*%(probs[i+1,]*foo)
      }
      lbeta[[mix]][i,] <- log(foo)+lscale
      sumfoo <- sum(foo)
      foo <- foo/sumfoo
      lscale <- lscale+log(sumfoo)
    }
  }
  
  return(lbeta)
}
