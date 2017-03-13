
#' Backward log-probabilities
#'
#' Used in \code{\link{stateProbs}}.
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return The matrix of backward log-probabilities.
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
  
  if(is.miSum(m)){
    beta<-m$Par$beta$beta$est
  } else {
    beta <- m$mle$beta
  }
  
  nbStates <- length(m$stateNames)
  nbObs <- nrow(m$data)
  lbeta <- matrix(NA,nbObs,nbStates)

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  
  allProbs <- allProbs(m,nbStates)
  
  trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))

  lscale <- log(nbStates)
  foo <- rep(1,nbStates)/nbStates
  lbeta[nbObs,] <- rep(0,nbStates)

  for(i in (nbObs-1):1) {
    gamma <- trMat[,,(i+1)]
    foo <- gamma%*%(allProbs[i+1,]*foo)
    lbeta[i,] <- log(foo)+lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }

  return(lbeta)
}
