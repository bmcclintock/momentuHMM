
#' Backward log-probabilities
#'
#' Used in \code{\link{stateProbs}}.
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return The matrix of backward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' lb <- logBeta(m)
#' }

logBeta <- function(m)
{
  nbStates <- length(m$stateNames)
  nbObs <- nrow(m$data)
  lbeta <- matrix(NA,nbObs,nbStates)

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  
  allProbs <- allProbs(m,nbStates)
  
  trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))

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
