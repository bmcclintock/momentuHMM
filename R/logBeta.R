
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
  data <- m$data
  nbStates <- length(m$stateNames)
  nbObs <- nrow(data)
  lbeta <- matrix(NA,nbObs,nbStates)

  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM
  par <- m$mle

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  allProbs <- allProbs(data,nbStates,dist,par,m$conditions$zeroInflation)
  
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
