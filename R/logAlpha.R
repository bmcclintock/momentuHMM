
#' Forward log-probabilities
#'
#' Used in \code{\link{stateProbs}} and \code{\link{pseudoRes}}.
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return The matrix of forward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' la <- logAlpha(m)
#' }

logAlpha <- function(m)
{
  nbStates <- length(m$stateNames)
  nbObs <- nrow(m$data)
  lalpha <- matrix(NA,nbObs,nbStates)

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  allProbs <- allProbs(m,nbStates)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  lscale <- 0
  foo <- m$mle$delta*allProbs[1,]
  lalpha[1,] <- log(foo)+lscale

  for(i in 2:nbObs) {
    gamma <- trMat[,,i-1]
    foo <- foo%*%gamma*allProbs[i,]
    lscale <- lscale+log(sum(foo))
    foo <- foo/sum(foo)
    lalpha[i,] <- log(foo)+lscale
  }

  return(lalpha)
}
