
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
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.momentuHMMMI(m))
    stop("'m' must be a momentuHMM, miHMM, or momentuHMMMI object (as output by fitHMM, MIfitHMM, or MI_summary)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  if(is.momentuHMMMI(m)){
    beta<-m$Par$beta$beta$est
    delta<-m$Par$real$delta$est
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
  }
  
  nbStates <- length(m$stateNames)
  nbObs <- nrow(m$data)
  lalpha <- matrix(NA,nbObs,nbStates)

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  allProbs <- allProbs(m,nbStates)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  lscale <- 0
  foo <- (delta%*%trMat[,,1])*allProbs[1,]
  lalpha[1,] <- log(foo)+lscale

  for(i in 2:nbObs) {
    gamma <- trMat[,,i]
    foo <- foo%*%gamma*allProbs[i,]
    lscale <- lscale+log(sum(foo))
    foo <- foo/sum(foo)
    lalpha[i,] <- log(foo)+lscale
  }

  return(lalpha)
}
