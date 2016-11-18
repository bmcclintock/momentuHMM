
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
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  lbeta <- matrix(NA,nbObs,nbStates)

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="omega" & names(data)!="dry" & names(data)!="dive" & names(data)!="ice" & names(data)!="land")
  covs <- model.matrix(m$conditions$formula,data)
  
  allProbs <- allProbs(data,nbStates,m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,m$mle$stepPar,m$mle$anglePar,m$mle$omegaPar,m$mle$dryPar,m$mle$divePar,m$mle$icePar,m$mle$landPar,m$conditions$zeroInflation)
  
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
