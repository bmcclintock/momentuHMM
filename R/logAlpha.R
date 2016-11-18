
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
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  lalpha <- matrix(NA,nbObs,nbStates)

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="omega" & names(data)!="dry" & names(data)!="dive" & names(data)!="ice" & names(data)!="land")
  covs <- model.matrix(m$conditions$formula,data)

  allProbs <- allProbs(data,nbStates,m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,m$mle$stepPar,m$mle$anglePar,m$mle$omegaPar,m$mle$dryPar,m$mle$divePar,m$mle$icePar,m$mle$landPar,m$conditions$zeroInflation)

  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  lscale <- 0
  foo <- m$mle$delta*allProbs[1,]
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
