
#' Forward log-probabilities
#'
#' Used in \code{\link{stateProbs}} and \code{\link{pseudoRes}}.
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#'
#' @return The matrix of forward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' la <- momentuHMM:::logAlpha(m)
#' }

logAlpha <- function(m)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  if(is.miSum(m)){
    beta<-m$Par$beta$beta$est
    delta<-m$Par$real$delta$est
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
  }
  
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  nbObs <- nrow(m$data)
  lalpha <- matrix(NA,nbObs,nbStates)

  # identify covariates
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  
  covs <- model.matrix(newformula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  probs <- allProbs(m,nbStates)
  
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(m$data$ID==unique(m$data$ID)[i])[1])
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

  k <- 1
  for(i in 1:nbObs) {
    gamma <- trMat[,,i]
    if(any(i==aInd)){
      foo <- (delta[k,] %*% gamma)*probs[i,]
      lscale <- 0
      k <- k + 1
    } else {
      foo <- (foo %*% gamma)*probs[i,]
    }
    lscale <- lscale+log(sum(foo))
    foo <- foo/sum(foo)
    lalpha[i,] <- log(foo)+lscale
  }

  return(lalpha)
}
