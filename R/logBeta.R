
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
    delta<-m$Par$real$delta$est
  } else {
    beta <- m$mle$beta
    delta <- m$mle$delta
  }
  
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  nbObs <- nrow(m$data)
  lbeta <- matrix(NA,nbObs,nbStates)

  # identify covariates
  formula<-m$conditions$formula
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  
  covs <- model.matrix(newformula,m$data)
  
  probs <- allProbs(m,nbStates)
  
  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))
  
  aInd <- NULL
  #aInd2 <- NULL
  for(i in 1:nbAnimals){
    aInd <- c(aInd,max(which(m$data$ID==unique(m$data$ID)[i])))
    #aInd2 <- c(aInd2,which(m$data$ID==unique(m$data$ID)[i])[1])
  }

  for(i in nbObs:1) {
    
    if(any(i==aInd)){
        foo <- rep(1,nbStates)
        lscale <- 0
    } else {
      gamma <- trMat[,,i+1]
      #if(any(i==aInd2)){
      #  gamma <- delta %*% gamma
      #}
      foo <- gamma%*%(probs[i+1,]*foo)
    }
    lbeta[i,] <- log(foo)+lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }

  return(lbeta)
}
