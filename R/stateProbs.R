
#' State probabilities
#'
#' For a given model, computes the probability of the process being in the different states
#' at each time point.
#'
#' @param m A \code{momentuHMM} object.
#'
#' @return The matrix of state probabilities, with element [i,j] the probability
#' of being in state j in observation i.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' sp <- stateProbs(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export

stateProbs <- function(m)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  data <- m$data
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))

  if(nbStates==1)
    stop("No states to decode (nbStates=1)")

  nbObs <- nrow(data)
  la <- logAlpha(m) # forward log-probabilities
  lb <- logBeta(m) # backward log-probabilities
  stateProbs <- matrix(NA,nbObs,nbStates)

  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,max(which(data$ID==unique(data$ID)[i])))
  
  mixtures <- m$conditions$mixtures
  if(mixtures>1) pie <- m$mle$pi
  else pie <- matrix(1,nbAnimals,1)
    
  # get probability that individual is in each mixture
  eta <- lnum <- matrix(0,nbAnimals,mixtures)
  for(i in 1:nbAnimals){
    for(mix in 1:mixtures){
      c <- max(la[[mix]][aInd[i],]+log(pie[i,mix]))
      lnum[i,mix] <- c + log(sum(exp(la[[mix]][aInd[i],]+log(pie[i,mix])-c)))
    }
    c <- max(lnum[i,])
    eta[i,] <- exp(lnum[i,] - c - log(sum(exp(lnum[i,]-c))))
  }
  
  k <- 1
  stateProb <- matrix(NA,mixtures,nbStates)
  for(i in nbObs:1){
    for(mix in 1:mixtures){
      if(any(i==aInd)){
        c <- max(la[[mix]][i,]) # cancels out below ; prevents numerical errors
        llk <- c + log(sum(exp(la[[mix]][i,]-c)))
        ieta <- eta[k,]
        if(mix==mixtures) k <- k + 1
      }
      stateProb[mix,] <- exp(la[[mix]][i,]+lb[[mix]][i,]-llk)
    }
    stateProbs[i,] <- colSums(stateProb/rowSums(stateProb) * ieta)
  }
  colnames(stateProbs) <- m$stateNames
  
  return(stateProbs)
}
