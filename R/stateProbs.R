
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
    
  for(i in nbObs:1){
    if(any(i==aInd)){
      c <- max(la[i,]) # cancels out below ; prevents numerical errors
      llk <- c + log(sum(exp(la[i,]-c)))
    }
    stateProbs[i,] <- exp(la[i,]+lb[i,]-llk)
  }
  colnames(stateProbs) <- m$stateNames
  
  return(stateProbs)
}
