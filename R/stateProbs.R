
#' State probabilities
#'
#' For a given model, computes the probability of the process being in the different states
#' at each time point.
#'
#' @param m A \code{momentuHMM} or \code{momentuHierHMM} object.
#' @param hierarchical Logical indicating whether or not to return a list of state probabilities for each level of a hierarchical HMM. Ignored unless \code{m} is a \code{\link{momentuHierHMM}} object.
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
#' @importFrom Brobdingnag sum brob
#' @export

stateProbs <- function(m, hierarchical=FALSE)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")
  
  m <- delta_bc(m)

  data <- m$data
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))

  if(nbStates==1)
    stop("No states to decode (nbStates=1)")

  nbObs <- nrow(data)
  la <- logAlpha(m) # forward log-probabilities
  lb <- logBeta(m) # backward log-probabilities
  stateProbs <- matrix(0,nbObs,nbStates)

  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,max(which(data$ID==unique(data$ID)[i])))
  
  mixtures <- m$conditions$mixtures
  if(mixtures>1) pie <- m$mle[["pi"]]
  else pie <- matrix(1,nbAnimals,1)
    
  eta <- mixtureProbs(m)
  
  k <- 1
  stateProb <- list()
  for(i in nbObs:1){
    for(mix in 1:mixtures){
      if(any(i==aInd)){
        c <- max(la[[mix]][i,]) # cancels out below ; prevents numerical errors
        llk <- c + log(sum(exp(la[[mix]][i,]-c)))
        ieta <- eta[k,]
        if(mix==mixtures) k <- k + 1
      }
      stateProb[[mix]] <- exp(la[[mix]][i,]+lb[[mix]][i,]-llk)
      if(any(!is.finite(stateProb[[mix]])) || !sum(stateProb[[mix]])) {
        stateProb[[mix]] <- Brobdingnag::brob(la[[mix]][i,]+lb[[mix]][i,]-llk)
        stateProbs[i,] <- stateProbs[i,] + as.numeric(stateProb[[mix]]/Brobdingnag::sum(stateProb[[mix]]) * ieta[mix])
      } else {
        stateProbs[i,] <- stateProbs[i,] + stateProb[[mix]] / sum(stateProb[[mix]]) * ieta[mix]
      }
    }
  }
  colnames(stateProbs) <- m$stateNames
  
  if(inherits(m,"momentuHierHMM") && hierarchical){
    return(hierStateProbs(m, stateProbs))
  } else return(stateProbs)
}

hierStateProbs <- function(m, stateProbs){
  
  installDataTree()
  
  out <- list()
  for(j in 1:(m$conditions$hierStates$height-1)){
    if(j==m$conditions$hierStates$height-1) ref <- m$conditions$hierStates$Get("state",filterFun=data.tree::isLeaf)
    else ref <- m$conditions$hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
    out[[paste0("level",j)]] <- stateProbs[which(m$data$level %in% c(j,paste0(j,"i"))),ref]
    colnames(out[[paste0("level",j)]]) <- names(ref)
  }
  class(out) <- append("hierarchical",class(out))
  out
}