
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
#' @importFrom Brobdingnag brob
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
  probs <- allProbs(m)
  la <- logAlpha(m, probs=probs) # forward log-probabilities
  lb <- logBeta(m, probs=probs) # backward log-probabilities

  if (!is.list(la)) la <- list(la)
  if (!is.list(lb)) lb <- list(lb)

  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,max(which(data$ID==unique(data$ID)[i])))
  
  mixtures <- m$conditions$mixtures
  if(mixtures>1) pie <- m$mle[["pi"]]
  else pie <- matrix(1,nbAnimals,1)
    
  eta <- mixtureProbs(m)
  
  stateProbs_mat <- stateProbs_cpp(
    nbObs = nbObs,
    nbStates = nbStates,
    mixtures = mixtures,
    aInd = aInd,
    eta = eta,
    laList = la,
    lbList = lb
  )
  colnames(stateProbs_mat) <- m$stateNames
  
  if(inherits(m,"momentuHierHMM") && hierarchical){
    return(hierStateProbs(m, stateProbs_mat))
  } else return(stateProbs_mat)
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