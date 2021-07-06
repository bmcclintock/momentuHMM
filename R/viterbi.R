
#' Viterbi algorithm
#'
#' For a given model, reconstructs the most probable states sequence,
#' using the Viterbi algorithm.
#'
#' @param m An object \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}}
#' @param hierarchical Logical indicating whether or not to return a list of Viterbi-decoded states for each level of a hierarchical HMM. Ignored unless \code{m} is a \code{\link{momentuHierHMM}} object.
#'
#' @return The sequence of most probable states. If \code{hierarchical} is \code{TRUE}, then a list of the most probable states for each level of the hierarchy is returned.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' # reconstruction of states sequence
#' states <- viterbi(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export

viterbi <- function(m, hierarchical=FALSE)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")
  
  m <- delta_bc(m)
  
  data <- m$data
  nbStates <- length(m$stateNames)
  nbAnimals <- length(unique(m$data$ID))
  beta <- m$mle$beta
  delta <- m$mle$delta
  
  if(nbStates==1)
    stop("No states to decode (nbStates=1)")
  
  # identify covariates
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,m$data,par=m$mle)
  covs <- reForm$covs
  aInd <- reForm$aInd
  
  probs <- allProbs(m)
  mixtures <- m$conditions$mixtures
  mixProbs <- mixtureProbs(m)
  
  trMat <- list()
  for(mix in 1:mixtures){
    trMat[[mix]] <- trMatrix_rcpp(nbStates,beta[(mix-1)*ncol(covs)+1:ncol(covs),,drop=FALSE],as.matrix(covs),m$conditions$betaRef)
  }
  
  allStates <- NULL
  
  for(zoo in 1:nbAnimals) {
    
    nbObs <- length(which(data$ID==unique(data$ID)[zoo])) # nb of observations for animal zoo
    tmxi <- matrix(0,nbObs,nbStates)
    xi_mix <- matrix(0,nbObs,nbStates)
    
    if(zoo!=nbAnimals) {
      tm <- lapply(trMat,function(x) x[,,aInd[zoo]:(aInd[zoo+1]-1)])
      p <- probs[aInd[zoo]:(aInd[zoo+1]-1),]
    } else {
      tm <- lapply(trMat,function(x) x[,,aInd[zoo]:nrow(probs)])
      p <- probs[aInd[zoo]:nrow(probs),]
    }
    foo <- rowSums(mapply(function(mix) mixProbs[zoo,mix]*delta[(mix-1)*nbAnimals+zoo,]%*%tm[[mix]][,,1]*p[1,],1:mixtures))
    xi_mix[1,] <- foo/sum(foo)
    
    for(i in 2:nbObs) {
      foo <- Reduce("+",mapply(function(mix) mixProbs[zoo,mix]*xi_mix[i-1,]*tm[[mix]][,,i],1:mixtures,SIMPLIFY = FALSE))
      foo <- apply(foo,2,max)*p[i,]
      xi_mix[i,] <- foo/sum(foo)
    }
    
    stSeq <- rep(NA,nbObs)
    tmxi[nbObs,] <- xi_mix[nbObs,]
    stSeq[nbObs] <- which.max(tmxi[nbObs,])
    for(i in (nbObs-1):1){
      for(mix in 1:mixtures){
        tmxi[i,] <- tmxi[i,] + xi_mix[i,]*mixProbs[zoo,mix]*tm[[mix]][,stSeq[i+1],i+1]
      }
      stSeq[i] <- which.max(tmxi[i,])
    }
    allStates <- c(allStates,stSeq)
  }
  
  if(inherits(m,"momentuHierHMM") && hierarchical){
    installDataTree()
    return(hierViterbi(m, allStates))
  } else return(allStates)
}

hierViterbi <- function(m, allStates, stateNames = TRUE){
    out <- list()
    for(j in 1:(m$conditions$hierStates$height-1)){
      if(j==m$conditions$hierStates$height-1) ref <- m$conditions$hierStates$Get("state",filterFun=data.tree::isLeaf)
      else ref <- m$conditions$hierStates$Get(function(x) data.tree::Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)
      out[[paste0("level",j)]] <- match(allStates[which(m$data$level %in% c(j,paste0(j,"i")))],ref,nomatch=0)
      if(stateNames) out[[paste0("level",j)]] <- names(ref)[out[[paste0("level",j)]]]
    }
    class(out) <- append("hierarchical",class(out))
    out
}
