
#' Matrix of all probabilities
#'
#' Used in functions \code{\link{viterbi}}, \code{\link{logAlpha}}, \code{\link{logBeta}}.
#'
#' @param data Object \code{momentuHMMData}.
#' @param nbStates Number of states of the HMM.
#' @param stepDist Name of the distribution of the step lengths.
#' @param angleDist Name of the distribution of the turning angles.
#' Set to "none" if the angle distribution should not be estimated.
#' @param stepPar Parameters of the step length distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of \code{stepDist}),
#' and one column for each state.
#' @param anglePar Parameters of the turning angle distribution. Must be provided in a
#' matrix with one row for each parameter (in the order expected by the pdf of \code{angleDist}),
#' and one column for each state. Default: \code{NULL} ; if the turning angles distribution
#' is not estimated.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default: \code{FALSE}.
#'
#' @return Matrix of all probabilities.
#'
#' @examples
#' \dontrun{
#' stepPar <- c(1,10,1,5,0.2,0.3)
#' anglePar <- c(0,pi,0.5,2)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                  anglePar=anglePar,nbCovs=2,zeroInflation=TRUE)
#' P <- allProbs(data=data,nbStates=2,stepDist=stepDist,angleDist=angleDist,
#'                stepPar=matrix(stepPar,ncol=2,byrow=TRUE),anglePar=matrix(anglePar,ncol=2,
#'                byrow=TRUE),zeroInflation=TRUE)
#' }

allProbs <- function(data,nbStates,dist,par,zeroInflation)
{
  
  Fun <- lapply(dist,function(x) paste("d",x,sep=""))
  
  nbObs <- nrow(data)
  allProbs <- matrix(1,nrow=nbObs,ncol=nbStates)
  
  for(i in names(dist)){
  
    genInd <- which(!is.na(data[[i]]))
    sp <- par[[i]]
  
    for(state in 1:nbStates) {
      genPar <- sp
      genProb <- rep(1,nbObs)
      genFun <- Fun[[i]]
      
      # Constitute the lists of state-dependent parameters for the step and angle
      genArgs <- list(data[[i]][genInd])
      if(zeroInflation[[i]]) {
        zeromass <- genPar[nrow(genPar),state]
        genPar <- genPar[-nrow(genPar),]
      }
  
      for(j in 1:nrow(genPar))
        genArgs[[j+1]] <- genPar[j,state]
  
      # conversion between mean/sd and shape/scale if necessary
      if(dist[[i]]=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
      }
      if(zeroInflation[[i]]) {
        genProb[genInd] <- ifelse(data[[i]][genInd]==0,
                                    zeromass, # if gen==0
                                    (1-zeromass)*do.call(genFun,genArgs)) # if gen != 0
      }
      else genProb[genInd] <- do.call(genFun,genArgs)
  
      allProbs[,state] <- allProbs[,state]*genProb;
    }
  }
  return(allProbs)
}
