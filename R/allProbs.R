
#' Matrix of all probabilities
#'
#' Used in functions \code{\link{viterbi}}, \code{\link{logAlpha}}, \code{\link{logBeta}}.
#'
#' @param data Object \code{moveData}.
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

allProbs <- function(data,nbStates,stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,stepPar,anglePar=NULL,omegaPar,dryPar,divePar,icePar,landPar,zeroInflation=FALSE)
{
  stepFun <- paste("d",stepDist,sep="")
  if(angleDist!="none") angleFun <- paste("d",angleDist,sep="")
  omegaFun <- paste("d",omegaDist,sep="")
  dryFun <- paste("d",dryDist,sep="")
  diveFun <- paste("d",diveDist,sep="")
  iceFun <- paste("d",iceDist,sep="")
  landFun <- paste("d",landDist,sep="")
  
  nbObs <- length(data$step)
  allProbs <- matrix(1,nrow=nbObs,ncol=nbStates)
  stepInd <- which(!is.na(data$step))
  if(angleDist!="none") angleInd <- which(!is.na(data$angle))
  omegaInd <- which(!is.na(data$omega))
  dryInd <- which(!is.na(data$dry))
  diveInd <- which(!is.na(data$dive))
  iceInd <- which(!is.na(data$ice))
  landInd <- which(!is.na(data$land))
  
  sp <- stepPar

  for(state in 1:nbStates) {
    stepPar <- sp
    stepProb <- rep(1,nbObs)
    angleProb <- rep(1,nbObs)
    omegaProb <- rep(1,nbObs)
    dryProb <- rep(1,nbObs)
    diveProb <- rep(1,nbObs)
    iceProb <- rep(1,nbObs)
    landProb <- rep(1,nbObs)
    
    # Constitute the lists of state-dependent parameters for the step and angle
    stepArgs <- list(data$step[stepInd])
    if(angleDist!="none") angleArgs <- list(data$angle[angleInd])
    omegaArgs <- list(data$omega[omegaInd])
    dryArgs <- list(data$dry[dryInd])
    diveArgs <- list(data$dive[diveInd])
    iceArgs <- list(data$ice[iceInd])
    landArgs <- list(data$land[landInd])
    if(zeroInflation) {
      zeromass <- stepPar[nrow(stepPar),state]
      stepPar <- stepPar[-nrow(stepPar),]
    }

    for(j in 1:nrow(stepPar))
      stepArgs[[j+1]] <- stepPar[j,state]

    # conversion between mean/sd and shape/scale if necessary
    if(stepDist=="gamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      scale <- stepArgs[[3]]^2/stepArgs[[2]]
      stepArgs[[2]] <- shape
      stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
    }
    if(zeroInflation) {
      stepProb[stepInd] <- ifelse(data$step[stepInd]==0,
                                  zeromass, # if step==0
                                  (1-zeromass)*do.call(stepFun,stepArgs)) # if step != 0
    }
    else stepProb[stepInd] <- do.call(stepFun,stepArgs)

    if(angleDist!="none") {
      for(j in 1:nrow(anglePar))
        angleArgs[[j+1]] <- anglePar[j,state]

      angleProb[angleInd] <- do.call(angleFun,angleArgs)
      allProbs[,state] <- stepProb*angleProb
    }
    else allProbs[,state] <- stepProb # model step length only

    if(omegaDist!="none") {
      for(j in 1:nrow(omegaPar))
        omegaArgs[[j+1]] <- omegaPar[j,state]
  
      omegaProb[omegaInd] <- do.call(omegaFun,omegaArgs)
      allProbs[,state] <- allProbs[,state]*omegaProb;
    }
    if(dryDist!="none") {
      for(j in 1:nrow(dryPar))
        dryArgs[[j+1]] <- dryPar[j,state]
  
      dryProb[dryInd] <- do.call(dryFun,dryArgs)
      allProbs[,state] <- allProbs[,state]*dryProb;    
    }
    
    if(diveDist!="none") {
      for(j in 1:nrow(divePar))
        diveArgs[[j+1]] <- divePar[j,state]
  
      diveProb[diveInd] <- do.call(diveFun,diveArgs)
      allProbs[,state] <- allProbs[,state]*diveProb;
    }
    
    if(iceDist!="none") {
      for(j in 1:nrow(icePar))
        iceArgs[[j+1]] <- icePar[j,state]
  
      iceProb[iceInd] <- do.call(iceFun,iceArgs)
      allProbs[,state] <- allProbs[,state]*iceProb;
    }
    
    if(landDist!="none") {
      for(j in 1:nrow(landPar))
        landArgs[[j+1]] <- landPar[j,state]
  
      landProb[landInd] <- do.call(landFun,landArgs)
      allProbs[,state] <- allProbs[,state]*landProb;
    }  
  }
  return(allProbs)
}
