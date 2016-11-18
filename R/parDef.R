
#' Parameters definition
#'
#' @param stepDist Name of the distribution of the step lengths.
#' @param angleDist Name of the distribution of the turning angles.
#' Set to \code{"none"} if the angle distribution should not be estimated.
#' @param nbStates Number of states of the HMM.
#' @param estAngleMean \code{TRUE} if the mean of the turning angles distribution is estimated,
#' \code{FALSE} otherwise.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#'
#' @return A list of:
#' \item{parSize}{Vector of two values: number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution}
#' \item{bounds}{Matrix with 2 columns and \code{sum(parSize)} rows - each row contains the lower and upper
#' bound for the correponding parameter)}
#' \item{parNames}{Names of parameters of step distribution (the names of the parameters of the
#' angle distribution are always the same).}

parDef <- function(stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,nbStates,estAngleMean,zeroInflation,userBounds=NULL,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM)
{
  parSize <- rep(0,7)

  switch(stepDist,
         "gamma"={
           parSize[1] <- 2
           #stepBounds <- matrix(rep(c(0,Inf),2*nbStates),ncol=2,byrow=TRUE)
           parNames <- c("mean","sd")
           stepBounds <- matrix(rep(c(0,Inf),ncol(stepDM)),ncol=2,byrow=TRUE)
           boundNames <-paste0("beta_",1:ncol(stepDM))
             
           # include zero-mass
           if(zeroInflation) {
            parSize[1] <- parSize[1]+1
            stepBounds[ncol(stepDM)-nbStates:1+1,2] <- 1
            parNames <- c(parNames,"zero-mass")
            #boundNames <- c(boundNames,paste0("zero-mass",1:nbStates))
           } 
  
  rownames(stepBounds)<-paste0(stepDist,boundNames)
         },
         "weibull"={
           parSize[1] <- 2
           stepBounds <- matrix(rep(c(0,Inf),2*nbStates),ncol=2,byrow=TRUE)
           parNames <- c("shape","scale")
         },
         "lnorm"={
           parSize[1] <- 2
           stepBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                ncol=2,byrow=TRUE)
           parNames <- c("location","scale")
         },
         "exp"={
           parSize[1] <- 1
           stepBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
           parNames <- c("rate")
         })


  switch(angleDist,
         "none"={
           parSize[2] <- 0
           angleBounds <- NULL
         },
         "vm"={
           if(estAngleMean) { # if the angle mean is estimated
             parSize[2] <- 2
             # bounds are chosen such that the parameters are not scaled
             # (already in the right intervals for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,Inf),nbStates)),
                                   ncol=2,byrow=TRUE)
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
           }
         },
         "wrpcauchy"={
           if(estAngleMean) {
             parSize[2] <- 2
             # bounds are chosen such that the mean is not scaled, but the concentration is
             # scaled from ]0,1[ to ]0,Inf[ (for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,1),nbStates)),
                                   ncol=2,byrow=TRUE)
             rownames(angleBounds)<-c(paste0(rep("mean",each=nbStates),1:nbStates),paste0("concbeta_",1:ncol(angleDM)))
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0, 1), nbStates), ncol = 2,byrow = TRUE)
             rownames(angleBounds)<-paste0("concbeta_",1:ncol(angleDM))
           }
         })
  
  switch(omegaDist,
          "none"={
           parSize[3] <- 0
           omegaBounds <- NULL
         },
         "beta"={
           parSize[3] <- 2
           omegaBounds <- matrix(rep(c(0,Inf),ncol(omegaDM)),ncol=2,byrow=TRUE)
           rownames(omegaBounds)<-paste0("omegabeta_",1:ncol(omegaDM))
         })
  
  switch(dryDist,
           "none"={
           parSize[4] <- 0
           dryBounds <- NULL
         },
         "beta"={
           parSize[4] <- 2
           dryBounds <- matrix(rep(c(0,Inf),ncol(dryDM)),ncol=2,byrow=TRUE)
           rownames(dryBounds)<-paste0("drybeta_",1:ncol(dryDM))
         })

  switch(diveDist,
          "none"={
           parSize[5] <- 0
           diveBounds <- NULL
         },
         "pois"={
           parSize[5] <- 1
           diveBounds <- matrix(rep(c(0,Inf),ncol(diveDM)),ncol=2,byrow=TRUE)
           rownames(diveBounds)<-paste0("divebeta_",1:ncol(diveDM))
         })
  
  switch(iceDist,
          "none"={
           parSize[6] <- 0
           iceBounds <- NULL
         },
         "beta"={
           parSize[6] <- 2
           iceBounds <- matrix(rep(c(0,Inf),ncol(iceDM)),ncol=2,byrow=TRUE)
           rownames(iceBounds)<-paste0("icebeta_",1:ncol(iceDM))
         })
  
  switch(landDist,
          "none"={
           parSize[7] <- 0
           landBounds <- NULL
         },         
         "beta"={
           parSize[7] <- 2
           landBounds <- matrix(rep(c(0,Inf),ncol(landDM)),ncol=2,byrow=TRUE)
           rownames(landBounds)<-paste0("landbeta_",1:ncol(landDM))
         })
  
  bounds <- rbind(stepBounds,angleBounds,omegaBounds,dryBounds,diveBounds,iceBounds,landBounds)
  if(!is.null(userBounds)) {
    #print(cbind(bounds,userBounds))
    rownames(userBounds)<-rownames(bounds)
    bounds <- userBounds
  }
  boundInd<-list(step=getboundInd(stepDM),angle=getboundInd(angleDM),omega=getboundInd(omegaDM),dry=getboundInd(dryDM),dive=getboundInd(diveDM),ice=getboundInd(iceDM),land=getboundInd(landDM))
  return(list(parSize=parSize,bounds=bounds,parNames=parNames,boundInd=boundInd))
}
