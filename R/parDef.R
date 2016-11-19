
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

parDef <- function(dist,nbStates,estAngleMean,zeroInflation,userBounds=NULL,DM)
{
  distnames<-names(dist)
  eval(parse(text=paste0(distnames,"Dist='",dist,"'")))
  parSize <- rep(0,length(distnames))
  otherDist <- distnames[-which(distnames %in% c("step","angle"))]
  parNames <- vector('list',length(dist))
  names(parNames) <- distnames

  switch(stepDist,
         "gamma"={
           parSize[1] <- 2
           parNames$step <- c("mean","sd")
           stepBounds <- matrix(rep(c(0,Inf),ncol(DM$step)),ncol=2,byrow=TRUE)
           rownames(stepBounds) <-paste0("stepbeta_",1:ncol(DM$step))
             
           # include zero-mass
           if(zeroInflation) {
            parSize[1] <- parSize[1]+1
            stepBounds[ncol(DM$step)-nbStates:1+1,2] <- 1
            parNames$step <- c(parNames,"zero-mass")
           } 
         },
         "weibull"={
           parSize[1] <- 2
           parNames$step <- c("shape","scale")
           stepBounds <- matrix(rep(c(0,Inf),ncol(DM$step)),ncol=2,byrow=TRUE)
           rownames(stepBounds) <-paste0("stepbeta_",1:ncol(DM$step))
           
           # include zero-mass
           if(zeroInflation) {
             parSize[1] <- parSize[1]+1
             stepBounds[ncol(DM$step)-nbStates:1+1,2] <- 1
             parNames$step <- c(parNames,"zero-mass")
           } 
         },
         "lnorm"={
           parSize[1] <- 2
           stepBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                ncol=2,byrow=TRUE)
           parNames$step <- c("location","scale")
         },
         "exp"={
           parSize[1] <- 1
           stepBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
           parNames$step <- c("rate")
         },
         "none"={
           parSize[1]<-0
           stepBounds <- NULL
           #parNames$step <- NULL
         })


  switch(angleDist,
         "none"={
           parSize[2] <- 0
           angleBounds <- NULL
           #parNames$angle <- NULL
         },
         "vm"={
           if(estAngleMean) { # if the angle mean is estimated
             parSize[2] <- 2
             # bounds are chosen such that the parameters are not scaled
             # (already in the right intervals for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,Inf),nbStates)),
                                   ncol=2,byrow=TRUE)
             rownames(angleBounds)<-c(paste0(rep("mean",each=nbStates),1:nbStates),paste0("anglebeta_",1:ncol(DM$angle)))
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
             rownames(angleBounds)<-paste0("anglebeta_",1:ncol(DM$angle))
             parNames$angle <- c("mean","sd")
           }
         },
         "wrpcauchy"={
           if(estAngleMean) {
             parSize[2] <- 2
             # bounds are chosen such that the mean is not scaled, but the concentration is
             # scaled from ]0,1[ to ]0,Inf[ (for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,1),nbStates)),
                                   ncol=2,byrow=TRUE)
             rownames(angleBounds)<-c(paste0(rep("mean",each=nbStates),1:nbStates),paste0("anglebeta_",1:ncol(DM$angle)))
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0, 1), nbStates), ncol = 2,byrow = TRUE)
             rownames(angleBounds)<-paste0("anglebeta_",1:ncol(DM$angle))
             parNames$angle <- c("mean","concentration")
           }
         })
  
  bounds <- list(step=stepBounds,angle=angleBounds)
  
  if(length(dist[otherDist])){
    for(i in 1:length(dist[otherDist])){
      Dist<-dist[otherDist][[i]]
      dataname<-names(dist[otherDist][i])
      dm<-DM[otherDist][[i]]
      switch(Dist,
             "beta"={
               parSize[2+i] <- 2
               tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
               parNames[[dataname]]<-c("shape1","shape2")
             },
             "pois"={
               parSize[2+i] <- 1
               tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
               parNames[[dataname]]<-"lambda"
             },
             "weibull"={
               parSize[2+i] <- 2
               tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
               parNames[[dataname]]<-c("shape","scale")
             },
             "gamma"={
               parSize[2+i] <- 2
               tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
               parNames[[dataname]]<-c("mean","sd")
             })
      rownames(tmpbounds)<-paste0(names(Dist),"beta_",1:ncol(dm))
      bounds[[dataname]] <- tmpbounds
    }
  }
  
  if(!is.null(userBounds)) {
    for(i in names(dist)){
      rownames(userBounds[[i]])<-rownames(bounds[[i]])
    }
    bounds <- userBounds
  }
  boundInd <- lapply(DM,getboundInd)
  return(list(parSize=parSize,bounds=bounds[distnames],parNames=parNames,boundInd=boundInd[distnames]))
}
