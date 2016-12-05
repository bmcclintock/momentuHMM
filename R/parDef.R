
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
  parSize <- rep(0,length(distnames))
  parNames <- vector('list',length(dist))
  names(parNames) <- distnames
  bounds <- vector('list',length(dist))
  names(bounds) <- distnames

  for(i in 1:length(dist)){
    Dist<-dist[[i]]
    dataname<-names(dist[i])
    dm<-DM[[i]]
    switch(Dist,
           "beta"={
             parSize[i] <- 2
             tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
             parNames[[dataname]]<-c("shape1","shape2")
           },
           "pois"={
             parSize[i] <- 1
             tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
             parNames[[dataname]]<-"lambda"
           },
           "weibull"={
             parSize[i] <- 2
             tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
             parNames[[dataname]]<-c("shape","scale")
           },
           "gamma"={
             parSize[i] <- 2
             tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
             parNames[[dataname]]<-c("mean","sd")
           },
           "lnorm"={
             parSize[i] <- 2
             tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                  ncol=2,byrow=TRUE)
             parNames[[dataname]] <- c("location","scale")
           },
           "exp"={
             parSize[i] <- 1
             tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
             parNames[[dataname]] <- c("rate")
           },
           "vm"={
             if(estAngleMean[[i]]) { # if the angle mean is estimated
               parSize[i] <- 2
               # bounds are chosen such that the parameters are not scaled
               # (already in the right intervals for computing x and y)
               #tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,Inf),nbStates)),
               #                       ncol=2,byrow=TRUE)
               meanind<-unique(unlist(apply(dm[1:nbStates,],1,function(x) which(x!=0))))
               sdind<-unique(unlist(apply(dm[nbStates+1:nbStates,],1,function(x) which(x!=0))))
               if(any(intersect(meanind,sdind))) stop("'DM' for ",distname," cannot have common parameters in common for mean and sd")
               tmpbounds <- rbind(matrix(rep(c(-pi,pi),length(meanind)),ncol=2,byrow=TRUE),matrix(rep(c(0,Inf),length(sdind)),ncol=2,byrow=TRUE))
               parNames[[dataname]] <- c("mean","sd") 
            }
             else {
               parSize[i] <- 1
               tmpbounds <- matrix(rep(c(0,Inf),ncol(dm)),ncol=2,byrow=TRUE)
               parNames[[dataname]] <- c("mean","sd")
             }
           },
           "wrpcauchy"={
             if(estAngleMean[[i]]) {
               parSize[i] <- 2
               # bounds are chosen such that the mean is not scaled, but the concentration is
               # scaled from ]0,1[ to ]0,Inf[ (for computing x and y)
               #tmpbounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,1),nbStates)),
               #                       ncol=2,byrow=TRUE)
               meanind<-unique(unlist(apply(dm[1:nbStates,],1,function(x) which(x!=0))))
               sdind<-unique(unlist(apply(dm[nbStates+1:nbStates,],1,function(x) which(x!=0))))
               if(any(intersect(meanind,sdind))) stop("'DM' for ",distname," cannot have common parameters in common for mean and concentration")
               tmpbounds <- rbind(matrix(rep(c(-pi,pi),length(meanind)),ncol=2,byrow=TRUE),matrix(rep(c(0,1),length(sdind)),ncol=2,byrow=TRUE))
               parNames[[dataname]] <- c("mean","concentration")
             }
             else {
               parSize[i] <- 1
               tmpbounds <- matrix(rep(c(0, 1), ncol(dm)), ncol = 2,byrow = TRUE)
               parNames[[dataname]] <- c("mean","concentration")
             }
           }
          )
    if(zeroInflation[[dataname]]) {
      parSize[i] <- parSize[i]+1
      tmpbounds[ncol(dm)-nbStates:1+1,2] <- 1
      parNames[[dataname]] <- c(parNames[[dataname]],"zero-mass")
    } 
    rownames(tmpbounds)<-paste0(dataname,"beta_",1:ncol(dm))
    bounds[[dataname]] <- tmpbounds
  }
  
  if(!is.null(userBounds)) {
    for(i in names(dist)){
      if(is.null(userBounds[[i]])) 
        userBounds[[i]]<-bounds[[i]]
      rownames(userBounds[[i]])<-rownames(bounds[[i]])
    }
    bounds <- userBounds
  }
  boundInd <- lapply(DM,getboundInd)
  return(list(parSize=parSize,bounds=bounds[distnames],parNames=parNames,boundInd=boundInd[distnames]))
}
