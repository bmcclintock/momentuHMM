getWorkBounds <- function(workBounds,distnames,wpar,parindex,parCount,DM,beta=NULL,delta=NULL){
  
  for(i in distnames){
    if(is.null(workBounds[[i]]) | is.null(DM[[i]])){
      workBounds[[i]] <- matrix(c(-Inf,Inf),parCount[[i]],2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds[[i]])) stop("workBounds$",i," must be a matrix")
      if(!all(dim(workBounds[[i]])==c(parCount[[i]],2))) stop("workBounds$",i," must be a 2 x ",parCount[[i]]," matrix")
    }
    tmpwpar<-wpar[parindex[[i]]+1:parCount[[i]]]
    
    ind1<-which(is.finite(workBounds[[i]][,1]) & is.infinite(workBounds[[i]][,2]))
    ind2<-which(is.finite(workBounds[[i]][,1]) & is.finite(workBounds[[i]][,2]))
    ind3<-which(is.infinite(workBounds[[i]][,1]) & is.finite(workBounds[[i]][,2]))
    
    tmpwpar[ind1] <- exp(tmpwpar[ind1])+workBounds[[i]][ind1,1]
    tmpwpar[ind2] <- (workBounds[[i]][ind2,2]-workBounds[[i]][ind2,1]) * boot::inv.logit(tmpwpar[ind2])+workBounds[[i]][ind2,1]
    tmpwpar[ind3] <- -(exp(-tmpwpar[ind3]) - workBounds[[i]][ind3,2])
    if(any(tmpwpar<workBounds[[i]][,1] | tmpwpar>workBounds[[i]][,2])) stop("starting values for ",i," data stream do not satisfy workBounds")
  }
  
  if(!is.null(beta)){
    if(is.null(workBounds$beta)){
      workBounds$beta <- matrix(c(-Inf,Inf),length(beta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$beta)) stop("workBounds$beta must be a matrix")
      if(!all(dim(workBounds$beta)==c(length(beta),2))) stop("workBounds$",i," must be a 2 x ",length(beta)," matrix")
    }
    ind1<-which(is.finite(workBounds$beta[,1]) & is.infinite(workBounds$beta[,2]))
    ind2<-which(is.finite(workBounds$beta[,1]) & is.finite(workBounds$beta[,2]))
    ind3<-which(is.infinite(workBounds$beta[,1]) & is.finite(workBounds$beta[,2]))
    
    beta[ind1] <- exp(beta[ind1])+workBounds$beta[ind1,1]
    beta[ind2] <- (workBounds$beta[ind2,2]-workBounds$beta[ind2,1]) * boot::inv.logit(beta[ind2])+workBounds$beta[ind2,1]
    beta[ind3] <- -(exp(-beta[ind3]) - workBounds$beta[ind3,2])
    if(any(c(beta)<workBounds$beta[,1] | c(beta)>workBounds$beta[,2])) stop("beta0 do not satisfy workBounds")
  }

  if(!is.null(delta)){
    if(is.null(workBounds$delta)){
      workBounds$delta <- matrix(c(-Inf,Inf),length(delta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$delta)) stop("workBounds$delta must be a matrix")
      if(!all(dim(workBounds$delta)==c(length(delta),2))) stop("workBounds$",i," must be a 2 x ",length(delta)," matrix")
    }
    ind1<-which(is.finite(workBounds$delta[,1]) & is.infinite(workBounds$delta[,2]))
    ind2<-which(is.finite(workBounds$delta[,1]) & is.finite(workBounds$delta[,2]))
    ind3<-which(is.infinite(workBounds$delta[,1]) & is.finite(workBounds$delta[,2]))
    
    delta[ind1] <- exp(delta[ind1])+workBounds$delta[ind1,1]
    delta[ind2] <- (workBounds$delta[ind2,2]-workBounds$delta[ind2,1]) * boot::inv.logit(delta[ind2])+workBounds$delta[ind2,1]
    delta[ind3] <- -(exp(-delta[ind3]) - workBounds$delta[ind3,2])
    if(any(c(delta)<workBounds$delta[,1] | c(delta)>workBounds$delta[,2])) stop("beta0 do not satisfy workBounds")
  }
  workBounds
}