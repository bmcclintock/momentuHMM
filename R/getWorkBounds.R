getWorkBounds <- function(workBounds,distnames,wpar,parindex,parCount,DM,beta0=NULL,delta0=NULL){
  
  lowerBounds <- upperBounds <- vector('list',length(distnames))
  for(i in distnames){
    if(is.null(workBounds[[i]]) | is.null(DM[[i]])){
      lowerBounds[[i]] <- rep(-Inf,parCount[[i]])
      upperBounds[[i]] <- rep(Inf,parCount[[i]])
    } else {
      if(!is.matrix(workBounds[[i]])) stop("workBounds$",i," must be a matrix")
      if(!all(dim(workBounds[[i]])==c(parCount[[i]],2))) stop("workBounds$",i," must be a 2 x ",parCount[[i]]," matrix")
      lowerBounds[[i]] <- workBounds[[i]][,1]
      upperBounds[[i]] <- workBounds[[i]][,2]
    }
    tmpwpar<-wpar[parindex[[i]]+1:parCount[[i]]]
    
    ind1<-which(is.finite(lowerBounds[[i]][parindex[[i]]+1:parCount[[i]]]) & is.infinite(upperBounds[[i]][parindex[[i]]+1:parCount[[i]]]))
    ind2<-which(is.finite(lowerBounds[[i]][parindex[[i]]+1:parCount[[i]]]) & is.finite(upperBounds[[i]][parindex[[i]]+1:parCount[[i]]]))
    ind3<-which(is.infinite(lowerBounds[[i]][parindex[[i]]+1:parCount[[i]]]) & is.finite(upperBounds[[i]][parindex[[i]]+1:parCount[[i]]]))
    
    tmpwpar[ind1] <- exp(tmpwpar[ind1])+lowerBounds[[i]][parindex[[i]]+ind1]
    tmpwpar[ind2] <- (upperBounds[[i]][parindex[[i]]+ind2]-lowerBounds[[i]][parindex[[i]]+ind2]) * boot::inv.logit(tmpwpar[ind2])+lowerBounds[[i]][parindex[[i]]+ind2]
    tmpwpar[ind3] <- -(exp(-tmpwpar[ind3]) - upperBounds[[i]][parindex[[i]]+ind3])
    if(any(tmpwpar[parindex[[i]]+1:parCount[[i]]]<lowerBounds[[i]] | tmpwpar[parindex[[i]]+1:parCount[[i]]]>upperBounds[[i]])) stop("starting values for ",i," data stream do not satisfy workBounds")
  }
  
  if(!is.null(beta0)){
    if(is.null(workBounds$beta)){
      lowerBounds$beta <- rep(-Inf,length(beta0))
      upperBounds$beta <- rep(Inf,length(beta0))
    } else {
      if(!is.matrix(workBounds$beta)) stop("workBounds$beta must be a matrix")
      if(!all(dim(workBounds$beta)==c(length(beta0),2))) stop("workBounds$",i," must be a 2 x ",length(beta0)," matrix")
      lowerBounds$beta <- workBounds$beta[,1]
      upperBounds$beta <- workBounds$beta[,2]
    }
    ind1<-which(is.finite(lowerBounds$beta) & is.infinite(upperBounds$beta))
    ind2<-which(is.finite(lowerBounds$beta) & is.finite(upperBounds$beta))
    ind3<-which(is.infinite(lowerBounds$beta) & is.finite(upperBounds$beta))
    
    beta0[ind1] <- exp(beta0[ind1])+lowerBounds$beta[ind1]
    beta0[ind2] <- (upperBounds$beta[ind2]-lowerBounds$beta[ind2]) * boot::inv.logit(beta0[ind2])+lowerBounds$beta[ind2]
    beta0[ind3] <- -(exp(-beta0[ind3]) - upperBounds$beta[ind3])
    if(any(c(beta0)<lowerBounds$beta | c(beta0)>upperBounds$beta)) stop("beta0 do not satisfy workBounds")
  }

  if(!is.null(delta0)){
    if(is.null(workBounds$delta)){
      lowerBounds$delta <- rep(-Inf,length(delta0))
      upperBounds$delta <- rep(Inf,length(delta0))
    } else {
      if(!is.matrix(workBounds$delta)) stop("workBounds$delta must be a matrix")
      if(!all(dim(workBounds$delta)==c(length(delta0),2))) stop("workBounds$",i," must be a 2 x ",length(delta0)," matrix")
      lowerBounds$delta <- workBounds$delta[,1]
      upperBounds$delta <- workBounds$delta[,2]
    }
    ind1<-which(is.finite(lowerBounds$delta) & is.infinite(upperBounds$delta))
    ind2<-which(is.finite(lowerBounds$delta) & is.finite(upperBounds$delta))
    ind3<-which(is.infinite(lowerBounds$delta) & is.finite(upperBounds$delta))
    
    delta0[ind1] <- exp(delta0[ind1])+lowerBounds$delta[ind1]
    delta0[ind2] <- (upperBounds$delta[ind2]-lowerBounds$delta[ind2]) * boot::inv.logit(delta0[ind2])+lowerBounds$delta[ind2]
    delta0[ind3] <- -(exp(-delta0[ind3]) - upperBounds$delta[ind3])
    if(any(c(delta0)<lowerBounds$delta | c(delta0)>upperBounds$delta)) stop("delta0 do not satisfy workBounds")
  }
  
  lowerBounds <- unlist(lowerBounds[c(distnames,"beta","delta")])
  upperBounds <- unlist(upperBounds[c(distnames,"beta","delta")])

  list(lower=lowerBounds,upper=upperBounds)
}