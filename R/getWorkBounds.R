getWorkBounds <- function(workBounds,distnames,wpar,parindex,parCount,DM,beta=NULL,delta=NULL){
  
  for(i in distnames){
    if(is.null(workBounds[[i]]) | is.null(DM[[i]])){
      workBounds[[i]] <- matrix(c(-Inf,Inf),parCount[[i]],2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds[[i]])) stop("workBounds$",i," must be a matrix")
      if(!all(dim(workBounds[[i]])==c(parCount[[i]],2))) stop("workBounds$",i," must be a 2 x ",parCount[[i]]," matrix")
    }
    #tmpwpar <- w2wn(wpar[parindex[[i]]+1:parCount[[i]]],workBounds[[i]])
    #if(any(tmpwpar<workBounds[[i]][,1] | tmpwpar>workBounds[[i]][,2])) stop("starting values for ",i," data stream do not satisfy workBounds")
  }
  
  if(length(beta)){
    if(is.null(workBounds$beta)){
      workBounds$beta <- matrix(c(-Inf,Inf),length(beta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$beta)) stop("workBounds$beta must be a matrix")
      if(!all(dim(workBounds$beta)==c(length(beta),2))) stop("workBounds$",i," must be a 2 x ",length(beta)," matrix")
    }
    #beta <- w2wn(c(beta),workBounds$beta)
    #if(any(beta<workBounds$beta[,1] | beta>workBounds$beta[,2])) stop("beta0 do not satisfy workBounds")
  }

  if(length(delta)){
    if(is.null(workBounds$delta)){
      workBounds$delta <- matrix(c(-Inf,Inf),length(delta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$delta)) stop("workBounds$delta must be a matrix")
      if(!all(dim(workBounds$delta)==c(length(delta),2))) stop("workBounds$",i," must be a 2 x ",length(delta)," matrix")
    }
    #delta <- w2wn(c(delta),workBounds$delta)
    #if(any(delta<workBounds$delta[,1] | delta>workBounds$delta[,2])) stop("delta0 do not satisfy workBounds")
  }
  workBounds
}