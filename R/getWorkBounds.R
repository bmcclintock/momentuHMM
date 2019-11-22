getWorkBounds <- function(workBounds,distnames,wpar,parindex,parCount,DM,beta=NULL,delta=NULL){
  
  for(i in distnames){
    if(is.null(workBounds[[i]]) | is.null(DM[[i]])){
      workBounds[[i]] <- matrix(c(-Inf,Inf),parCount[[i]],2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds[[i]])) stop("workBounds$",i," must be a matrix")
      if(!all(dim(workBounds[[i]])==c(parCount[[i]],2))) stop("workBounds$",i," must be a ",parCount[[i]]," x 2 matrix")
    }
    #tmpwpar <- w2wn(wpar[parindex[[i]]+1:parCount[[i]]],workBounds[[i]])
    #if(any(tmpwpar<workBounds[[i]][,1] | tmpwpar>workBounds[[i]][,2])) stop("starting values for ",i," data stream do not satisfy workBounds")
  }
  outBounds <- workBounds[distnames]
  
  if(length(beta$beta)){
    if(is.null(workBounds$beta)){
      workBounds$beta <- matrix(c(-Inf,Inf),length(beta$beta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$beta)) stop("workBounds$beta must be a matrix")
      if(!all(dim(workBounds$beta)==c(length(beta$beta),2))) stop("workBounds$beta must be a ",length(beta$beta)," x 2 matrix")
    }
    #beta <- w2wn(c(beta),workBounds$beta)
    #if(any(beta<workBounds$beta[,1] | beta>workBounds$beta[,2])) stop("beta0 do not satisfy workBounds")
    outBounds <- workBounds[c(distnames,"beta")]
  }
  if(length(beta$pi)){
    if(is.null(workBounds$pi)){
      workBounds$pi <- matrix(c(-Inf,Inf),length(beta$pi),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$pi)) stop("workBounds$pi must be a matrix")
      if(!all(dim(workBounds$pi)==c(length(beta$pi),2))) stop("workBounds$pi must be a ",length(beta$pi)," x 2 matrix")
    }
    outBounds <- workBounds[c(distnames,"beta","pi")]
  }
  if(length(delta)){
    if(is.null(workBounds$delta)){
      workBounds$delta <- matrix(c(-Inf,Inf),length(delta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$delta)) stop("workBounds$delta must be a matrix")
      if(!all(dim(workBounds$delta)==c(length(delta),2))) stop("workBounds$delta must be a ",length(delta)," x 2 matrix")
    }
    #delta <- w2wn(c(delta),workBounds$delta)
    #if(any(delta<workBounds$delta[,1] | delta>workBounds$delta[,2])) stop("delta0 do not satisfy workBounds")
    outBounds <- workBounds[c(distnames,"beta","pi","delta")]
  }
  if(length(beta$theta)){
    if(is.null(workBounds$g0)){
      workBounds$g0 <- matrix(c(-Inf,Inf),length(beta$g0),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$g0)) stop("workBounds$g0 must be a matrix")
      if(!all(dim(workBounds$g0)==c(length(beta$g0),2))) stop("workBounds$g0 must be a ",length(beta$g0)," x 2 matrix")
    }
    if(is.null(workBounds$theta)){
      workBounds$theta <- matrix(c(-Inf,Inf),length(beta$theta),2,byrow=TRUE)
    } else {
      if(!is.matrix(workBounds$theta)) stop("workBounds$theta must be a matrix")
      if(!all(dim(workBounds$theta)==c(length(beta$theta),2))) stop("workBounds$theta must be a ",length(beta$theta)," x 2 matrix")
    }
    outBounds <- workBounds[c(distnames,"beta","pi","delta","g0","theta")]
  }
  outBounds
}