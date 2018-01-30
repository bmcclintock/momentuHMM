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
    if(any(wpar[parindex[[i]]+1:parCount[[i]]]<lowerBounds[[i]] | wpar[parindex[[i]]+1:parCount[[i]]]>upperBounds[[i]])) stop("starting values for ",i," data stream are not within workBounds")
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
    if(any(c(beta0)<lowerBounds$beta | c(beta0)>upperBounds$beta)) stop("beta0 are not within workBounds")
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
    if(any(c(delta0)<lowerBounds$delta | c(delta0)>upperBounds$delta)) stop("delta0 are not within workBounds")
  }
  
  lowerBounds <- unlist(lowerBounds[c(distnames,"beta","delta")])
  upperBounds <- unlist(upperBounds[c(distnames,"beta","delta")])

  list(lower=lowerBounds,upper=upperBounds)
}