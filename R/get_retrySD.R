get_retrySD <- function(retryFits,retrySD,wpar,parmInd,distnames,parCount,nbStates,beta0,mixtures,recharge,delta0){
  if(retryFits<0) stop("retryFits must be non-negative")
  if(retryFits>=1){
    if(is.null(retrySD)){
      retrySD <- rep(10,length(wpar))
      retrySD[1:parmInd] <- 1
    } else if(!is.list(retrySD)){
      if(length(retrySD)>1) stop('retrySD must be a scalar or list')
      retrySD <- rep(retrySD,length(wpar))
    } else {
      for(i in distnames){
        if(is.null(retrySD[[i]])){
          retrySD[[i]] <- rep(1,parCount[[i]])
        } else {
          if(any(!is.numeric(retrySD[[i]]))) stop("retrySD$",i," must be numeric")
          if(any(retrySD[[i]]<0)) stop("retrySD$",i," must be non-negative")
          if(length(retrySD[[i]])==1){
            retrySD[[i]] <- rep(retrySD[[i]],parCount[[i]])
          } else if(length(retrySD[[i]])!=parCount[[i]]){
            stop("retrySD$",i," must be of length 1 or ",parCount[[i]])
          }
        }
      }
      if(nbStates>1){
        if(is.null(retrySD$beta)){
          retrySD$beta <- rep(10,length(beta0$beta))
        } else {
          if(any(!is.numeric(retrySD$beta))) stop("retrySD$beta must be numeric")
          if(any(retrySD$beta<0)) stop("retrySD$beta must be non-negative")
          if(length(retrySD$beta)==1){
            retrySD$beta <- rep(retrySD$beta,length(beta0$beta))
          } else if(length(retrySD$beta)!=length(beta0$beta)){
            stop("retrySD$beta must be of length 1 or ",length(beta0$beta))
          }
        }
        if(mixtures>1){
          if(is.null(retrySD$pi)) {
            retrySD$pi <- rep(10,length(beta0$pi))
          } else {
            if(any(!is.numeric(retrySD$pi))) stop("retrySD$pi must be numeric")
            if(any(retrySD$pi<0)) stop("retrySD$pi must be non-negative")
            if(length(retrySD$pi)==1){
              retrySD$pi <- rep(retrySD$pi,length(beta0$pi))
            } else if(length(retrySD$pi)!=length(beta0$pi)){
              stop("retrySD$pi must be of length 1 or ",length(beta0$pi))
            }              
          }
        }
        if(!is.null(recharge)){
          for(j in c("g0","theta")){
            if(is.null(retrySD[[j]])) {
              retrySD[[j]] <- rep(10,length(beta0[[j]]))
            } else {
              if(any(!is.numeric(retrySD[[j]]))) stop("retrySD$",j," must be numeric")
              if(any(retrySD[[j]]<0)) stop("retrySD$",j," must be non-negative")
              if(length(retrySD[[j]])==1){
                retrySD[[j]] <- rep(retrySD[[j]],length(beta0[[j]]))
              } else if(length(retrySD[[j]])!=length(beta0[[j]])){
                stop("retrySD$",j," must be of length 1 or ",length(beta0[[j]]))
              }              
            }
          }
        }
        if(is.null(retrySD$delta)){
          retrySD$delta <- rep(10,length(delta0))
        } else {
          if(any(!is.numeric(retrySD$delta))) stop("retrySD$delta must be numeric")
          if(any(retrySD$delta<0)) stop("retrySD$delta must be non-negative")
          if(length(retrySD$delta)==1){
            retrySD$delta <- rep(retrySD$delta,length(delta0))
          } else if(length(retrySD$delta)!=length(delta0)){
            stop("retrySD$delta must be of length 1 or ",length(delta0))
          }
        }
      }
      retrySD <- unlist(retrySD[c(distnames,"beta","pi","delta","g0","theta")])
    }
  }
  retrySD
}