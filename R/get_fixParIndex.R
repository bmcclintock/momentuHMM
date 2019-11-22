get_fixParIndex <- function(Par0,beta0,delta0,fixPar,distnames,inputs,p,nbStates,DMinputs,recharge,nbG0covs,nbRecovs,workBounds,mixtures,newformula,formulaStates,nbCovs,betaCons,betaRef,deltaCons,stationary,nbCovsDelta,formulaDelta,formulaPi,nbCovsPi,data,newdata){
  
  if(nbStates>1){
    if(!is.null(betaCons)){
      if(!is.matrix(betaCons)) stop("betaCons must be a matrix")
      if(nrow(betaCons)!=(nbCovs+1)*mixtures | ncol(betaCons)!=(nbStates*(nbStates-1))) stop("betaCons must be a ",(nbCovs+1)*mixtures,"x",nbStates*(nbStates-1)," matrix")
      if(any(abs(as.integer(betaCons)-betaCons)!=0)) stop("betaCons must be a matrix composed of integers")
      if(min(betaCons)<1 | max(betaCons)>(nbCovs+1)*mixtures*(nbStates*(nbStates-1))) stop("betaCons must be composed of integers between 1 and ",(nbCovs+1)*mixtures*(nbStates*(nbStates-1)))
    }
    if(!is.null(deltaCons)){
      if(is.null(formulaDelta)) stop("deltaCons cannot be specified unless formulaDelta includes a formula")
      if(!is.matrix(deltaCons)) stop("deltaCons must be a matrix")
      if(nrow(deltaCons)!=(nbCovsDelta+1)*mixtures | ncol(deltaCons)!=(nbStates-1)) stop("deltaCons must be a ",(nbCovsDelta+1)*mixtures,"x",(nbStates-1)," matrix")
      if(any(abs(as.integer(deltaCons)-deltaCons)!=0)) stop("deltaCons must be a matrix composed of integers")
      if(min(deltaCons)<1 | max(deltaCons)>(nbCovsDelta+1)*mixtures*(nbStates-1)) stop("deltaCons must be composed of integers between 1 and ",(nbCovsDelta+1)*mixtures**(nbStates-1))
    }
  }
  
  if(!is.null(betaRef)){
    if(length(betaRef)!=nbStates) stop("betaRef must be a vector of length ",nbStates)
    if(!is.numeric(betaRef)) stop("betaRef must be a numeric vector")
    if(min(betaRef)<1 | max(betaRef)>nbStates) stop("betaRef elements must be between 1 and ",nbStates)
  } else {
    betaRef <- 1:nbStates
  }
  betaRef <- as.integer(betaRef)
  
  if(!is.null(beta0$beta)) {
    if(is.null(dim(beta0$beta)))
      stop(paste(paste0("beta",ifelse(is.null(recharge),"0","0$beta"))," has wrong dimensions: it should have",(nbCovs+1)*mixtures,"rows and",
                 nbStates*(nbStates-1),"columns."))
    if(ncol(beta0$beta)!=nbStates*(nbStates-1) | nrow(beta0$beta)!=(nbCovs+1)*mixtures)
      stop(paste(paste0("beta",ifelse(is.null(recharge),"0","0$beta"))," has wrong dimensions: it should have",(nbCovs+1)*mixtures,"rows and",
                 nbStates*(nbStates-1),"columns."))
    if(!is.null(betaCons)) {
      if(!all(beta0$beta == beta0$beta[c(betaCons)])) warning(paste0("beta",ifelse(is.null(recharge),"0","0$beta"))," is not consistent with betaCons; values for ",paste0("beta",ifelse(is.null(recharge),"0","0$beta"))," will be assigned based on the first duplicated element(s) in betaCons")
      beta0$beta <- matrix(beta0$beta[c(betaCons)],nrow(beta0$beta),ncol(beta0$beta))
    }
  } else {
    if(nbStates>1){
      beta0$beta <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                           nrow=(nbCovs+1)*mixtures,ncol=nbStates*(nbStates-1),byrow=TRUE)
      for(i in 1:nbStates){
        if(betaRef[i]!=i){
          beta0$beta[1,(i-1)*(nbStates-1)+i-(betaRef[i]<i)] <- 1.5
        }
      }
      if(!is.null(betaCons)) beta0$beta <- matrix(beta0$beta[c(betaCons)],nrow(beta0$beta),ncol(beta0$beta))
    }
  }
  if(!is.null(recharge)){
    if(is.null(beta0$g0)) beta0$g0 <- nw2w(rep(0,nbG0covs+1),workBounds$g0)
    else if(length(beta0$g0)!=(nbG0covs+1) | any(!is.numeric(beta0$g0))) stop("beta0$g0 must be a numeric vector of length ",nbG0covs+1)
    if(is.null(beta0$theta)) beta0$theta <- nw2w(c(-1,rep(1,nbRecovs)),workBounds$theta)
    else if(length(beta0$theta)!=(nbRecovs+1) | any(!is.numeric(beta0$theta))) stop("beta0$theta must be a numeric vector of length ",nbRecovs+1)
  }
  
  if(mixtures<=1){
    beta0$pi <- NULL
  } else {
    if(!is.null(beta0$pi)){
      if(!nbCovsPi & is.null(formulaPi)){
        if(length(beta0$pi) != (nbCovsPi+1)*mixtures)
          stop(paste("beta0$pi has the wrong length: it should have",mixtures,"elements."))
      } else {
        if(is.null(dim(beta0$pi)) || (ncol(beta0$pi)!=mixtures-1 | nrow(beta0$pi)!=(nbCovsPi+1)))
          stop(paste("beta0$pi has wrong dimensions: it should have",(nbCovsPi+1),"rows and",
                     mixtures-1,"columns."))
      }
    } 
  }
  
  if(!is.null(delta0)){
    if(!nbCovsDelta & is.null(formulaDelta)){
      if(mixtures==1){
        if(length(delta0) != (nbCovsDelta+1)*nbStates)
          stop(paste("delta0 has the wrong length: it should have",nbStates*mixtures,"elements."))
      } else {
        if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates | nrow(delta0)!=mixtures))
          stop(paste("delta0 has wrong dimensions: it should have",mixtures,"rows and",
                     nbStates,"columns."))
      }
    } else {
      if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates-1 | nrow(delta0)!=(nbCovsDelta+1)*mixtures))
        stop(paste("delta0 has wrong dimensions: it should have",(nbCovsDelta+1)*mixtures,"rows and",
                   nbStates-1,"columns."))
      if(!is.null(deltaCons)) {
        if(!all(delta0 == delta0[c(deltaCons)])) warning("delta0 is not consistent with deltaCons; values for delta0 will be assigned based on the first duplicated element(s) in deltaCons")
        delta0 <- matrix(delta0[c(deltaCons)],nrow(delta0),ncol(delta0))
      }
    }
  }
  
  if(stationary) 
    delta0 <- NULL
  else if(!nbCovsDelta & is.null(formulaDelta)){
    if(is.null(delta0)){
      delta0 <- matrix(1/nbStates,(nbCovsDelta+1)*mixtures,nbStates)
    } else {
      delta0 <- matrix(delta0,(nbCovsDelta+1)*mixtures,nbStates)
    }
    delta0 <- t(apply(delta0,1,function(x) log(x[-1]/x[1])))
  } else if(is.null(delta0)){
    delta0 <- matrix(0,nrow=(nbCovsDelta+1)*mixtures,ncol=nbStates-1)
  }
  
  if(!nbCovsPi & is.null(formulaPi)){
    if(is.null(beta0$pi)){
      beta0$pi <- matrix(1/mixtures,(nbCovsPi+1),mixtures)
    } else {
      beta0$pi <- matrix(beta0$pi,(nbCovsPi+1),mixtures)
    }
    beta0$pi <- log(beta0$pi[-1]/beta0$pi[1])
  } else if(is.null(beta0$pi)){
    beta0$pi <- matrix(0,nrow=(nbCovsPi+1),ncol=mixtures-1)
  }
  
  # build the vector of initial working parameters
  wparIndex <- numeric()
  if(is.null(fixPar)){
    fixPar <- vector('list',length(distnames))
    names(fixPar) <- distnames
    for(i in distnames){
      fixPar[[i]] <- rep(NA,length(Par0[[i]]))
    }
  }
  
  if(is.null(fixPar$beta)) {
    fixPar$beta <- rep(NA,length(beta0$beta))
  }
  if(!is.null(recharge)){
    if(is.null(fixPar$g0)) {
      fixPar$g0 <- rep(NA,nbG0covs+1)
    } else if(length(fixPar$g0)!=(nbG0covs+1)) stop("fixPar$g0 must be a vector of length ",nbG0covs+1)
    if(is.null(fixPar$theta)){
      fixPar$theta <- rep(NA,nbRecovs+1)
    } else if(length(fixPar$theta)!=(nbRecovs+1)) stop("fixPar$theta must be a vector of length ",nbRecovs+1)
  } else {
    fixPar$g0 <- fixPar$theta <- NULL
  }
  
  if(nbStates>1){
    if(!is.null(newdata)) data <- cbind(data,newdata)
    for(state in 1:(nbStates*(nbStates-1))){
      noBeta<-which(match(colnames(model.matrix(newformula,data)),colnames(model.matrix(formulaStates[[state]],data)),nomatch=0)==0)
      if(length(noBeta)) {
        for(mix in 1:mixtures){
          beta0$beta[noBeta+(nbCovs+1)*(mix-1),state] <- NA
          fixPar$beta[noBeta+(state-1)*(nbCovs+1)+(mix-1)*(nbCovs+1)*(nbStates*(nbStates-1))] <- 0
        }
      }
    }
  }
  
  parindex <- c(0,cumsum(unlist(lapply(Par0,length))))
  names(parindex) <- c(distnames,"beta")
  ofixPar <- fixPar
  for(i in distnames){
    if(!is.null(fixPar[[i]])){
      if(length(fixPar[[i]])!=length(Par0[[i]])) stop("fixPar$",i," must be of length ",length(Par0[[i]]))
      tmp <- which(!is.na(fixPar[[i]]))
      Par0[[i]][tmp]<-fixPar[[i]][tmp]
      
      # if DM is not specified convert fixPar from real to working scale
      if(is.null(inputs$DM[[i]])){
        fixPar[[i]][tmp]<-n2w(Par0[i],p$bounds,NULL,NULL,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind,inputs$dist)[tmp]
      } 
      wparIndex <- c(wparIndex,parindex[[i]]+tmp)
    } else {
      fixPar[[i]] <- ofixPar[[i]] <- rep(NA,length(Par0[[i]]))
      #wparIndex <- c(wparIndex,parindex[[i]]+1:length(Par0[[i]]))
    }
  }
  if(nbStates>1){
    beta0$beta[which(is.na(beta0$beta))] <- 0
    if(length(fixPar$beta)!=length(beta0$beta)) stop("fixPar$beta must be of length ",length(beta0$beta))
    tmp <- which(!is.na(fixPar$beta))
    if(!all(fixPar$beta == fixPar$beta[betaCons],na.rm=TRUE)) stop("fixPar$beta not consistent with betaCons")
    beta0$beta[tmp]<-fixPar$beta[tmp]
    wparIndex <- c(wparIndex,parindex[["beta"]]+tmp)
    
    if(mixtures>1){
      if(any(!is.na(fixPar$pi))){
        tmp <- which(!is.na(fixPar$pi))
        beta0$pi <- fixPar$pi
        if(length(tmp)!=length(beta0$pi) | sum(beta0$pi)!=1) stop("fixPar$pi must sum to 1")
        beta0$pi <- matrix(beta0$pi,1,mixtures)
        beta0$pi <- log(beta0$pi[-1]/beta0$pi[1])
        fixPar$pi <- as.vector(beta0$pi)
        wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+1:length(beta0$pi))
      }
    }
  }
  if(!is.null(fixPar$pi)){
    if(!nbCovsPi & is.null(formulaPi)){
      if(length(fixPar$pi)!=mixtures) stop("fixPar$pi must be of length ",mixtures)
    } else if(length(fixPar$pi)!=length(beta0$pi)) stop("fixPar$pi must be of length ",length(beta0$pi))
    if(any(!is.na(fixPar$pi))){
      tmp <- which(!is.na(fixPar$pi))
      if(!nbCovsPi & is.null(formulaPi)){
        beta0$pi <- matrix(fixPar$pi,1,mixtures)
        if(length(tmp)!=length(beta0$pi) | !all(rowSums(beta0$pi)==1)) stop("fixPar$pi must sum to 1")
        beta0$pi <- matrix(log(beta0$pi[-1]/beta0$pi[1]),(nbCovsPi+1),mixtures-1)
        fixPar$pi <- as.vector(beta0$pi)
        wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+1:(length(beta0$pi)))
      } else {
        beta0$pi[tmp] <- fixPar$pi[tmp]
        wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+tmp)
      }
    } else fixPar$pi <- rep(NA,length(beta0$pi))
  } else {
    if(!nbCovsPi & is.null(formulaPi)) fixPar$pi <- ofixPar$pi <- rep(NA,ifelse(!length(beta0$pi),0,mixtures-1))
    else fixPar$pi <- ofixPar$pi <- rep(NA,length(beta0$pi))
  }
  if(!is.null(fixPar$delta)){
    if(stationary & any(!is.na(fixPar$delta))) stop("delta cannot be fixed when stationary=TRUE")
    else if(!stationary) {
      if(!nbCovsDelta & is.null(formulaDelta)){
        if(length(fixPar$delta)!=nbStates*mixtures) stop("fixPar$delta must be of length ",nbStates*mixtures)
      } else if(length(fixPar$delta)!=length(delta0)) stop("fixPar$delta must be of length ",length(delta0))
      if(any(!is.na(fixPar$delta))){
        tmp <- which(!is.na(fixPar$delta))
        if(!nbCovsDelta & is.null(formulaDelta)){
          delta0 <- matrix(fixPar$delta,mixtures,nbStates)
          if(length(tmp)!=length(delta0) | !all(mapply(function(x) isTRUE(all.equal(x,1)),rowSums(delta0)))) stop("fixPar$delta",ifelse(mixtures>1," rows "," "),"must sum to 1")
          delta0 <- matrix(apply(delta0,1,function(x) log(x[-1]/x[1])),(nbCovsDelta+1)*mixtures,nbStates-1)
          fixPar$delta <- as.vector(delta0)
          wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+length(beta0$pi)+1:(length(delta0)))
        } else {
          if(!all(fixPar$delta == fixPar$delta[deltaCons],na.rm=TRUE)) stop("fixPar$delta not consistent with deltaCons")
          delta0[tmp] <- fixPar$delta[tmp]
          wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+length(beta0$pi)+tmp)
        }
      } else fixPar$delta <- rep(NA,length(delta0))
    }
  } else {
    if(!nbCovsDelta & is.null(formulaDelta)) fixPar$delta <- ofixPar$delta <- rep(NA,ifelse(!length(delta0),0,(nbStates-1)*mixtures))
    else fixPar$delta <- ofixPar$delta <- rep(NA,length(delta0))
  }
  if(!is.null(recharge)){
    if(any(!is.na(fixPar$g0))){
      tmp <- which(!is.na(fixPar$g0))
      beta0$g0[tmp] <- fixPar$g0[tmp]
      wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+length(beta0$pi)+length(delta0)+tmp)
    }
    if(any(!is.na(fixPar$theta))){
      tmp <- which(!is.na(fixPar$theta))
      beta0$theta[tmp] <- fixPar$theta[tmp]
      wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0$beta)+length(beta0$pi)+length(delta0)+length(beta0$g0)+tmp)
    }
  }
  
  if(any(!(names(fixPar) %in% c(distnames,"beta","pi","delta","g0","theta")))) stop("fixPar names can only include ",paste0(c(distnames,"beta","pi","delta","g0"),sep=", "),"or theta")
  fixPar <- fixPar[c(distnames,"beta","pi","delta","g0","theta")]
  ofixPar <- ofixPar[c(distnames,"beta","pi","delta","g0","theta")]
  
  return(list(fixPar=fixPar,ofixPar=ofixPar,wparIndex=wparIndex,Par0=Par0,beta0=beta0,delta0=delta0,betaRef=betaRef))
}