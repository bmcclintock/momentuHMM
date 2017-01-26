
#' Confidence intervals
#'
#' Computes the standard errors and confidence intervals on the real (i.e., natural) scale of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations in the computation of the CIs for the angle parameters.
#' Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{stepPar}{Standard errors and confidence intervals for the parameters of the step lengths distribution}
#' \item{anglePar}{Standard errors and confidence intervals for the parameters of the turning angles distribution}
#' \item{omegaPar}{Standard errors and confidence intervals for the parameters of the omega distribution}
#' \item{dryPar}{Standard errors and confidence intervals for the parameters of the dry distribution}
#' \item{divePar}{Standard errors and confidence intervals for the parameters of the dive distribution}
#' \item{icePar}{Standard errors and confidence intervals for the parameters of the ice distribution}
#' \item{landPar}{Standard errors and confidence intervals for the parameters of the land distribution}
#' \item{gamma}{Standard errors and confidence intervals for the transition probabilities. If covariates are included in TPM formula, then the covariate values for the first observation in the data are used.}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CI_real(m)
#'
#' @export
#' @importFrom MASS ginv
#' @importFrom numDeriv grad
#' @importFrom utils tail

CI_real <- function(m,alpha=0.95,nbSims=10^6)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- length(m$stateNames)

  dist <- m$conditions$dist
  distnames <- names(dist)
  fullDM <- m$conditions$fullDM
  DMind <- m$conditions$DMind

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)

  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds)
  bounds <- p$bounds
  #if(!all(unlist(lapply(p$bounds,is.numeric)))){
  #  for(i in distnames){
  #    if(!is.numeric(bounds[[i]])){
  #      bounds[[i]] <- gsub(i,"",bounds[[i]],fixed=TRUE)
  #    }
  #  }
  #}
  
  parindex <- c(0,cumsum(unlist(lapply(fullDM,ncol)))[-length(fullDM)])
  names(parindex) <- distnames
  
  Par <- list()
  lower<-list()
  upper<-list()
  se<-list()
  
  for(i in distnames){
    if(!m$conditions$DMind[[i]]){
      tmpDM<-fullDM[[i]]
      k <- which(matrix(mapply(length,fullDM[[i]])>1,nrow(fullDM[[i]]),ncol(fullDM[[i]])),arr.ind=TRUE)
      if(length(k)){
        for(j in 1:nrow(k)){
          tmpDM[[k[j,1],k[j,2]]]<-mean(fullDM[[i]][[k[j,1],k[j,2]]],na.rm=TRUE)
        }
      }
      fullDM[[i]]<-matrix(as.numeric(tmpDM),nrow(tmpDM),ncol(tmpDM))
      DMind[[i]]<-TRUE
      par <- c(w2n(m$mod$estimate,bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$stationary,m$conditions$cons,fullDM,DMind,m$conditions$workcons,1,dist[i],m$conditions$Bndind)[[i]])
    } else {
      par <- as.vector(t(m$mle[[i]]))
    }
    if(!(dist[[i]] %in% angledists) | (dist[[i]] %in% angledists & m$conditions$estAngleMean[[i]] & !m$conditions$Bndind[[i]])) {
      Par[[i]] <- get_CI(m$mod$estimate,par,m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],p$parSize[[i]],Sigma,nbStates,alpha,p$parNames[[i]],m$stateNames)
    } else {
      if(!m$conditions$estAngleMean[[i]])
        Par[[i]] <- get_CI(m$mod$estimate,par[-c(1:nbStates)],m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],p$parSize[[i]],Sigma,nbStates,alpha,p$parNames[[i]],m$stateNames)
      else {
        if(m$conditions$Bndind[[i]]){
          Par[[i]] <- CI_angle(m$mod$estimate,par,m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],p$parSize[[i]],Sigma,nbStates,alpha,p$parNames[[i]],m$stateNames)
        }
      }
    }
  }

  if(nbStates>1 & !is.null(m$mle$gamma)) {
    # identify parameters of interest
    i2 <- tail(cumsum(unlist(lapply(fullDM,ncol))),1)+1
    i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)-1
    wpar <- m$mle$beta
    quantSup <- qnorm(1-(1-alpha)/2)
    lower<-upper<-se<-matrix(NA,nbStates,nbStates)
    for(i in 1:nbStates){
      for(j in 1:nbStates){
        dN<-numDeriv::grad(get_gamma,wpar,covs=covs,nbStates=nbStates,i=i,j=j)
        se[i,j]<-suppressWarnings(sqrt(dN%*%Sigma[i2:i3,i2:i3]%*%dN))
        lower[i,j]<-1/(1+exp(-(log(m$mle$gamma[i,j]/(1-m$mle$gamma[i,j]))-quantSup*(1/(m$mle$gamma[i,j]-m$mle$gamma[i,j]^2))*se[i,j])))#m$mle$gamma[i,j]-quantSup*se[i,j]
        upper[i,j]<-1/(1+exp(-(log(m$mle$gamma[i,j]/(1-m$mle$gamma[i,j]))+quantSup*(1/(m$mle$gamma[i,j]-m$mle$gamma[i,j]^2))*se[i,j])))#m$mle$gamma[i,j]+quantSup*se[i,j]
      }
    }
    Par$gamma <- list(est=m$mle$gamma,se=se,lower=lower,upper=upper)
    rownames(Par$gamma$est) <- m$stateNames
    rownames(Par$gamma$se) <- m$stateNames
    rownames(Par$gamma$lower) <- m$stateNames
    rownames(Par$gamma$upper) <- m$stateNames
    colnames(Par$gamma$est) <- m$stateNames
    colnames(Par$gamma$se) <- m$stateNames
    colnames(Par$gamma$lower) <- m$stateNames
    colnames(Par$gamma$upper) <- m$stateNames
  }
  
  wpar<-m$mod$estimate
  foo <- length(wpar)-nbStates+2
  quantSup <- qnorm(1-(1-alpha)/2)
  lower<-upper<-se<-rep(NA,nbStates)
  for(i in 1:nbStates){
    dN<-numDeriv::grad(get_delta,wpar[foo:length(wpar)],i=i)
    se[i]<-suppressWarnings(sqrt(dN%*%Sigma[foo:length(wpar),foo:length(wpar)]%*%dN))
    lower[i]<-1/(1+exp(-(log(m$mle$delta[i]/(1-m$mle$delta[i]))-quantSup*(1/(m$mle$delta[i]-m$mle$delta[i]^2))*se[i])))#m$mle$delta[i]-quantSup*se[i]
    upper[i]<-1/(1+exp(-(log(m$mle$delta[i]/(1-m$mle$delta[i]))+quantSup*(1/(m$mle$delta[i]-m$mle$delta[i]^2))*se[i])))#m$mle$delta[i]+quantSup*se[i]
  }
  est<-matrix(m$mle$delta,nrow=1,ncol=nbStates,byrow=TRUE)
  lower<-matrix(lower,nrow=1,ncol=nbStates,byrow=TRUE)
  upper<-matrix(upper,nrow=1,ncol=nbStates,byrow=TRUE)  
  se<-matrix(se,nrow=1,ncol=nbStates,byrow=TRUE)
  Par$delta <- list(est=est,se=se,lower=lower,upper=upper)  
  colnames(Par$delta$est) <- m$stateNames
  colnames(Par$delta$se) <- m$stateNames
  colnames(Par$delta$lower) <- m$stateNames
  colnames(Par$delta$upper) <- m$stateNames

  return(Par)
}

get_gamma <- function(beta,covs,nbStates,i,j){
  gamma <- trMatrix_rcpp(nbStates,beta,covs)[,,1]
  gamma[i,j]
}

get_delta <- function(delta,i){
  delta <- exp(c(0,delta))
  delta <- delta/sum(delta)
  delta[i]
}

parm_list<-function(est,se,lower,upper,rnames,cnames){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  rownames(Par$est) <- rnames
  rownames(Par$se) <- rnames
  rownames(Par$lower) <- rnames
  rownames(Par$upper) <- rnames
  colnames(Par$est) <- cnames
  colnames(Par$se) <- cnames
  colnames(Par$lower) <- cnames
  colnames(Par$upper) <- cnames
  Par
}

get_CI<-function(wpar,Par,m,ind,DM,DMind,Bounds,cons,workcons,parSize,Sigma,nbStates,alpha,rnames,cnames){

  w<-wpar[ind]
  lower<-upper<-se<-numeric(nrow(DM))
  for(k in 1:nrow(DM)){
    dN<-numDeriv::grad(w2nDM,w,bounds=Bounds,DM=DM,DMind=DMind,cons=cons,workcons=workcons,nbObs=1,parSize=parSize,k=k)
    se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
    lower[k] <- Par[k] - qnorm(1-(1-alpha)/2) * se[k]
    upper[k] <- Par[k] + qnorm(1-(1-alpha)/2) * se[k]
    #cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
    #lower[k]<-Par[k]/cn
    #upper[k]<-Par[k]*cn
  }
  est<-matrix(Par,ncol=nbStates,byrow=TRUE)
  l<-matrix(lower,ncol=nbStates,byrow=TRUE)
  u<-matrix(upper,ncol=nbStates,byrow=TRUE)  
  s<-matrix(se,ncol=nbStates,byrow=TRUE)
  out <- parm_list(est,s,l,u,rnames,cnames)
  out
}

CI_angle<-function(wpar,Par,m,ind,DM,DMind,Bounds,cons,workcons,parSize,Sigma,nbStates,alpha,rnames,cnames){
  
  w<-wpar[ind]
  lower<-upper<-se<-numeric(nrow(DM))
  for(k in 1:nrow(DM)){
    dN<-numDeriv::grad(w2nDMangle,w,bounds=Bounds,DM=DM,DMind=DMind,cons=cons,workcons=workcons,nbObs=1,parSize=parSize,nbStates=nbStates,k=k)
    se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
    lower[k] <- Par[k] - qnorm(1-(1-alpha)/2) * se[k]
    upper[k] <- Par[k] + qnorm(1-(1-alpha)/2) * se[k]
    #cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
    #lower[k]<-Par[k]/cn
    #upper[k]<-Par[k]*cn
  }
  est<-matrix(Par,ncol=nbStates,byrow=TRUE)
  l<-matrix(lower,ncol=nbStates,byrow=TRUE)
  u<-matrix(upper,ncol=nbStates,byrow=TRUE)  
  s<-matrix(se,ncol=nbStates,byrow=TRUE)
  out <- parm_list(est,s,l,u,rnames,cnames)
  out
}
  
w2nDMangle<-function(w,bounds,DM,DMind,cons,workcons,nbObs,parSize,nbStates,k){
  
  bounds[,1] <- -Inf
  bounds[which(bounds[,2]!=1),2] <- Inf
    
  foo <- length(w) - nbStates + 1
  x <- w[(foo - nbStates):(foo - 1)]
  y <- w[foo:length(w)]
  angleMean <- Arg(x + (0+1i) * y)
  kappa <- sqrt(x^2 + y^2)
  w[(foo - nbStates):(foo - 1)] <- angleMean
  w[foo:length(w)] <- kappa
  
  w2nDM(w,bounds,DM,DMind,cons,workcons,nbObs,parSize,k)
}