
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
  DM <- m$conditions$DM

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)

  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$bounds,DM)
  bounds <- p$bounds
  if(!all(unlist(lapply(p$bounds,is.numeric)))){
    for(i in distnames){
      if(!is.numeric(bounds[[i]])){
        bounds[[i]] <- gsub(i,"",bounds[[i]],fixed=TRUE)
      }
    }
  }
  
  parindex <- c(0,cumsum(unlist(lapply(DM,ncol)))[-length(DM)])
  names(parindex) <- distnames
  
  Par <- list()
  lower<-list()
  upper<-list()
  se<-list()
  
  for(i in distnames){
    if(!(dist[[i]] %in% c("wrpcauchy","vm"))) {
      Par[[i]] <- get_CI(as.vector(t(m$mle[[i]])),m,parindex[[i]]+1:ncol(DM[[i]]),DM[[i]],bounds[[i]],m$conditions$cons[[i]],p$boundInd[[i]],logitcons=m$conditions$logitcons[[i]],Sigma,nbStates,alpha,m$mle[[i]])
    } else {
      wpar<-m$mod$estimate[parindex[[i]]+1:ncol(DM[[i]])]
      lower<-upper<-se<-numeric(nrow(DM[[i]]))
      if(m$conditions$estAngleMean[[i]]){
        for(k in 1:nrow(DM[[i]])){
          dN<-numDeriv::grad(w2nDM,wpar,bounds=bounds[[i]],DM=DM[[i]],cons=m$conditions$cons[[i]],boundInd=p$boundInd[[i]],logitcons=m$conditions$logitcons[[i]],k=k)
          se[k]<-suppressWarnings(sqrt(dN%*%Sigma[parindex[[i]]+1:ncol(DM[[i]]),parindex[[i]]+1:ncol(DM[[i]])]%*%dN))
          lower[k]<-t(m$mle[[i]])[k]-se[k]*qnorm(1-(1-alpha)/2)
          upper[k]<-t(m$mle[[i]])[k]+se[k]*qnorm(1-(1-alpha)/2)
          #if(k <= nbStates){
          #  l <- 2*atan(tan(lower[k]/2))
          #  u <- 2*atan(tan(upper[k]/2))
          #  if(l>u){
          #    lower[k]<-u
          #    upper[k]<-l
          #  } else {
          #    lower[k]<-l
          #    upper[k]<-u
          #  }
          #}
        }
        lower<-matrix(lower,ncol=nbStates,byrow=T)
        upper<-matrix(upper,ncol=nbStates,byrow=T)  
        se<-matrix(se,ncol=nbStates,byrow=T)
        Par[[i]] <- parm_list(se,lower,upper,m$mle[[i]])
      } else {
        for(k in 1:nrow(DM[[i]])){
          dN<-numDeriv::grad(w2nDM,wpar,bounds=bounds[[i]],DM=DM[[i]],cons=m$conditions$cons[[i]],boundInd=p$boundInd[[i]],logitcons=m$conditions$logitcons[[i]],k=k)
          se[k]<-suppressWarnings(sqrt(dN%*%Sigma[parindex[[i]]+1:ncol(DM[[i]]),parindex[[i]]+1:ncol(DM[[i]])]%*%dN))
          lower[k]<-m$mle[[i]][-1,k]-se[k]*qnorm(1-(1-alpha)/2)
          upper[k]<-m$mle[[i]][-1,k]+se[k]*qnorm(1-(1-alpha)/2)
        }
        lower<-matrix(c(rep(NA,nbStates),lower),ncol=nbStates,byrow=T)
        upper<-matrix(c(rep(NA,nbStates),upper),ncol=nbStates,byrow=T)  
        se<-matrix(c(rep(NA,nbStates),se),ncol=nbStates,byrow=T)
        Par[[i]] <- parm_list(se,lower,upper,m$mle[[i]])
      }
    }
    
    #if(check)
    #  warning(paste("Some of the parameter estimates seem to lie close to the boundaries of",
    #                "their parameter space.\n  The associated CIs are probably unreliable",
    #                "(or might not be computable)."))
  }

  if(nbStates>1 & !is.null(m$mle$gamma)) {
    # identify parameters of interest
    i2 <- tail(cumsum(unlist(lapply(DM,ncol))),1)+1
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
    Par$gamma <- list(se=se,lower=lower,upper=upper)
    rownames(Par$gamma$se) <- m$stateNames
    rownames(Par$gamma$lower) <- m$stateNames
    rownames(Par$gamma$upper) <- m$stateNames
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
    lower[i]<-m$mle$delta[i]-quantSup*se[i]
    upper[i]<-m$mle$delta[i]+quantSup*se[i]
  }
  lower<-matrix(lower,nrow=1,ncol=nbStates,byrow=T)
  upper<-matrix(upper,nrow=1,ncol=nbStates,byrow=T)  
  se<-matrix(se,nrow=1,ncol=nbStates,byrow=T)
  Par$delta <- list(se=se,lower=lower,upper=upper)    

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

parm_list<-function(se,lower,upper,m){
  Par <- list(se=se,lower=lower,upper=upper)
  rownames(Par$se) <- rownames(m)
  rownames(Par$lower) <- rownames(m)
  rownames(Par$upper) <- rownames(m)
  colnames(Par$se) <- colnames(m)
  colnames(Par$lower) <- colnames(m)
  colnames(Par$upper) <- colnames(m)
  Par
}

get_CI<-function(Par,m,ind,DM,Bounds,cons,boundInd,logitcons,Sigma,nbStates,alpha,mmlePar){

  w<-m$mod$estimate[ind]
  lower<-upper<-se<-numeric(nrow(DM))
  for(k in 1:nrow(DM)){
    dN<-numDeriv::grad(w2nDM,w,bounds=Bounds,DM=DM,cons=cons,boundInd=boundInd,logitcons=logitcons,k=k)
    se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
    cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
    lower[k]<-Par[k]/cn
    upper[k]<-Par[k]*cn
  }
  l<-matrix(lower,ncol=nbStates,byrow=T)
  u<-matrix(upper,ncol=nbStates,byrow=T)  
  s<-matrix(se,ncol=nbStates,byrow=T)
  out <- parm_list(s,l,u,mmlePar)
  out
}