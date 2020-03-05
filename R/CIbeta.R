
#' Confidence intervals for working (i.e., beta) parameters
#'
#' Computes the standard errors and confidence intervals on the beta (i.e., working) scale of the data stream probability distribution parameters,
#' as well as for the transition probabilities regression parameters. Working scale depends on the real (i.e., natural) scale of the parameters. For 
#' non-circular distributions or for circular distributions with \code{estAngleMean}=FALSE:
#' 
#' 1) if both lower and upper bounds are finite then logit is the working scale;
#' 2) if lower bound is finite and upper bound is infinite then log is the working scale.
#'
#' For circular distributions with \code{estAngleMean}=TRUE and no constraints imposed by a design matrix (DM) or bounds (userBounds), then the working parameters 
#' are complex functions of both the angle mean and concentrations/sd natural parameters (in this case, it's probably best just to focus on the real parameter
#' estimates!).  However, if constraints are imposed by DM or userBounds on circular distribution parameters with \code{estAngleMean}=TRUE and \code{circularAngleMean}=FALSE:
#' 
#' 1) if the natural bounds are (-pi,pi] then tangent is the working scale, otherwise if both lower and upper bounds are finite then logit is the working scale;
#' 2) if lower bound is finite and upper bound is infinite then log is the working scale.
#' 
#' When circular-circular regression is specified using \code{circularAngleMean}, the working scale 
#' for the mean turning angle is not as easily interpretable, but the 
#' link function is atan2(sin(X)*B,1+cos(X)*B), where X are the angle covariates and B the angle coefficients. 
#' Under this formulation, the reference turning angle is 0 (i.e., movement in the same direction as the previous time step). 
#' In other words, the mean turning angle is zero when the coefficient(s) B=0.
#' 
#' @param m A \code{momentuHMM} object
#' @param alpha Significance level of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#'
#' @return A list of the following objects:
#' \item{...}{List(s) of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the working parameters of the data streams}
#' \item{beta}{List of estimates ('est'), standard errors ('se'), and confidence intervals ('lower', 'upper') for the working parameters of the transition probabilities}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CIbeta(m)
#'
#' @export
#' @importFrom utils tail

CIbeta <- function(m,alpha=0.95)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(!is.null(m$conditions$fit) && !m$conditions$fit)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- length(m$stateNames)
  
  dist <- m$conditions$dist
  distnames <- names(dist)
  fullDM <- m$conditions$fullDM
  
  m <- delta_bc(m)

  # identify covariates
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,m$data,par=m$mle)
  recharge <- reForm$recharge
  covs <- reForm$covs
  nbCovs <- reForm$nbCovs

  # inverse of Hessian
  if(!is.null(m$mod$hessian) && !inherits(m$mod$Sigma,"error")) Sigma <- m$mod$Sigma
  else Sigma <- NULL

  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,m$conditions$DM,m$conditions$userBounds)
  bounds <- p$bounds
  
  ncmean <- get_ncmean(distnames,fullDM,m$conditions$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  tmPar <- lapply(m$mle[distnames],function(x) c(t(x)))
  parCount<- lapply(fullDM,ncol)
  for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount))[-length(fullDM)])
  names(parindex) <- distnames
  
  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)
  
  wpar <- m$mod$estimate
  
  Par <- list()
  for(i in distnames){
    est <- w2wn(wpar[parindex[[i]]+1:parCount[[i]]]^m$conditions$cons[[i]]+m$conditions$workcons[[i]],m$conditions$workBounds[[i]])
    
    pnames <- colnames(fullDM[[i]])
    if(!isFALSE(m$conditions$circularAngleMean[[i]])) pnames <- unique(gsub("cos","",gsub("sin","",pnames)))
    Par[[i]] <- get_CIwb(wpar[parindex[[i]]+1:parCount[[i]]],est,parindex[[i]]+1:parCount[[i]],Sigma,alpha,m$conditions$workBounds[[i]],cnames=pnames,cons=m$conditions$cons[[i]])

  }
  
  mixtures <- m$conditions$mixtures

  # group CIs for t.p. coefficients
  if(nbStates>1){
    
    est <- w2wn(wpar[tail(cumsum(unlist(parCount)),1)+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures)],m$conditions$workBounds$beta)
    
    Par$beta <- get_CIwb(wpar[tail(cumsum(unlist(parCount)),1)+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures)],est,tail(cumsum(unlist(parCount)),1)+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures),Sigma,alpha,m$conditions$workBounds$beta,rnames=rownames(m$mle$beta),cnames=colnames(m$mle$beta),cons=rep(1,length(est)))
    
    # fill in constraints based on betaCons
    Par$beta <- lapply(Par$beta,function(x) matrix(x[c(m$conditions$betaCons)],dim(x),dimnames=list(rownames(x),colnames(x))))
    
    nbCovsDelta <- ncol(m$covsDelta)-1
    
  }
  
  # group CIs for pi
  if(mixtures>1){
    
    nbCovsPi <- ncol(m$covsPi)-1
    piInd <- tail(cumsum(unlist(parCount)),1) + ((nbCovs+1)*nbStates*(nbStates-1)*mixtures) + 1:((nbCovsPi+1)*(mixtures-1))
    
    est <- w2wn(wpar[piInd],m$conditions$workBounds$pi)
    
    Par$pi <- get_CIwb(wpar[piInd],est,piInd,Sigma,alpha,m$conditions$workBounds$pi,rnames=colnames(m$covsPi),cnames=paste0("mix",2:mixtures),cons=rep(1,length(est)))
  }
  
  # group CIs for initial distribution
  if(nbStates>1 & !m$conditions$stationary){

    dInd <- length(wpar)-ifelse(reForm$nbRecovs,(reForm$nbRecovs+1)+(reForm$nbG0covs+1),0)
    foo <- dInd -(nbCovsDelta+1)*(nbStates-1)*mixtures+1
    
    est <- w2wn(wpar[foo:dInd],m$conditions$workBounds$delta)
    
    rnames <- rep(colnames(m$covsDelta),mixtures)
    if(mixtures>1) rnames <- paste0(rep(colnames(m$covsDelta),mixtures),"_mix",rep(1:mixtures,each=length(colnames(m$covsDelta))))
    Par$delta <- get_CIwb(wpar[foo:dInd],est,foo:dInd,Sigma,alpha,m$conditions$workBounds$delta,rnames=rnames,cnames=m$stateNames[-1],cons=rep(1,length(est)))
    
  }
  
  if(!is.null(recharge)){
    
    ind <- tail(cumsum(unlist(parCount)),1)+(nbCovs+1)*nbStates*(nbStates-1)+(nbCovsDelta+1)*(nbStates-1)+1:(reForm$nbG0covs+1)
    est <- w2wn(wpar[ind],m$conditions$workBounds$g0)
    
    Par$g0 <- get_CIwb(wpar[ind],est,ind,Sigma,alpha,m$conditions$workBounds$g0,rnames="[1,]",cnames=colnames(reForm$g0covs),cons=rep(1,length(est)))
    
    ind <- tail(cumsum(unlist(parCount)),1)+(nbCovs+1)*nbStates*(nbStates-1)+(nbCovsDelta+1)*(nbStates-1)+reForm$nbG0covs+1+1:(reForm$nbRecovs+1)
    est <- w2wn(wpar[ind],m$conditions$workBounds$theta)
    
    Par$theta <- get_CIwb(wpar[ind],est,ind,Sigma,alpha,m$conditions$workBounds$theta,rnames="[1,]",cnames=colnames(reForm$recovs),cons=rep(1,length(est)))
    
  }
  return(Par)
}

get_gradwb<-function(wpar,workBounds,cons){
  ind1<-which(is.finite(workBounds[,1]) & is.infinite(workBounds[,2]))
  ind2<-which(is.finite(workBounds[,1]) & is.finite(workBounds[,2]))
  ind3<-which(is.infinite(workBounds[,1]) & is.finite(workBounds[,2]))
  
  dN <- diag(cons*(wpar^(cons-1)))
  dN[cbind(ind1,ind1)] <- exp(wpar[ind1])
  dN[cbind(ind2,ind2)] <- (workBounds[ind2,2]-workBounds[ind2,1])*exp(wpar[ind2])/(1+exp(wpar[ind2]))^2 
  dN[cbind(ind3,ind3)] <- exp(-wpar[ind3])
  dN
}

get_CIwb<-function(wpar,Par,ind,Sigma,alpha,workBounds,rnames="[1,]",cnames,cons){
  
  npar <- length(wpar)
  bRow <- (rnames=="[1,]")
  lower<-upper<-se<-rep(NA,npar)
  
  if(!is.null(Sigma)){
    dN <- get_gradwb(wpar,workBounds,cons)
    
    se <- suppressWarnings(sqrt(diag(dN%*%Sigma[ind,ind]%*%t(dN))))
    lower <- Par - qnorm(1-(1-alpha)/2) * se
    upper <- Par + qnorm(1-(1-alpha)/2) * se
  }  
  
  est<-matrix(Par,ncol=length(cnames),byrow=bRow)
  l<-matrix(lower,ncol=length(cnames),byrow=bRow)
  u<-matrix(upper,ncol=length(cnames),byrow=bRow)  
  s<-matrix(se,ncol=length(cnames),byrow=bRow)
  
  beta_parm_list(est,s,l,u,rnames,cnames)
}

beta_parm_list<-function(est,se,lower,upper,rnames,cnames){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  rownames(Par$est) <- rownames(Par$se) <- rownames(Par$lower) <- rownames(Par$upper) <- rnames
  colnames(Par$est) <- cnames
  colnames(Par$se) <- cnames
  colnames(Par$lower) <- cnames
  colnames(Par$upper) <- cnames
  Par
}
