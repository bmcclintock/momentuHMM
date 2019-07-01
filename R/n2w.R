
#' Scaling function: natural to working parameters.
#'
#' Scales each data stream probability distribution parameter from its natural interval to the set of real numbers, to allow for
#' unconstrained optimization. Used during the optimization of the log-likelihood. Parameters of
#' any data streams for which a design matrix is specified are ignored.
#'
#' @param par Named list of vectors containing the initial parameter values for each data stream.
#' @param bounds Named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream.
#' @param beta List of regression coefficients for the transition probabilities.
#' @param delta Initial distribution. Default: \code{NULL} ; if the initial distribution is not estimated.
#' @param nbStates The number of states of the HMM.
#' @param estAngleMean Named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy').
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a matrix.  
#' @param cons Named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. 
#' @param workcons Named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. 
#' @param Bndind Named list indicating whether \code{DM} is NULL with default parameter bounds for each data stream.
#' @param dist A named list indicating the probability distributions of the data streams.
#' 
#' @return A vector of unconstrained parameters.
#'
#' @examples
#' \dontrun{
#' m<-example$m
#' nbStates <- 2
#' nbCovs <- 2
#' parSize <- list(step=2,angle=2)
#' par <- list(step=c(t(m$mle$step)),angle=c(t(m$mle$angle)))
#' bounds <- m$conditions$bounds
#' beta <- matrix(rnorm(6),ncol=2,nrow=3)
#' delta <- c(0.6,0.4)
#' 
#' #working parameters
#' wpar <- momentuHMM:::n2w(par,bounds,list(beta=beta),log(delta[-1]/delta[1]),nbStates,
#' m$conditions$estAngleMean,NULL,m$conditions$cons,m$conditions$workcons,m$conditions$Bndind,
#' m$conditions$dist)
#' 
#' #natural parameter
#' p <-   momentuHMM:::w2n(wpar,bounds,parSize,nbStates,nbCovs,m$conditions$estAngleMean,
#' m$conditions$circularAngleMean,lapply(m$conditions$dist,function(x) x=="vmConsensus"),
#' m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,
#' m$conditions$DMind,m$conditions$workcons,1,m$conditions$dist,m$conditions$Bndind,
#' matrix(1,nrow=length(unique(m$data$ID)),ncol=1),covsDelta=m$covsDelta,
#' workBounds=m$conditions$workBounds)
#' }
#'
#' @importFrom stats dunif

n2w <- function(par,bounds,beta,delta=NULL,nbStates,estAngleMean,DM,cons,workcons,Bndind,dist)
{
  wpar <- NULL
  for(i in names(par)){
    p <- par[[i]]
    if(is.null(DM[[i]])){
      if(estAngleMean[[i]] & Bndind[[i]]){
        if(length(which(par[[i]][1:nbStates]<=bounds[[i]][1:nbStates,1] | par[[i]][1:nbStates]>bounds[[i]][1:nbStates,2]))>0)
          stop(paste0("Check the parameter bounds for ",i," (the initial parameters should be ",
                      "strictly between the bounds of their parameter space). The angle mean should be in (-pi,pi]."))
        if(length(which(par[[i]][nbStates+1:nbStates]<=bounds[[i]][nbStates+1:nbStates,1] | par[[i]][nbStates+1:nbStates]>=bounds[[i]][nbStates+1:nbStates,2]))>0)
          stop(paste0("Check the parameter bounds for ",i," (the initial parameters should be ",
                      "strictly between the bounds of their parameter space)."))
      } else {
        if(length(which(par[[i]]<=bounds[[i]][,1] | par[[i]]>=bounds[[i]][,2]))>0)
          stop(paste0("Check the parameter bounds for ",i," (the initial parameters should be ",
                      "strictly between the bounds of their parameter space)."))
      }
      if(estAngleMean[[i]] & Bndind[[i]]){ 
        bounds[[i]][,1] <- -Inf
        bounds[[i]][which(bounds[[i]][,2]!=1),2] <- Inf
        
        p<-n2wDM(bounds[[i]],diag(length(par[[i]])),par[[i]],cons[[i]],workcons[[i]],nbStates)
        
        foo <- length(p) - nbStates + 1
        angleMean <- p[(foo - nbStates):(foo - 1)]
        kappa <- p[foo:length(p)]
        x <- kappa * cos(angleMean)
        y <- kappa * sin(angleMean)
        p[(foo - nbStates):(foo - 1)] <- x
        p[foo:length(p)] <- y
      } else if(grepl("cat",dist[[i]])){
        dimCat <- length(par[[i]])/nbStates
        for(j in 1:nbStates){
          p[seq(j,dimCat*nbStates,nbStates)] <- log(par[[i]][seq(j,dimCat*nbStates,nbStates)]/(1-sum(par[[i]][seq(j,dimCat*nbStates,nbStates)])))
        }
      } else p<-n2wDM(bounds[[i]],diag(length(par[[i]])),par[[i]],cons[[i]],workcons[[i]],nbStates)
    }
    wpar <- c(wpar,p)
  }

  wbeta <- as.vector(beta$beta) # if beta is NULL, wbeta is NULL as well
  wpi <- as.vector(beta$pi)
  wg0 <- as.vector(beta$g0)
  wtheta <- as.vector(beta$theta)
  wdelta <- as.vector(delta)#log(delta[,-1,drop=FALSE]/delta[,1,drop=FALSE]) # if delta is NULL, wdelta is NULL as well
  return(c(wpar,wbeta,wpi,wdelta,wg0,wtheta))
}

#' @importFrom stats qlogis
nw2w <-function(wpar,workBounds){
  
  ind1<-which(is.finite(workBounds[,1]) & is.infinite(workBounds[,2]))
  ind2<-which(is.finite(workBounds[,1]) & is.finite(workBounds[,2]))
  ind3<-which(is.infinite(workBounds[,1]) & is.finite(workBounds[,2]))
  
  wpar[ind1] <- log(wpar[ind1]-workBounds[ind1,1])
  wpar[ind2] <- stats::qlogis((wpar[ind2]-workBounds[ind2,1])/(workBounds[ind2,2]-workBounds[ind2,1]))
  wpar[ind3] <- -log(-wpar[ind3] + workBounds[ind3,2])
  
  wpar
  
}

#' @importFrom stats qlogis
n2wDM<-function(bounds,DM,par,cons,workcons,nbStates){
  
  a<-bounds[,1]
  b<-bounds[,2]
  
  zeroInflation <- any(grepl("zeromass",rownames(bounds)))
  oneInflation <- any(grepl("onemass",rownames(bounds)))
  
  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  zoInd <- as.logical((grepl("zeromass",rownames(bounds)) | grepl("onemass",rownames(bounds)))*(zeroInflation*oneInflation))
  ind2<-which(zoInd)
  ind3<-which(!piInd & !zoInd)

  p<-numeric(nrow(DM))
  
  if(length(ind1)) p[ind1] <- (tan(par[ind1]/2)-workcons[ind1])^(1/cons[ind1])
  if(length(ind2)){
    for(j in 1:nbStates){
      zoParInd <- which(grepl(paste0("zeromass_",j),rownames(bounds)) | grepl(paste0("onemass_",j),rownames(bounds)))
      zoPar <- c(par[zoParInd],1.-sum(par[zoParInd]))
      if(any(!stats::dunif(zoPar)) | zoPar[3]==0) stop("the sum of zero-mass and one-mass parameters cannot be >=1")
      p[zoParInd] <- log(zoPar[-3]/zoPar[3])
    }
  }
  
  p[ind3] <- par[ind3]
  
  ind31<-ind3[which(is.finite(a[ind3]) & is.infinite(b[ind3]))]
  ind32<-ind3[which(is.finite(a[ind3]) & is.finite(b[ind3]))]
  ind33<-ind3[which(is.infinite(a[ind3]) & is.finite(b[ind3]))]
  
  p[ind31]<-(log(par[ind31]-a[ind31])-workcons[ind31])^(1/cons[ind31])
  p[ind32]<-(stats::qlogis((par[ind32]-a[ind32])/(b[ind32]-a[ind32]))-workcons[ind32])^(1/cons[ind32])
  p[ind33]<-(-log(-par[ind33]+b[ind33])-workcons[ind33])^(1/cons[ind33])
  
  p
}     
