
#' Scaling function: natural to working parameters.
#'
#' Scales each parameter from its natural interval to the set of real numbers, to allow for
#' unconstrained optimization. Used during the optimization of the log-likelihood. Parameters of
#' any data streams for which a design matrix is specified are ignored.
#'
#' @param par Named list of vectors containing the initial parameter values for each data stream.
#' @param bounds Named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream.
#' @param beta Matrix of regression coefficients for the transition probabilities.
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
#' wpar <- momentuHMM:::n2w(par,bounds,beta,delta,nbStates,m$conditions$estAngleMean,NULL,
#' m$conditions$cons,m$conditions$workcons,m$conditions$Bndind)
#' 
#' #natural parameter
#' p <-   momentuHMM:::w2n(wpar,bounds,parSize,nbStates,nbCovs,m$conditions$estAngleMean,
#' m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,m$conditions$fullDM,
#' m$conditions$DMind,m$conditions$workcons,1,m$conditions$dist,m$conditions$Bndind)
#' }
#'
#' @importFrom boot logit

n2w <- function(par,bounds,beta,delta=NULL,nbStates,estAngleMean,DM,cons,workcons,Bndind)
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
        
        p<-n2wDM(bounds[[i]],diag(length(par[[i]])),par[[i]],cons[[i]],workcons[[i]])
        
        foo <- length(p) - nbStates + 1
        angleMean <- p[(foo - nbStates):(foo - 1)]
        kappa <- p[foo:length(p)]
        x <- kappa * cos(angleMean)
        y <- kappa * sin(angleMean)
        p[(foo - nbStates):(foo - 1)] <- x
        p[foo:length(p)] <- y
      }
      else p<-n2wDM(bounds[[i]],diag(length(par[[i]])),par[[i]],cons[[i]],workcons[[i]])
    }
    wpar <- c(wpar,p)
  }

  wbeta <- as.vector(beta) # if beta is NULL, wbeta is NULL as well
  wdelta <- log(delta[-1]/delta[1]) # if delta is NULL, wdelta is NULL as well
  return(c(wpar,wbeta,wdelta))
}

n2wDM<-function(bounds,DM,par,cons,workcons){
  
  a<-bounds[,1]
  b<-bounds[,2]
  
  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  p<-numeric(nrow(DM))
  
  if(length(ind1)) p[ind1] <- (tan(par[ind1]/2)-workcons[ind1])^(1/cons[ind1])
  
  p[ind2] <- par[ind2]
  
  ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
  ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
  ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
  
  p[ind21]<-(log(par[ind21]-a[ind21])-workcons[ind21])^(1/cons[ind21])
  p[ind22]<-(logit((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))-workcons[ind22])^(1/cons[ind22])
  p[ind23]<-(-log(-par[ind23]+b[ind23])-workcons[ind23])^(1/cons[ind23])
  
  p
}     
