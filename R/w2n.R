
#' Scaling function: working to natural parameters
#'
#' Scales each parameter from the set of real numbers, back to its natural interval.
#' Used during the optimization of the log-likelihood.
#'
#' @param wpar Vector of working parameters.
#' @param bounds Named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream.
#' @param parSize Named list indicating the number of natural parameters of the data stream probability distributions
#' @param nbStates The number of states of the HMM.
#' @param nbCovs The number of covariates.
#' @param estAngleMean Named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy').
#' @param circularAngleMean Named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  
#' @param stationary \code{FALSE} if there are covariates. If TRUE, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param cons Named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. 
#' @param fullDM Named list containing the full (i.e. not shorthand) design matrix for each data stream.
#' @param DMind Named list indicating whether \code{fullDM} includes individual- and/or temporal-covariates for each data stream
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workcons Named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. 
#' @param nbObs Number of observations in the data.
#' @param dist Named list indicating the probability distributions of the data streams. 
#' @param Bndind Named list indicating whether \code{DM} is NULL with default parameter bounds for each data stream.
#'
#' @return A list of:
#' \item{...}{Matrices containing the natural parameters for each data stream (e.g., 'step', 'angle', etc.)}
#' \item{beta}{Matrix of regression coefficients of the transition probabilities}
#' \item{delta}{Initial distribution}
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
#'
#' @importFrom boot inv.logit
#' @importFrom Brobdingnag as.brob sum

w2n <- function(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean,circularAngleMean,stationary,cons,fullDM,DMind,workcons,nbObs,dist,Bndind)
{

  # identify initial distribution parameters
  if(!stationary & nbStates>1) {
    foo <- length(wpar)-nbStates+2
    delta <- wpar[foo:length(wpar)]
    expdelta <- exp(c(0,delta))
    if(!is.finite(sum(expdelta))){
      delta <- exp(Brobdingnag::as.brob(c(0,delta)))
      delta <- as.numeric(delta/Brobdingnag::sum(delta))
    } else {
      delta <- expdelta/sum(expdelta)
    }
    wpar <- wpar[-(foo:length(wpar))]
  }
  else delta <- NULL

  # identify regression coefficients for the transition probabilities
  if(nbStates>1) {
    foo <- length(wpar)-(nbCovs+1)*nbStates*(nbStates-1)+1
    beta <- wpar[foo:length(wpar)]
    beta <- matrix(beta,nrow=nbCovs+1)
    wpar <- wpar[-(foo:length(wpar))]
  }
  else beta <- NULL
  
  distnames <- names(dist)
  parindex <- c(0,cumsum(unlist(lapply(fullDM,ncol)))[-length(fullDM)])
  names(parindex) <- names(fullDM)

  parlist<-list()
  
  for(i in distnames){
    tmpwpar<-wpar[parindex[[i]]+1:ncol(fullDM[[i]])]
    if(estAngleMean[[i]] & Bndind[[i]]){ 
      bounds[[i]][,1] <- -Inf
      bounds[[i]][which(bounds[[i]][,2]!=1),2] <- Inf

      foo <- length(tmpwpar) - nbStates + 1
      x <- tmpwpar[(foo - nbStates):(foo - 1)]
      y <- tmpwpar[foo:length(tmpwpar)]
      angleMean <- Arg(x + (0+1i) * y)
      kappa <- sqrt(x^2 + y^2)
      tmpwpar[(foo - nbStates):(foo - 1)] <- angleMean
      tmpwpar[foo:length(tmpwpar)] <- kappa
    }
    parlist[[i]]<-w2nDM(tmpwpar,bounds[[i]],fullDM[[i]],DMind[[i]],cons[[i]],workcons[[i]],nbObs,circularAngleMean[[i]],nbStates)

    if((dist[[i]] %in% angledists) & !estAngleMean[[i]]){
      tmp<-matrix(0,nrow=(parSize[[i]]+1)*nbStates,ncol=nbObs)
      tmp[nbStates+1:nbStates,]<-parlist[[i]]
      parlist[[i]] <- tmp
    }
    
  }

  parlist[["beta"]]<-beta
  parlist[["delta"]]<-delta

  return(parlist)
}

w2nDM<-function(wpar,bounds,DM,DMind,cons,workcons,nbObs,circularAngleMean,nbStates,k=0){
  
  a<-bounds[,1]
  b<-bounds[,2]
  
  piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
  ind1<-which(piInd)
  ind2<-which(!piInd)
  
  XB <- p <- getXB(DM,nbObs,wpar,cons,workcons,DMind,circularAngleMean,nbStates)
  
  if(length(ind1) & !circularAngleMean)
    p[ind1,] <- (2*atan(XB[ind1,]))
  
  ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
  ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
  ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
  
  p[ind21,] <- (exp(XB[ind21,,drop=FALSE])+a[ind21])
  p[ind22,] <- ((b[ind22]-a[ind22])*boot::inv.logit(XB[ind22,,drop=FALSE])+a[ind22])
  p[ind23,] <- -(exp(-XB[ind23,,drop=FALSE]) - b[ind23])
  
  if(any(p<a | p>b))
    stop("Scaling error. Check initial values and bounds.")
  
  if(k) {
    p <- p[k]
  } else if(DMind) {
    p <- matrix(p,length(ind1)+length(ind2),nbObs)
  }
  return(p)
}   
