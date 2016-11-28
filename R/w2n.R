
#' Scaling function: working to natural parameters
#'
#' Scales each parameter from the set of real numbers, back to its natural interval.
#' Used during the optimization of the log-likelihood.
#'
#' @param wpar Vector of state-dependent distributions unconstrained parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in \code{wpar}. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values: number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution.
#' @param nbStates The number of states of the HMM.
#' @param nbCovs The number of covariates.
#' @param estAngleMean \code{TRUE} if the angle mean is estimated, \code{FALSE} otherwise.
#' @param stationary \code{FALSE} if there are covariates. If TRUE, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param cons Logical indicating whether to use parameter constraints (cons=\code{1}) or not (cons=\code{0}).
#'
#' @return A list of:
#' \item{stepPar}{Matrix of natural parameters of the step length distribution}
#' \item{anglePar}{Matrix of natural parameters of the turning angle distribution}
#' \item{beta}{Matrix of regression coefficients of the transition probabilities}
#' \item{delta}{Initial distribution}
#'
#' @examples
#' \dontrun{
#' nbStates <- 3
#' nbCovs <- 2
#' par <- c(0.001,0.999,0.5,0.001,1500.3,7.1)
#' parSize <- c(1,1)
#' bounds <- matrix(c(0,1,0,1,0,1,
#'                    0,Inf,0,Inf,0,Inf),
#'                  byrow=TRUE,ncol=2)
#' beta <- matrix(rnorm(18),ncol=6,nrow=3)
#' delta <- c(0.6,0.3,0.1)
#' wpar <- n2w(par,bounds,beta,delta,nbStates,FALSE)
#' print(w2n(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean=FALSE,stationary=FALSE))
#' }
#'
#'
#' @importFrom boot inv.logit

w2n <- function(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,DM,boundInd,logitcons)
{
  if(nbStates<1)
    stop("Number of states must be 1 at least.")

  # identify initial distribution parameters
  if(!stationary & nbStates>1) {
    foo <- length(wpar)-nbStates+2
    delta <- wpar[foo:length(wpar)]
    delta <- exp(c(0,delta))
    delta <- delta/sum(delta)
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

  if(estAngleMean) {
    # identify working parameters for the angle distribution (x and y)
    foo <- length(wpar)-nbStates+1
    x <- wpar[(foo-nbStates):(foo-1)]
    y <- wpar[foo:length(wpar)]

    # compute natural parameters for the angle distribution
    angleMean <- Arg(x+1i*y)
    kappa <- sqrt(x^2+y^2)
    # to scale them if necessary (see parDef)
    wpar[(foo-nbStates):(foo-1)] <- angleMean
    wpar[foo:length(wpar)] <- kappa
  }
  
  distnames <- names(bounds)
  parindex <- c(0,cumsum(unlist(lapply(bounds,nrow)))[-length(bounds)])
  names(parindex) <- distnames

  par <- NULL
  nbounds<-NULL
  #Par<-NULL
  parlist<-NULL
  
  for(i in distnames){
    Par<-w2nDM(wpar[parindex[[i]]+1:ncol(DM[[i]])],bounds[[i]],DM[[i]],cons[[i]],boundInd[[i]],logitcons[[i]])
    nbounds <- rbind(nbounds, matrix(sapply(bounds[[i]],function(x) eval(parse(text=x))),ncol=2)[boundInd[[i]],])
    par <- c(par,Par)
    Par<-matrix(Par,ncol=nbStates,byrow=TRUE)
    parlist[[i]]<-Par
  }
    
  if(length(which(par<nbounds[,1] | par>nbounds[,2]))>0)
    stop("Scaling error.",rownames(nbounds)[which(par<nbounds[,1] | par>nbounds[,2])],par[which(par<nbounds[,1] | par>nbounds[,2])],nbounds[which(par<nbounds[,1] | par>nbounds[,2]),])
  
  parlist[["beta"]]<-beta
  parlist[["delta"]]<-delta

  return(parlist)
}

