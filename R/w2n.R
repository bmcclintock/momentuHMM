
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

w2n <- function(wpar,stepDist,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,boundInd,logitcons)
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
  
  #tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
  #tmpbounds<-bounds
  
  #if(!is.numeric(bounds)){
 #   bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
  #}
  
  #for(i in 1:nbPar) {
    #index <- (i-1)*nbStates+1
  switch(stepDist,
    "gamma"={wstep<-which(substr(rownames(bounds),1,6)=="gammab")},
    "weibull"={wstep<-which(substr(rownames(bounds),1,9)=="weibullsc")}
  )
  wangle <-which(substr(rownames(bounds),1,2)=="co")
  womega <-which(substr(rownames(bounds),1,2)=="om")
  wdry  <-which(substr(rownames(bounds),1,2)=="dr") 
  wdive <-which(substr(rownames(bounds),1,2)=="di") 
  wice  <-which(substr(rownames(bounds),1,2)=="ic") 
  wland <-which(substr(rownames(bounds),1,2)=="la") 

#nbPar <- length(wpar)/nbStates
  lwpar <- (parSize>0)*c(ncol(stepDM),ncol(angleDM),ncol(omegaDM),ncol(dryDM),ncol(diveDM),ncol(iceDM),ncol(landDM))
  lwparSize <- sum(lwpar)
  lwpar <- lwpar[which(lwpar>0)]
  lwpar <- 1+c(0,cumsum(lwpar[-length(lwpar)]))

  if(lwparSize!=length(wpar))
    stop("Wrong number of parameters.")

  par <- NULL
  nbounds<-NULL
  
  for(index in lwpar) {
    #a <- bounds[index,1]
    #b <- bounds[index,2]
    #p <- wpar[index]
    
    #if(is.finite(a) & is.finite(b)) { # R -> [a,b]
    #  p <- sum((b-a)*inv.logit(p)+a)
    #}
    #else if(is.infinite(a) & is.finite(b)) { # R -> ]-Inf,b]
    #  p <- sum(-(exp(-p)-b))
    #}
    #else if(is.finite(a) & is.infinite(b)) { # R -> [a,Inf[
    #  p <- sum(exp(p)+a)
    #}
    #else if(index %in% wstep){
    if(index %in% wstep){
      Par<-w2nDM(wpar[wstep],bounds[wstep,],stepDM,cons$step,boundInd$step)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wstep,],function(x) eval(parse(text=x))),ncol=2)[boundInd$step,])
    }
    else if(index %in% wangle){
      Par<-w2nDM(wpar[wangle],bounds[wangle,],angleDM,cons$angle,boundInd$angle,logitcons)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wangle,],function(x) eval(parse(text=x))),ncol=2)[boundInd$angle,])
    }
    else if(index %in% womega){
      Par<-w2nDM(wpar[womega],bounds[womega,],omegaDM,cons$omega,boundInd$omega)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[womega,],function(x) eval(parse(text=x))),ncol=2)[boundInd$omega,])
    }
    else if(index %in% wdry){
      Par<-w2nDM(wpar[wdry],bounds[wdry,],dryDM,cons$dry,boundInd$dry)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wdry,],function(x) eval(parse(text=x))),ncol=2)[boundInd$dry,])
    }
    else if(index %in% wdive){
      Par<-w2nDM(wpar[wdive],bounds[wdive,],diveDM,cons$dive,boundInd$dive)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wdive,],function(x) eval(parse(text=x))),ncol=2)[boundInd$dive,])
    }
    else if(index %in% wice){
      Par<-w2nDM(wpar[wice],bounds[wice,],iceDM,cons$ice,boundInd$ice)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wice,],function(x) eval(parse(text=x))),ncol=2)[boundInd$ice,])
    }
    else if(index %in% wland){
      Par<-w2nDM(wpar[wland],bounds[wland,],landDM,cons$land,boundInd$land)
      nbounds <- rbind(nbounds, matrix(sapply(bounds[wland,],function(x) eval(parse(text=x))),ncol=2)[boundInd$land,])
    }
    par <- c(par,Par)
  }
    
  if(length(which(par<nbounds[,1] | par>nbounds[,2]))>0)
    stop("Scaling error.",rownames(nbounds)[which(par<nbounds[,1] | par>nbounds[,2])],par[which(par<nbounds[,1] | par>nbounds[,2])],nbounds[which(par<nbounds[,1] | par>nbounds[,2]),])

  # identify parameters related to step dist
  if(parSize[1])
    stepPar <- matrix(par[1:ncol(stepDM)],ncol=nbStates,byrow=T)
  
  # identify parameters related to angle dist
  if(parSize[2]) {
    anglePar <- matrix(par[parSize[1]*nbStates+1:(parSize[2]*nbStates)],ncol=nbStates,byrow=T)
    par <- par[-(parSize[1]*nbStates+1:(parSize[2]*nbStates))] # remove pars related to angle dist
  }
  else anglePar <- NULL
  
  omegaPar <- dryPar <- divePar <-  icePar <- landPar <-  NULL
  
  if(parSize[3]) omegaPar<- matrix(par[parSize[1]*nbStates+1:(parSize[3]*nbStates)],ncol=nbStates,byrow=T)
  if(parSize[4]) dryPar <-  matrix(par[parSize[1]*nbStates+parSize[3]*nbStates+1:(parSize[4]*nbStates)],ncol=nbStates,byrow=T) 
  if(parSize[5]) divePar <-  matrix(par[parSize[1]*nbStates+parSize[3]*nbStates+parSize[4]*nbStates+1:(parSize[5]*nbStates)],ncol=nbStates,byrow=T) 
  if(parSize[6]) icePar <-  matrix(par[parSize[1]*nbStates+parSize[3]*nbStates+parSize[4]*nbStates+parSize[5]*nbStates+1:(parSize[6]*nbStates)],ncol=nbStates,byrow=T) 
  if(parSize[7]) landPar <-  matrix(par[parSize[1]*nbStates+parSize[3]*nbStates+parSize[4]*nbStates+parSize[5]*nbStates+parSize[6]*nbStates+1:(parSize[7]*nbStates)],ncol=nbStates,byrow=T) 
  
  return(list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar,beta=beta,delta=delta))
}

