
#' Scaling function: natural to working parameters.
#'
#' Scales each parameter from its natural interval to the set of real numbers, to allow for
#' unconstrained optimization. Used during the optimization of the log-likelihood.
#'
#' @param par Vector of state-dependent distributions parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in \code{par}. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param beta Matrix of regression coefficients for the transition probabilities.
#' @param delta Initial distribution. Default: \code{NULL} ; if the initial distribution is not estimated.
#' @param nbStates The number of states of the HMM.
#' @param estAngleMean \code{TRUE} if the angle mean is estimated, \code{FALSE} otherwise.
#' @param cons Logical indicating whether to use parameter constraints (cons=\code{1}) or not (cons=\code{0}).
#' 
#' @return A vector of unconstrained parameters.
#'
#' @examples
#' \dontrun{
#' nbStates <- 3
#' par <- c(0.001,0.999,0.5,0.001,1500.3,7.1)
#' bounds <- matrix(c(0,1, # bounds for first parameter
#'                    0,1, # bounds for second parameter
#'                    0,1, # ...
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf),
#'                  byrow=TRUE,ncol=2)
#' beta <- matrix(rnorm(18),ncol=6,nrow=3)
#' delta <- c(0.6,0.3,0.1)
#'
#' # vector of working parameters
#' wpar <- n2w(par=par,bounds=bounds,beta=beta,delta=delta,nbStates=nbStates,
#'            estAngleMean=FALSE)
#' }
#'
#' @importFrom boot logit

n2w <- function(par,stepDist,bounds,beta,delta=NULL,nbStates,estAngleMean,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,cons,logitcons)
{
  #nbPar <- length(par)/nbStates
  wpar <- NULL
  #for(i in 1:nbPar) {
  #index <- (i-1)*nbStates+1
  switch(stepDist,
    "gamma"={wstep<-which(substr(rownames(bounds),1,6)=="gammab")},
    "weibull"={wstep<-which(substr(rownames(bounds),1,9)=="weibullsc")}
  )
  wangle <-which(substr(rownames(bounds),1,2)=="co")
  womega <-which(substr(rownames(bounds),1,2)=="om")
  wdry <- which(substr(rownames(bounds),1,2)=="dr")
  wdive <- which(substr(rownames(bounds),1,2)=="di")
  wice <- which(substr(rownames(bounds),1,2)=="ic")
  wland <- which(substr(rownames(bounds),1,2)=="la")
  remParms <- NULL
  lwpar<- (1:length(par))
  if(length(wstep)){
    remParms <- c(remParms,wstep[(2:length(wstep))])      
  }
  if(length(wangle)){
    remParms <- c(remParms,wangle[(2:length(wangle))])      
  }
  if(length(womega)){
    remParms <- c(remParms,womega[(2:length(womega))])
  }
  if(length(wdry)){
    remParms <- c(remParms,wdry[(2:length(wdry))])      
  }
  if(length(wdive)){
    remParms <- c(remParms,wdive[(2:length(wdive))])      
  }
  if(length(wice)){
    remParms <- c(remParms,wice[(2:length(wice))])      
  }
  if(length(wland)){
    remParms <- c(remParms,wland[(2:length(wland))])      
  }
  if(length(remParms)){
    lwpar<- (1:length(par))[-remParms]
  }
  
  #tmpind<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)!=bounds
  #tmpbounds<-bounds
  
  #if(!is.numeric(bounds)){
  #  bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(bounds)))
  #}

  #par <- par - bounds[,1]*tmpind[,1]
  
  for(index in lwpar) {
    #a <- bounds[index,1]
    #b <- bounds[index,2]
    #p <- par[index]
    
    #if(is.finite(a) & is.finite(b)) { # [a,b] -> R
    #  p <- logit((p-a)/(b-a))
    #}
    #else if(is.infinite(a) & is.finite(b)) { # ]-Inf,b] -> R
    #  p <- -log(-p+b)
    #}
    #else if(is.finite(a) & is.infinite(b)) { # [a,Inf[ -> R
    #  p <- log(p-a)
    #} 
    if(index==wstep[1]){
      p<-n2wDM(bounds[wstep,],stepDM,par[wstep],cons$step)
    }
    else if(length(wangle) & index==wangle[1]){
      p<-n2wDM(bounds[wangle,],angleDM,par[wangle],cons$angle,logitcons)
    }
    else if(length(womega) & index==womega[1]){
      p<-n2wDM(bounds[womega,],omegaDM,par[womega],cons$omega,)
    }
    else if(length(wdry) & index==wdry[1]){
      p<-n2wDM(bounds[wdry,],dryDM,par[wdry],cons$dry)
    }
    else if(length(wdive) & index==wdive[1]){
      p<-n2wDM(bounds[wdive,],diveDM,par[wdive],cons$dive)
    }
    else if(length(wice) & index==wice[1]){
      p<-n2wDM(bounds[wice,],iceDM,par[wice],cons$ice)
    }
    else if(length(wland) & index==wland[1]){
      p<-n2wDM(bounds[wland,],landDM,par[wland],cons$land)
    }
    wpar <- c(wpar,p)
  }

  if(estAngleMean) {
    # identify angle distribution parameters
    foo <- length(wpar)-nbStates+1
    angleMean <- wpar[(foo-nbStates):(foo-1)]
    kappa <- wpar[foo:length(wpar)]

    # compute the working parameters for the angle distribution
    x <- kappa*cos(angleMean)
    y <- kappa*sin(angleMean)

    wpar[(foo-nbStates):(foo-1)] <- x
    wpar[foo:length(wpar)] <- y
  }

  wbeta <- as.vector(beta) # if beta is NULL, wbeta is NULL as well
  wdelta <- log(delta[-1]/delta[1]) # if delta is NULL, wdelta is NULL as well
  return(c(wpar,wbeta,wdelta))
}
