#' Expand vector of free working parameters to vector of all working parameters including any fixed parameters (used in fitHMM.R and nLogLike.R)
#' 
#' @param optPar vector of free working parameters
#' @param optInd indices of constrained parameters
#' @param fixPar Vector of working parameters which are assumed known prior to fitting the model (NA indicates parameters is to be estimated)
#' @param wparIndex Vector of indices for the elements of \code{fixPar} that are not NA
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters.
#' @param deltaCons Matrix of the same dimension as \code{delta0} composed of integers identifying any equality constraints among the initial distribution working scale parameters. 
#' @param nbStates Number of states of the HMM
#' @param nbCovsDelta Number of initial distribution covariates
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param nbCovs Number of t.p.m. covariates
#' @param nbRecovs Number of recharge covariates
#' @param mixtures Number of mixtures for the state transition probabilities
#' @param nbCovsPi Number of mixture probability covariates
#' 
#' @return A vector of all working parameters including any fixed parameters
#' 
#' @examples
#' \dontrun{
#' nbStates <- 2
#' stepDist <- "gamma" # step distribution
#' angleDist <- "vm" # turning angle distribution
#' 
#' # extract data from momentuHMM example
#' data <- example$m$data
#'
#' ### 1. fit the model to the simulated data
#' # define initial values for the parameters
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
#' anglePar <- kappa0 # not estimating angle mean, so not included
#' formula <- ~cov1+cos(cov2)
#' 
#' # constrain cov1 effect to state 1 -> 2 and cov2 effect to state 2 -> 1
#' fixPar <- list(beta=c(NA,NA,0,NA,0,NA))
#'
#' m <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'             Par0=list(step=stepPar,angle=anglePar),formula=formula,fixPar=fixPar)
#'
#' # convert free parameter vector (m$mod$wpar) to full set of working parameters (m$mod$estimate)
#' est <- momentuHMM:::expandPar(m$mod$wpar,m$conditions$optInd,unlist(m$conditions$fixPar),
#'                               m$conditions$wparIndex,m$conditions$betaCons,m$conditions$deltaCons,
#'                               nbStates,
#'                               ncol(m$covsDelta)-1,m$conditions$stationary,nrow(m$mle$beta)-1)
#'
#' all(est==m$mod$estimate)
#' }
expandPar <- function(optPar,optInd,fixPar,wparIndex,betaCons,deltaCons,nbStates,nbCovsDelta,stationary,nbCovs,nbRecovs=0,mixtures=1,nbCovsPi=0){
  if(length(optInd)){
    wpar <- numeric(length(fixPar))
    wpar[-optInd] <- optPar
    if(length(wparIndex)) wpar[wparIndex] <- fixPar[wparIndex]
    if(nbStates>1){
      if(!is.null(betaCons)){
        foo <- length(wpar)-(nbCovsDelta+1)*(nbStates-1)*(!stationary)*mixtures-ifelse(nbRecovs,nbRecovs+2,0)-(nbCovsPi+1)*(mixtures-1)-((nbCovs+1)*nbStates*(nbStates-1)*mixtures-1):0
        wpar[foo] <- wpar[foo][betaCons]
      }
      if(!is.null(deltaCons)){
        foo <- length(wpar)-ifelse(nbRecovs,nbRecovs+2,0)-((nbCovsDelta+1)*(nbStates-1)*mixtures-1):0
        wpar[foo] <- wpar[foo][deltaCons]
      }
    }
  } else {
    wpar <- optPar
  }
  wpar
}