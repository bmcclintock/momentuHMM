
#' Negative log-likelihood function
#'
#' @param wpar Vector of working parameters.
#' @param nbStates Number of states of the HMM.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in \code{wpar}. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values: number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution.
#' @param data An object \code{moveData}.
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are: gamma, weibull, lnorm, exp. Default: gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are: vm, wrpcauchy. Set to \code{"none"} if the angle distribution should
#' not be estimated. Default: vm.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Default: \code{NULL} (the angle mean is estimated).
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default: \code{FALSE}. If \code{TRUE}, initial values for the zero-mass parameters should be
#' included in \code{stepPar0}.
#' @param stationary \code{FALSE} if there are covariates. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#'
#' @return The negative log-likelihood of the parameters given the data.
#'
#' @examples
#' \dontrun{
#' # data is a moveData object (as returned by prepData), automatically loaded with the package
#' data <- example$data
#' simPar <- example$simPar
#' par0 <- example$par0
#'
#' estAngleMean <- is.null(simPar$angleMean)
#' bounds <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
#'                  estAngleMean,simPar$zeroInflation)$bounds
#' parSize <- parDef(simPar$stepDist,simPar$angleDist,simPar$nbStates,
#'                   estAngleMean,simPar$zeroInflation)$parSize
#'
#' par <- c(par0$stepPar0,par0$anglePar0)
#' wpar <- n2w(par,bounds,par0$beta0,par0$delta0,simPar$nbStates,FALSE)
#'
#' l <- nLogLike(wpar=wpar,nbStates=simPar$nbStates,bounds=bounds,parSize=parSize,data=data,
#'              stepDist=simPar$stepDist,angleDist=simPar$angleDist,angleMean=simPar$angleMean,
#'              zeroInflation=simPar$zeroInflation)
#' }
#'
#' @export

nLogLike <- function(wpar,nbStates,formula,bounds,parSize,data,stepDist=c("gamma","weibull","lnorm","exp"),
                     angleDist=c("wrpcauchy","vm","none"),omegaDist=c("beta","none"),dryDist=c("beta","none"),diveDist=c("pois","none"),iceDist=c("beta","none"),landDist=c("beta","none"),
                     angleMean=NULL,zeroInflation=FALSE,
                     stationary=FALSE,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,boundInd,logitcons)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  omegaDist <- match.arg(omegaDist)
  dryDist <- match.arg(dryDist)
  diveDist <- match.arg(diveDist)
  iceDist <- match.arg(iceDist)
  landDist <- match.arg(landDist)
  if(nbStates<1)
    stop("nbStates must be at least 1.")

  covs <- model.matrix(formula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  initParms<-sum((parSize>0)*c(ncol(stepDM),ncol(angleDM),ncol(omegaDM),ncol(dryDM),ncol(diveDM),ncol(iceDM),ncol(landDM)))
  
  if(!stationary & (length(wpar)!=initParms+nbStates*(nbStates-1)*(nbCovs+1)+nbStates-1))
    stop("Wrong number of parameters in wpar.")
  if(stationary & (length(wpar)!=initParms+nbStates*(nbStates-1)*(nbCovs+1)))
    stop("Wrong number of parameters in wpar.")
  if(length(data)<1)
    stop("The data input is empty.")

  if(is.null(data$step))
    stop("Missing field(s) in data.")

  estAngleMean <- (is.null(angleMean) & angleDist!="none")

  # convert the parameters back to their natural scale
  par <- w2n(wpar,stepDist,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,boundInd,logitcons)
  
  if(!is.null(angleMean) & angleDist!="none") # if the turning angles' mean is not estimated
    par$anglePar <- rbind(angleMean,par$anglePar)

  nbObs <- length(data$step)

  nbAnimals <- length(unique(data$ID))

  # aInd = list of indices of first observation for each animal
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

  # NULL arguments don't suit C++
  if(angleDist=="none")
    par$anglePar <- matrix(NA)
  if(omegaDist=="none")
    par$omegaPar <- matrix(NA)
  if(dryDist=="none")
    par$dryPar <- matrix(NA)
  if(diveDist=="none")
    par$divePar <- matrix(NA)
  if(iceDist=="none")
    par$icePar <- matrix(NA)
  if(landDist=="none")
    par$landPar <- matrix(NA)
  if(stationary)
    par$delta <- c(NA)
  if(nbStates==1) {
    par$beta <- matrix(NA)
    par$delta <- c(NA)
    par$stepPar <- as.matrix(par$stepPar)
    par$anglePar <- as.matrix(par$anglePar)
    par$omegaPar <- as.matrix(par$omegaPar)
    par$dryPar <- as.matrix(par$dryPar)
    par$divePar <- as.matrix(par$divePar)
    par$icePar <- as.matrix(par$icePar)
    par$landPar <- as.matrix(par$landPar)
  }

  nllk <- nLogLike_rcpp(nbStates,par$beta,as.matrix(covs),data,stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,
                        par$stepPar,par$anglePar,par$omegaPar,par$dryPar,par$divePar,par$icePar,par$landPar,
                        par$delta,aInd,zeroInflation,stationary)

  return(nllk)
}
