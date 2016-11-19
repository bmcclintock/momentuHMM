
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

nLogLike <- function(wpar,nbStates,formula,bounds,parSize,data,dist,
                     angleMean=NULL,zeroInflation=FALSE,
                     stationary=FALSE,cons,DM,boundInd,logitcons)
{
  # check arguments
  distnames<-names(dist)
  eval(parse(text=paste0(distnames,"Dist='",dist,"'")))
  
  #####################
  ## Check arguments ##
  #####################
  stepDist <- match.arg(stepDist,c("gamma","weibull","none"))
  angleDist <- match.arg(angleDist,c("wrpcauchy","vm","none"))
  otherDist <- distnames[-which(distnames %in% c("step","angle"))]
  if(length(otherDist)){
    eval(parse(text=paste0(otherDist,"=match.arg(",otherDist,"Dist,c('gamma','weibull','beta','pois'))")))
  }
  if(nbStates<1)
    stop("nbStates must be at least 1.")

  covs <- model.matrix(formula,data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  initParms<-sum((parSize>0)*unlist(lapply(DM,ncol)))
  
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
  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,DM,boundInd,logitcons)
  
  if(!is.null(angleMean) & angleDist!="none") # if the turning angles' mean is not estimated
    par$angle <- rbind(angleMean,par$angle)

  nbObs <- length(data$step)

  nbAnimals <- length(unique(data$ID))

  # aInd = list of indices of first observation for each animal
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

  # NULL arguments don't suit C++
  if(any(unlist(lapply(dist,is.null)))){
    par[which(unlist(lapply(dist,is.null)))]<-matrix(NA)
  }

  if(stationary)
    par$delta <- c(NA)
  if(nbStates==1) {
    par$beta <- matrix(NA)
    par$delta <- c(NA)
    par[paste0(distnames,"Par")] <- lapply(par[paste0(distnames,"Par")],as.matrix)
  }

  nllk <- nLogLike_rcpp(nbStates,as.matrix(covs),data,names(dist),dist,
                        par,
                        aInd,zeroInflation,stationary)

  return(nllk)
}
