
#' Check parameter length and order for a \code{\link{fitHMM}} (or \code{\link{MIfitHMM}}) model
#' 
#' Verify initial parameters (\code{Par0}, \code{beta0}, and \code{delta0}) are of correct length and print the parameters with labels based on \code{DM}, \code{formula}, and/or \code{formulaDelta}.  See \code{\link{fitHMM}} for 
#' further argument details.
#' 
#' @param data An object \code{momentuHMMData}.
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. 
#' @param Par0 A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. 
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities. 
#' @param delta0 Initial value for the initial distribution of the HMM. 
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). 
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  
#' @param formula Regression formula for the transition probability covariates. 
#' @param formulaDelta Regression formula for the initial distribution. 
#' @param stationary \code{FALSE} if there are covariates in \code{formula} or \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream.
#' @param cons An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. 
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. 
#' @param workcons An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. 
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. 
#' 
#' @seealso \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#' 
#' @export

checkPar0 <- function(data,nbStates,dist,Par0,beta0=NULL,delta0=NULL,estAngleMean=NULL,circularAngleMean=NULL,formula=~1,formulaDelta=~1,stationary=FALSE,DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL,stateNames=NULL,fixPar=NULL)
{
  
  m<-suppressMessages(fitHMM(data=data,nbStates=nbStates,dist=dist,
                             Par0=Par0,beta0=beta0,delta0=delta0,
                             estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                             formula=formula,formulaDelta=formulaDelta,stationary=stationary,
                             DM=DM,cons=cons,userBounds=userBounds,workcons=workcons,fit=FALSE,
                             stateNames=stateNames,fixPar=fixPar))
  
  inputs <- checkInputs(nbStates,dist,Par0,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,stateNames,checkInflation = TRUE)
  p <- inputs$p
  DMinputs<-getDM(data,inputs$DM,dist,nbStates,p$parNames,p$bounds,Par0,inputs$cons,inputs$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,inputs$circularAngleMean)
  
  if(is.null(beta0) & nbStates>1) {
    beta0 <- m$mle$beta
  }
  
  nbCovsDelta <- ncol(m$covsDelta)-1
  if(stationary)
    delta0 <- NULL
  else if(!nbCovsDelta){
    if(is.null(delta0)){
      delta0 <- matrix(1/nbStates,nbCovsDelta+1,nbStates)
    } else {
      delta0 <- matrix(delta0,nbCovsDelta+1,nbStates)
    }
    delta0 <- log(delta0[-1]/delta0[1])
  } else if(is.null(delta0)){
    delta0 <- matrix(0,nrow=nbCovsDelta+1,ncol=nbStates-1)
  }
  
  wpar <- n2w(Par0,p$bounds,beta0,delta0,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  
  m$mod <- list()
  m$mod$estimate <- wpar
  m$mod$hessian <- matrix(0,length(wpar),length(wpar))
  
  m$CIreal<-CIreal(m)
  m$CIbeta<-CIbeta(m)
  
  distnames <- names(dist)
  
  for(i in distnames){
    cat("\n")
    if(is.null(inputs$DM[[i]])) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      if((dist[[i]] %in% angledists) & !m$conditions$estAngleMean[[i]])
        print(m$mle[[i]][-1,])
      else 
        print(m$mle[[i]])
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      print(m$CIbeta[[i]]$est)
    }
  }
  
  if(length(m$stateNames)>1){
    #if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("---------------------------------------------------\n")
    print(m$mle$beta)
    #}
    
    if(!stationary){
      cat("\n")
      m <- delta_bc(m)
      if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
        tmp <- m$mle$delta[1,]
        rownames(tmp)<-NULL
        cat("Initial distribution:\n")
        cat("---------------------\n")
        print(tmp)
      } else {
        cat("Regression coeffs for the initial distribution:\n")
        cat("---------------------------------------------------\n")
        print(m$CIbeta$delta$est)
      }
    }
  }
}