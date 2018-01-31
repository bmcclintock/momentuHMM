
#' Check parameter length and order for a \code{\link{fitHMM}} (or \code{\link{MIfitHMM}}) model
#' 
#' Prints parameters with labels based on \code{DM}, \code{formula}, and/or \code{formulaDelta}.  See \code{\link{fitHMM}} for 
#' further argument details.
#' 
#' @param data \code{\link{momentuHMMData}} object or a data frame containing the data stream and covariate values
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. 
#' @param Par0 Optional named list containing vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}.  If \code{Par0} is not provided, then ordered parameter indices are returned. 
#' @param beta0 Optional matrix of regression coefficients for the transition probabilities. If \code{beta0} is not provided, then ordered parameter indices are returned. 
#' @param delta0 Optional values or regression coefficients for the initial distribution of the HMM. If \code{delta0} is not provided, then ordered parameter indices are returned. 
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
#' @param cons Deprecated: please use \code{workBounds} instead. An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. 
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. 
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. 
#' @param workcons Deprecated: please use \code{workBounds} instead. An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. 
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. 
#' 
#' @seealso \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#' 
#' @examples
#' m <- example$m
#' checkPar0(data=m$data, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula)
#' 
#' par <- getPar(m)
#' checkPar0(data=m$data, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula,
#'           Par0=par$Par, beta0=par$beta, delta0=par$delta)
#'           
#' dummyDat <- data.frame(step=0,angle=0,cov1=0,cov2=0)
#' checkPar0(data=dummyDat, nbStates=2, dist=m$conditions$dist,
#'           estAngleMean = m$conditions$estAngleMean,
#'           formula = m$conditions$formula)
#'
#' \dontrun{
#' simDat <- simData(nbStates=2, dist=m$conditions$dist, Par = par$Par,
#'                   spatialCovs = list(forest=forest),
#'                   centers = matrix(0,1,2),
#'                   nbCovs = 2)
#' checkPar0(data = simDat, nbStates=2, dist=m$conditions$dist,
#'           formula = ~forest,
#'           DM = list(step=list(mean=~cov1, sd=~cov2),
#'                     angle=list(mean=~center1.angle,concentration=~1)),
#'           estAngleMean=list(angle=TRUE),
#'           circularAngleMean=list(angle=TRUE))
#'           
#' par <- list(step=rnorm(8),angle=rnorm(4))
#' beta0 <- matrix(rnorm(4),2,2)
#' delta0 <- c(0.5,0.5)
#' checkPar0(data = simDat, nbStates=2, dist=m$conditions$dist,
#'           Par0 = par, beta0 = beta0, delta0 = delta0,
#'           formula = ~forest,
#'           DM = list(step=list(mean=~cov1, sd=~cov2),
#'                     angle=list(mean=~center1.angle,concentration=~1)),
#'           estAngleMean=list(angle=TRUE),
#'           circularAngleMean=list(angle=TRUE))                
#' }
#' @export

checkPar0 <- function(data,nbStates,dist,Par0=NULL,beta0=NULL,delta0=NULL,estAngleMean=NULL,circularAngleMean=NULL,formula=~1,formulaDelta=~1,stationary=FALSE,DM=NULL,cons=NULL,userBounds=NULL,workBounds=NULL,workcons=NULL,stateNames=NULL,fixPar=NULL)
{
  
  ## check that the data is a momentuHMMData object or valid data frame
  if(!is.momentuHMMData(data)){
    if(!is.data.frame(data)) stop('data must be a data.frame')
    if(is.null(data$ID)) data$ID <- rep(1,nrow(data))
    data <- momentuHMMData(data)
  }
  
  if(is.null(Par0)){
    
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    distnames<-names(dist)
    if(any(is.na(match(distnames,names(data))))) stop(paste0(distnames[is.na(match(distnames,names(data)))],collapse=", ")," not found in data")
    for(i in distnames){
      dist[[i]]<-match.arg(dist[[i]],momentuHMMdists)
      if(!is.null(fixPar[[i]])) stop("fixPar$",i," cannot be specified unless Par0 is specified")
    }
    
    if(is.null(estAngleMean)){
      estAngleMean <- vector('list',length(distnames))
      names(estAngleMean) <- distnames
    } else {
      if(!is.list(estAngleMean) | is.null(names(estAngleMean))) stop("'estAngleMean' must be a named list")
    }
    for(i in distnames[which(!(dist %in% angledists))]){
      estAngleMean[[i]] <- FALSE
    }
    for(i in distnames){
      if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
      if(!is.logical(estAngleMean[[i]])) stop("estAngleMean$",i," must be logical")
    }
    estAngleMean<-estAngleMean[distnames]
    
    if(is.null(circularAngleMean)){
      circularAngleMean <- vector('list',length(distnames))
      names(circularAngleMean) <- distnames
    } else {
      if(!is.list(circularAngleMean) | is.null(names(circularAngleMean))) stop("'circularAngleMean' must be a named list")
    }
    consensus <- vector('list',length(distnames))
    names(consensus) <- distnames
    
    for(i in distnames){
      if(dist[[i]]=="vmConsensus"){
        consensus[[i]] <- TRUE
        estAngleMean[[i]] <- circularAngleMean[[i]] <- TRUE
        if(is.null(DM[[i]])) stop("DM$",i," must be specified when dist$",i,"=vmConsensus")
        #dist[[i]] <- gsub("Consensus","",dist[[i]])
      } else consensus[[i]] <- FALSE
      if(is.null(circularAngleMean[[i]]) | !estAngleMean[[i]]) circularAngleMean[[i]] <- FALSE
      if(!is.logical(circularAngleMean[[i]])) stop("circularAngleMean$",i," must be logical")
      if(circularAngleMean[[i]] & is.null(DM[[i]])) stop("DM$",i," must be specified when circularAngleMean$",i,"=TRUE")
    }
    circularAngleMean<-circularAngleMean[distnames]
    consensus<-consensus[distnames]
    
    # determine whether zero-inflation or one-inflation should be included
    zeroInflation <- oneInflation <- vector('list',length(distnames))
    names(zeroInflation) <- names(oneInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% zeroInflationdists){
        if(length(which(data[[i]]==0))>0) {
          zeroInflation[[i]]<-TRUE
        }
        else 
          zeroInflation[[i]]<-FALSE
      }
      else zeroInflation[[i]]<-FALSE
      if(dist[[i]] %in% oneInflationdists){
        if(length(which(data[[i]]==1))>0) {
          oneInflation[[i]]<-TRUE
        }
        else 
          oneInflation[[i]]<-FALSE
      }
      else oneInflation[[i]]<-FALSE
    }
    p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds)
    nPar <- lapply(p$bounds,function(x) runif(nrow(x),ifelse(is.finite(x[,1]),x[,1],0),ifelse(is.finite(x[,2]),x[,2],max(x[,1]+1,1.e+6))))
    par <- getParDM(data,nbStates,dist,nPar,zeroInflation,oneInflation,estAngleMean,circularAngleMean,DM,cons,userBounds,workBounds,workcons)
    
  } else par <- Par0
  m<-suppressMessages(fitHMM(data=data,nbStates=nbStates,dist=dist,
                             Par0=par,beta0=beta0,delta0=delta0,
                             estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                             formula=formula,formulaDelta=formulaDelta,stationary=stationary,
                             DM=DM,cons=cons,userBounds=userBounds,workBounds=workBounds,workcons=workcons,fit=FALSE,
                             stateNames=stateNames,fixPar=fixPar))
  
  inputs <- checkInputs(nbStates,dist,par,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$oneInflation,DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,stateNames,checkInflation = TRUE)
  p <- inputs$p
  DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,par,inputs$cons,inputs$workcons,m$conditions$zeroInflation,m$conditions$oneInflation,inputs$circularAngleMean)
  
  beta <- beta0
  if(is.null(beta0) & nbStates>1) {
    beta <- m$mle$beta
    if(!is.null(fixPar$beta)) stop("fixPar$beta cannot be specified unless beta0 is specified")
  }
  
  delta <- delta0
  nbCovsDelta <- ncol(m$covsDelta)-1
  if(stationary)
    delta <- NULL
  else if(!nbCovsDelta){
    if(is.null(delta0)){
      if(!is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
      delta <- matrix(1/nbStates,nbCovsDelta+1,nbStates)
    } else {
      delta <- matrix(delta,nbCovsDelta+1,nbStates)
    }
    delta <- log(delta[-1]/delta[1])
  } else if(is.null(delta0)){
    if(!is.null(fixPar$delta)) stop("fixPar$delta cannot be specified unless delta0 is specified")
    delta <- matrix(0,nrow=nbCovsDelta+1,ncol=nbStates-1)
  }
  
  wpar <- n2w(par,p$bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  
  m$mod <- list()
  m$mod$estimate <- wpar
  m$mod$hessian <- matrix(0,length(wpar),length(wpar))
  
  m$CIreal<-CIreal(m)
  m$CIbeta<-CIbeta(m)
  
  distnames <- names(inputs$dist)
  
  for(i in distnames){
    cat("\n")
    if(is.null(inputs$DM[[i]])) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      tmpPar <- m$mle[[i]]
      if((inputs$dist[[i]] %in% angledists) & !m$conditions$estAngleMean[[i]])
        tmpPar <- tmpPar[-1,,drop=FALSE]
      if(is.null(Par0))
        tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),byrow=TRUE,dimnames=list(rownames(tmpPar),colnames(tmpPar)))
      print(tmpPar)
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      tmpPar <- m$CIbeta[[i]]$est
      if(is.null(Par0))
        tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
      print(tmpPar)
    }
  }
  
  if(length(m$stateNames)>1){
    #if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("---------------------------------------------------\n")
    tmpPar <- m$mle$beta
    if(is.null(beta0))
      tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
    print(tmpPar)
    #}
    
    if(!stationary){
      cat("\n")
      m <- delta_bc(m)
      if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
        tmpPar <- m$mle$delta[1,]
        rownames(tmpPar)<-NULL
        if(is.null(delta0))
          tmpPar[1:nbStates] <- 1:length(tmpPar)
        cat("Initial distribution:\n")
        cat("---------------------\n")
        print(tmpPar)
      } else {
        cat("Regression coeffs for the initial distribution:\n")
        cat("---------------------------------------------------\n")
        tmpPar <- m$CIbeta$delta$est
        if(is.null(delta0))
          tmpPar <- matrix(1:length(tmpPar),nrow(tmpPar),ncol(tmpPar),dimnames=list(rownames(tmpPar),colnames(tmpPar)))
        print(tmpPar)
      }
    }
  }
}