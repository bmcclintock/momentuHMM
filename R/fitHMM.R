
#' Fit an HMM to the data
#'
#' Fit an hidden Markov model to the data provided, using numerical optimization of the log-likelihood
#' function.
#'
#' @param data An object \code{moveData}.
#' @param nbStates Number of states of the HMM.
#' @param stepPar Vector of initial state-dependent step length distribution parameters.
#' The parameters should be in the order expected by the pdf of \code{stepDist}, and the zero-mass
#' parameter should be the last. Note that zero-mass parameters are mandatory if there are steps of
#' length zero in the data.
#' For example, for a 2-state model using the gamma distribution and
#' including zero-inflation, the vector of initial parameters would be something like:
#' \code{c(mu1,mu2,sigma1,sigma2,zeromass1,zeromass2)}.
#' @param anglePar Vector of initial state-dependent turning angle distribution parameters.
#' The parameters should be in the order expected by the pdf of \code{angleDist}. For example, for a 2-state
#' model using the Von Mises (vm) distribution, the vector of initial parameters would be something like:
#' \code{c(mu1,mu2,kappa1,kappa2)}.
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities (more
#' information in "Details").
#' Default: \code{NULL}. If not specified, \code{beta0} is initialized such that the diagonal elements
#' of the transition probability matrix are dominant.
#' @param delta0 Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}.
#' @param formula Regression formula for the covariates. Default: \code{~1} (no covariate effect).
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are: gamma, weibull, lnorm, exp. Default: gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are: vm, wrpcauchy. Set to \code{"none"} if the angle distribution should
#' not be estimated. Default: vm.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Default: \code{NULL} (the angle mean is estimated).
#' @param stationary \code{FALSE} if there are covariates. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param verbose Determines the print level of the optimizer. The default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed.
#' @param nlmPar List of parameters to pass to the optimization function \code{nlm} (which should be either
#' '\code{gradtol}', '\code{stepmax}', '\code{steptol}', or '\code{iterlim}' -- see \code{nlm}'s documentation
#' for more detail)
#' @param fit \code{TRUE} if an HMM should be fitted to the data, \code{FALSE} otherwise.
#' If fit=\code{FALSE}, a model is returned with the MLE replaced by the initial parameters given in
#' input. This option can be used to assess the initial parameters. Default: \code{TRUE}.
#' @param userBounds User specified bounds for all parameters in order '\code{rbind(stepBounds,angleBounds,omegaBounds,dryBounds,diveBounds,iceBounds,landBounds)}'
#' @param stepDM stepPar design matrix
#' @param angleDM anglePar design matrix
#' @param omegaDM omegaPar design matrix
#' @param dryDM dryPar design matrix
#' @param diveDM divePar design matrix
#' @param iceDM icePar design matrix
#' @param landDM landPar design matrix
#' @param cons Optional list specifying power to raise parameters corresponding to each column of the design matrix for each data stream. While there could be other uses, primarily intended to constrain specific parameters to be positive by setting cons=2.  Default=NULL, which simply raises all parameters to the power of 1.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param logitcons Lower bound for wrpcauchy turning angle concentration parameter (on logit scale)
#'
#' @return A \code{momentuHMM} object, i.e. a list of:
#' \item{mle}{The maximum likelihood estimates of the parameters of the model (if the numerical algorithm
#' has indeed identified the global maximum of the likelihood function), which is a list
#' of: \code{stepPar} (step distribution parameters), \code{anglePar} (angle distribution
#' parameters), \code{beta} (transition probabilities regression coefficients - more information
#' in "Details"), and \code{delta} (initial distribution).}
#' \item{CI_real}{Standard errors and 95\% confidence intervals on the real (i.e., natural) scale of parameters}
#' \item{CI_beta}{Standard errors and 95\% confidence intervals on the beta (i.e., working) scale of parameters}
#' \item{data}{The movement data}
#' \item{mod}{The object returned by the numerical optimizer \code{nlm}}
#' \item{conditions}{Conditions used to fit the model (e.g., parameter bounds, distributions, \code{zeroInflation}, \code{estAngleMean},
#' \code{stationary}, and \code{formula})}
#' \item{rawCovs}{Raw covariate values, as found in the data (if any). Used in \code{\link{plot.momentuHMM}}.}
#'
#' @details
#' \itemize{
#' \item The matrix \code{beta} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix -
#' filled in row-wise).
#' In a covariate-free model (default), \code{beta} has one row, for the intercept.
#'
#' \item The choice of initial parameters is crucial to fit a model. The algorithm might not find
#' the global optimum of the likelihood function if the initial parameters are poorly chosen.
#' }
#'
#' @examples
#' ### 1. simulate data
#' # define all the arguments of simData
#' nbAnimals <- 2
#' nbStates <- 2
#' nbCovs <- 2
#' mu<-c(15,50)
#' sigma<-c(10,20)
#' angleMean <- c(pi,0)
#' kappa <- c(0.7,1.5)
#' stepPar <- c(mu,sigma)
#' anglePar <- c(angleMean,kappa)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' zeroInflation <- FALSE
#' obsPerAnimal <- c(50,100)
#'
#' data <- simData(nbAnimals=nbAnimals,nbStates=nbStates,stepDist=stepDist,angleDist=angleDist,
#'                  stepPar=stepPar,anglePar=anglePar,nbCovs=nbCovs,zeroInflation=zeroInflation,
#'                  obsPerAnimal=obsPerAnimal)
#'
#' ### 2. fit the model to the simulated data
#' # define initial values for the parameters
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
#' anglePar <- kappa0 # the angle mean is not estimated, so only the concentration parameter is needed
#' formula <- ~cov1+cos(cov2)
#'
#' m <- fitHMM(data=data,nbStates=nbStates,stepPar=stepPar,anglePar=anglePar,formula=formula,
#'               stepDist=stepDist,angleDist=angleDist,angleMean=angleMean)
#'
#' print(m)
#'
#' @references
#' Patterson T.A., Basson M., Bravington M.V., Gunn J.S. 2009.
#' Classifying movement behaviour in relation to environmental conditions using hidden Markov models.
#' Journal of Animal Ecology, 78 (6), 1113-1123.
#'
#' Langrock R., King R., Matthiopoulos J., Thomas L., Fortin D., Morales J.M. 2012.
#' Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions.
#' Ecology, 93 (11), 2336-2342.
#' 
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats model.matrix nlm terms
#' @import CircStats
#'
#' @useDynLib momentuHMM

fitHMM <- function(data,nbStates,dist,
                   Par,beta0=NULL,delta0=NULL,
                   estAngleMean=NULL,
                   formula=~1,
                   stationary=FALSE,verbose=0,nlmPar=NULL,fit=TRUE,
                   DM=NULL,cons=NULL,userBounds=NULL,logitcons=NULL,stateNames=NULL)
{
  # check that the data is a moveData object
  if(!is.moveData(data))
    stop("'data' must be a moveData object (as output by prepData or simData)")

  # check that the formula is a formula
  is.formula <- function(x)
    tryCatch(inherits(x,"formula"),error= function(e) {FALSE})

  if(!is.formula(formula))
    stop("Check the argument 'formula'.")

  # check that there is no response varibale in the formula
  if(attr(terms(formula),"response")!=0)
    stop("The response variable should not be specified in the formula.")
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
  distnames<-names(dist)
  if(any(is.na(match(distnames,names(data))))) stop(paste0(distnames[is.na(match(distnames,names(data)))],collapse=", ")," not found in data")
  if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
  Par <- Par[distnames]
  
  # build design matrix for t.p.m.
  covsCol <- seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
  covs <- model.matrix(formula,data)
  if(nrow(covs)!=nrow(data)) stop("covariates cannot contain missing values")
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  # check that stationary==FALSE if there are covariates
  if(nbCovs>0 & stationary==TRUE)
    stop("stationary can't be set to TRUE if there are covariates.")

  if(length(covsCol)>0) {
    rawCovs <- data[covsCol]
  }
  else {
    rawCovs <- NULL
  }
  
  mHind <- (is.null(DM) & is.null(userBounds)) # indicator for moveHMMwrap below
  
  
  #####################
  ## Check arguments ##
  #####################
  for(i in distnames){
    dist[[i]]<-match.arg(dist[[i]],momentuHMMdists)
  }

  # determine whether zero-inflation should be included
  zeroInflation <- vector('list',length(distnames))
  names(zeroInflation) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% zeroInflationdists){
      if(length(which(data[[i]]==0))>0) {
        zeroInflation[[i]]<-TRUE
        # check that zero-mass is in the open interval (0,1)
        zm0 <- Par[[i]][(length(Par[[i]])-nbStates+1):length(Par[[i]])]
        zm0[which(zm0==0)] <- 1e-8
        zm0[which(zm0==1)] <- 1-1e-8
        Par[[i]][(length(Par[[i]])-nbStates+1):length(Par[[i]])] <- zm0
      }
      else 
        zeroInflation[[i]]<-FALSE
    }
    else zeroInflation[[i]]<-FALSE
  }
  
  if(nbStates<0)
    stop("nbStates should be at least 1.")
  if(length(data)<1 | any(dim(data)<1))
    stop("The data input is empty.")
  
  if(is.null(estAngleMean)){
    estAngleMean <- vector('list',length(distnames))
    names(estAngleMean) <- distnames
  } else {
    if(!is.list(estAngleMean) | is.null(names(estAngleMean))) stop("'estAngleMean' must be a named list")
  }
  for(i in distnames){
    if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
  }
  for(i in distnames[which(!(dist %in% angledists))]){
    estAngleMean[[i]] <- FALSE
  }
  estAngleMean<-estAngleMean[distnames]
  
  
  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)

  par0 <- unlist(Par)
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,userBounds)
  parSize <- p$parSize

  if(!is.null(beta0)) {
    if(ncol(beta0)!=nbStates*(nbStates-1) | nrow(beta0)!=nbCovs+1) {
      error <- paste("beta0 has wrong dimensions: it should have",nbCovs+1,"rows and",
                     nbStates*(nbStates-1),"columns.")
      stop(error)
    }
  }

  if(!is.null(delta0))
    if(length(delta0)!=nbStates)
      stop(paste("delta0 has the wrong length: it should have",nbStates,"elements."))
  
  # check that verbose is in {0,1,2}
  if(!(verbose %in% c(0,1,2)))
    stop("verbose must be in {0,1,2}")

  # check that observations are within expected bounds
  for(i in which(unlist(lapply(dist,function(x) x %in% nonnegativedists))))
    if(length(which(data[[distnames[[i]]]]<0))>0)
      stop(distnames[[i]]," data should be non-negative")
  
  for(i in which(unlist(lapply(dist,function(x) x %in% angledists))))
    if(length(which(data[[distnames[[i]]]] < -pi | data[[distnames[[i]]]] > pi))>0)
      stop(distnames[[i]]," data should be between -pi and pi")
  
  for(i in which(unlist(lapply(dist,function(x) x %in% "beta"))))
    if(length(which(data[[distnames[[i]]]]<0 | data[[distnames[[i]]]]>=1))>0)
      stop(distnames[[i]]," data should be between 0 and 1")

  for(i in which(unlist(lapply(dist,function(x) x %in% "pois"))))
    if(!isTRUE(all.equal(data[[distnames[[i]]]],as.integer(data[[distnames[[i]]]]))))
      stop(distnames[[i]]," data should be non-negative integers")

  # check elements of nlmPar
  lsPars <- c("gradtol","stepmax","steptol","iterlim")
  if(length(which(!(names(nlmPar) %in% lsPars)))>0)
    stop("Check the names of the element of 'nlmPar'; they should be in
         ('gradtol','stepmax','steptol','iterlim')")

  ####################################
  ## Prepare initial values for nlm ##
  ####################################
  if(is.null(beta0) & nbStates>1) {
    beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                    nrow=nbCovs+1,byrow=TRUE)
  }

  if(is.null(delta0))
    delta0 <- rep(1,nbStates)/nbStates
  if(stationary)
    delta0 <- NULL

  bounds <- p$bounds

  if(is.null(DM)){
    DM <- cons <- logitcons <- vector('list',length(distnames))
    names(DM) <- names(cons) <- names(logitcons) <- distnames
  } else {
    if(!is.list(DM) | is.null(names(DM))) stop("'DM' must be a named list")
    if(!any(names(DM) %in% distnames)) stop("DM names must include at least one of: ",paste0(distnames,collapse=", "))
  }
  
  for(i in distnames){
    if(is.null(DM[[i]]) & length(Par[[i]])!=(parSize[[i]]*nbStates)){
      error<-paste0("Wrong number of initial parameters -- there should be ",parSize[[i]]*nbStates," initial ",i," parameters")
      if(zeroInflation[[i]]) error<-paste0(error," -- zero-mass parameters should be included")
      stop(error)
    }
  }
  
  fullDM<-getDM(data,DM,dist,nbStates,p$parNames,bounds,Par)
  DMind <- lapply(fullDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))
  
  for(i in distnames){
    if(nrow(fullDM[[i]])!=(parSize[[i]]*nbStates)){
      error<- paste0("DM for ",i," must have ",(parSize[[i]]*nbStates)," rows")
      if(zeroInflation[[i]])
        stop(paste0(error,". Should zero inflation parameters be included?"))
      else stop(error)
    }
  }
  
  if(any(unlist(lapply(Par,length))!=unlist(lapply(fullDM,ncol))))
    stop("Dimension mismatch between Par and DM for: ",paste(names(which(unlist(lapply(Par,length))!=unlist(lapply(fullDM,ncol)))),collapse=", "))
  
  if(sum((unlist(parSize)>0)*unlist(lapply(fullDM,ncol)))!=length(par0)) {
    error <- "Wrong number of initial parameters"
    stop(error)
  }
  
  if(is.null(cons)){
    cons <- vector('list',length(distnames))
    names(cons) <- distnames
  } else {
    if(!is.list(cons) | is.null(names(cons))) stop("'cons' must be a named list")
  }
  for(i in distnames){
    if(is.null(cons[[i]])) cons[[i]] <- rep(1,ncol(fullDM[[i]]))
  }
  cons<-cons[distnames]
  if(any(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and cons for: ",paste(names(which(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  
  if(is.null(logitcons)){
    logitcons <- vector('list',length(distnames))
    names(logitcons) <- distnames
  } else {
    if(!is.list(logitcons) | is.null(names(logitcons))) stop("'logitcons' must be a named list")
  }
  for(i in distnames){
    if(is.null(logitcons[[i]])) logitcons[[i]] <- rep(0,ncol(fullDM[[i]]))
  }
  for(i in which(!(dist %in% "wrpcauchy"))){
    logitcons[[distnames[i]]]<-rep(0,ncol(fullDM[[distnames[i]]]))
  }
  logitcons<-logitcons[distnames]
  if(any(unlist(lapply(logitcons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and logitcons for: ",paste(names(which(unlist(lapply(logitcons,length))!=unlist(lapply(Par,length)))),collapse=", "))

  # build the vector of initial working parameters
  wpar <- n2w(Par,bounds,beta0,delta0,nbStates,estAngleMean,DM,cons,logitcons)
  
  if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")
  
  ##################
  ## Optimization ##
  ##################
  # this function is used to muffle the warning "NA/Inf replaced by maximum positive value" in nlm and "value out of range in 'lgamma'" in nLogLike_rcpp
  h <- function(w) {
    if(any(grepl("NA/Inf replaced by maximum positive value",w)) | any(grepl("value out of range in 'lgamma'",w)))
      invokeRestart("muffleWarning")
  }

  # just use moveHMM if simpler models are specified
  if(all(distnames %in% c("step","angle")) & mHind){
    out<-moveHMMwrap(data,nbStates,dist,Par,beta0,delta0,estAngleMean,formula,stationary,verbose,nlmPar,fit)
    mod<-out$mod
    mle<-out$mle
  } 
  else if(fit) {
    # check additional parameters for nlm
    gradtol <- ifelse(is.null(nlmPar$gradtol),1e-6,nlmPar$gradtol)
    typsize = rep(1, length(wpar))
    defStepmax <- max(1000 * sqrt(sum((wpar/typsize)^2)),1000)
    stepmax <- ifelse(is.null(nlmPar$stepmax),defStepmax,nlmPar$stepmax)
    steptol <- ifelse(is.null(nlmPar$steptol),1e-6,nlmPar$steptol)
    iterlim <- ifelse(is.null(nlmPar$iterlim),1000,nlmPar$iterlim)

    # call to optimizer nlm
    withCallingHandlers(mod <- nlm(nLogLike,wpar,nbStates,formula,bounds,parSize,data,dist,covs,
                                   estAngleMean,zeroInflation,
                                   stationary,cons,fullDM,DMind,logitcons,
                                   print.level=verbose,gradtol=gradtol,
                                   stepmax=stepmax,steptol=steptol,
                                   iterlim=iterlim,hessian=TRUE),
                        warning=h) # filter warnings using function h

    # convert the parameters back to their natural scale
    mle <- w2n(mod$estimate,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,fullDM,DMind,logitcons,nrow(data),dist)
  }
  else {
    mod <- NA
    mle <- w2n(wpar,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,fullDM,DMind,logitcons,nrow(data),dist)
  }

  ####################
  ## Prepare output ##
  ####################
  if(is.null(stateNames)){
    for(i in 1:nbStates)
      stateNames[i] <- paste("state",i)
  }
  
  # name columns and rows of MLEs
  parindex <- c(0,cumsum(unlist(lapply(fullDM,ncol)))[-length(fullDM)])
  names(parindex) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% angledists)
      if(!estAngleMean[[i]]){
        p$parNames[[i]] <- c("mean",p$parNames[[i]])
      }
    if(DMind[[i]]){
      mle[[i]]<-matrix(mle[[i]][,1],nrow=length(p$parNames[[i]]),ncol=nbStates,byrow=TRUE)
      rownames(mle[[i]]) <- p$parNames[[i]]
      colnames(mle[[i]]) <- stateNames
    } else {
      mle[[i]]<-mod$estimate[parindex[[i]]+1:ncol(fullDM[[i]])]
      names(mle[[i]])<-colnames(fullDM[[i]])
      if(is.null(names(mle[[i]]))) warning("No names for the regression coeffs were provided in DM$",i)
    }
  }

  if(!is.null(mle$beta)) {
    rownames(mle$beta) <- c("(Intercept)",attr(terms(formula),"term.labels"))
    columns <- NULL
    for(i in 1:nbStates)
      for(j in 1:nbStates) {
        if(i<j)
          columns[(i-1)*nbStates+j-i] <- paste(i,"->",j)
        if(j<i)
            columns[(i-1)*(nbStates-1)+j] <- paste(i,"->",j)
      }
    colnames(mle$beta) <- columns
  }

  # compute stationary distribution
  if(stationary) {
    gamma <- trMatrix_rcpp(nbStates,mle$beta,covs)[,,1]

    # error if singular system
    tryCatch(
      mle$delta <- solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)),
      error = function(e) {
        stop(paste("A problem occurred in the calculation of the stationary",
                   "distribution. You may want to try different initial values",
                   "and/or the option stationary=FALSE."))
      }
    )
  }

  if(nbStates==1)
    mle$delta <- 1

  # compute t.p.m. if no covariates
  if(nbCovs==0 & nbStates>1) {
    trMat <- trMatrix_rcpp(nbStates,mle$beta,covs)
    mle$gamma <- trMat[,,1]
    colnames(mle$gamma)<-stateNames
    rownames(mle$gamma)<-stateNames
  }

  # conditions of the fit
  conditions <- list(dist=dist,zeroInflation=zeroInflation,
                     estAngleMean=estAngleMean,stationary=stationary,formula=formula,cons=cons,bounds=p$bounds,DM=fullDM,logitcons=logitcons)

  mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames)
  
  #CI_real<-CI_real(momentuHMM(mh))
  #CI_beta<-CI_beta(momentuHMM(mh))
  
  #mh <- list(data=data,mle=mle,CI_real=CI_real,CI_beta=CI_beta,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames)
  mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames)
  
  return(momentuHMM(mh))
}
