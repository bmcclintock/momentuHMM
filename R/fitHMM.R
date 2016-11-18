
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
#' @return A \code{moveHMM} object, i.e. a list of:
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
#' \item{rawCovs}{Raw covariate values, as found in the data (if any). Used in \code{\link{plot.moveHMM}}.}
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

fitHMM <- function(data,nbStates,stepPar,anglePar,omegaPar=NULL,dryPar=NULL,divePar=NULL,icePar=NULL,landPar=NULL,beta0=NULL,delta0=NULL,formula=~1,
                   stepDist=c("gamma","weibull"),angleDist=c("wrpcauchy","none"),
                   omegaDist=c("none","beta"),dryDist=c("none","beta"),diveDist=c("none","pois"),iceDist=c("none","beta"),landDist=c("none","beta"),
                   angleMean=NULL,stationary=FALSE,verbose=0,nlmPar=NULL,fit=TRUE,userBounds=NULL,
                   stepDM=NULL,angleDM=NULL,omegaDM=NULL,dryDM=NULL,diveDM=NULL,iceDM=NULL,landDM=NULL,cons=NULL,stateNames=NULL,logitcons=0)
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

  # build design matrix
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="omega" & names(data)!="dry" & names(data)!="dive" & names(data)!="ice" & names(data)!="land")
  covs <- model.matrix(formula,data)

  if(length(covsCol)>0) {
    rawCovs <- data[covsCol]
    data <- cbind(data[-covsCol],covs)
  }
  else {
    rawCovs <- NULL
    data <- cbind(data,covs)
  }

  nbCovs <- ncol(covs)-1 # substract intercept column

  # determine whether zero-inflation should be included
  if(length(which(data$step==0))>0)
    zeroInflation <- TRUE
  else
    zeroInflation <- FALSE

  # check that zero-mass is in the open interval (0,1)
  if(zeroInflation) {
    zm0 <- stepPar[(length(stepPar)-nbStates+1):length(stepPar)]
    zm0[which(zm0==0)] <- 1e-8
    zm0[which(zm0==1)] <- 1-1e-8
    stepPar[(length(stepPar)-nbStates+1):length(stepPar)] <- zm0
  }

  #####################
  ## Check arguments ##
  #####################
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  omegaDist <- match.arg(omegaDist)
  dryDist <- match.arg(dryDist)
  diveDist <- match.arg(diveDist)
  iceDist <- match.arg(iceDist)
  landDist <- match.arg(landDist)
  
  if(nbStates<0)
    stop("nbStates should be at least 1.")
  if(length(data)<1)
    stop("The data input is empty.")
  if(is.null(data$step))
    stop("Missing field in data: step.")
  
  if(is.null(stepDM)) stepDM <- diag((2+zeroInflation)*nbStates)
  if(is.null(angleDM)) angleDM <- diag((1+is.null(angleMean))*nbStates)
  if(is.null(omegaDM)) omegaDM <- diag(2*nbStates)
  if(is.null(dryDM)) dryDM <- diag(2*nbStates)
  if(is.null(diveDM)) diveDM <- diag(nbStates)
  if(is.null(iceDM)) iceDM <- diag(2*nbStates)
  if(is.null(landDM)) landDM <- diag(2*nbStates)
  if(is.null(cons)) cons <- list(step=rep(1,ncol(stepDM)),angle=rep(1,ncol(angleDM)),omega=rep(1,ncol(omegaDM)),dry=rep(1,ncol(dryDM)),dive=rep(1,ncol(diveDM)),ice=rep(1,ncol(iceDM)),land=rep(1,ncol(landDM)))
  
  #if(length(stepPar[1:(length(stepPar)-nbStates*zeroInflation)])!=ncol(stepDM))
  #  stop("Dimension mismatch between stepPar and stepDM")
  #if(length(anglePar[1:(length(anglePar)-nbStates*(is.null(angleMean)))])!=ncol(angleDM))
  #  stop("Dimension mismatch between anglePar and angleDM")
  if(omegaDist!="none" & length(omegaPar)!=ncol(omegaDM))
    stop("Dimension mismatch between omegaPar and omegaDM")
  else if (omegaDist!="none" & length(cons$omega)!=ncol(omegaDM))
    stop("Dimension mismatch between omega constraints (cons$omega) and omegaDM")
  if(dryDist!="none" & length(dryPar)!=ncol(dryDM))
    stop("Dimension mismatch between dryPar and dryDM")
  else if (dryDist!="none" & length(cons$dry)!=ncol(dryDM))
    stop("Dimension mismatch between dry constraints (cons$dry) and dryDM")
  if(diveDist!="none" & length(divePar)!=ncol(diveDM))
    stop("Dimension mismatch between divePar and diveDM")
  else if (diveDist!="none" & length(cons$dive)!=ncol(diveDM))
    stop("Dimension mismatch between dive constraints (cons$dive) and diveDM")
  if(iceDist!="none" & length(icePar)!=ncol(iceDM))
    stop("Dimension mismatch between icePar and iceDM")
  else if (iceDist!="none" & length(cons$ice)!=ncol(iceDM))
    stop("Dimension mismatch between ice constraints (cons$ice) and iceDM")
  if(landDist!="none" & length(landPar)!=ncol(landDM))
    stop("Dimension mismatch between landPar and landDM")
  else if (landDist!="none" & length(cons$land)!=ncol(landDM))
    stop("Dimension mismatch between land constraints (cons$land) and landDM")
  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)

  par0 <- c(stepPar,anglePar,omegaPar,dryPar,divePar,icePar,landPar)  
  #if(!is.null(userBounds) & !is.numeric(userBounds)){
  #  evalBounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2)
  #}
  #p <- parDef(stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,nbStates,is.null(angleMean),zeroInflation,evalBounds,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM)
  p <- parDef(stepDist,angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,nbStates,is.null(angleMean),zeroInflation,userBounds,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM)
  bounds <- p$bounds
  if(!is.numeric(bounds)){
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(p$bounds)))
  }
  parSize <- p$parSize
  if(length(stepPar%*%stepDM)!=(parSize[1]*nbStates))
    stop("zero inflation parameters must be included in 'stepPar'")
  if(sum((parSize>0)*c(ncol(stepDM),ncol(angleDM),ncol(omegaDM),ncol(dryDM),ncol(diveDM),ncol(iceDM),ncol(landDM)))!=length(par0)) {
    error <- "Wrong number of initial parameters"
    if(ncol(stepDM)!=length(stepPar)) {
      error <- paste(error,"-- there should be",(parSize[1]-2)*nbStates+ncol(stepDM),"initial step parameters")
      if(zeroInflation)
        error <- paste(error,"-- zero-mass parameters should be included")
    }
    if(angleDist!="none" & ncol(angleDM)!=length(anglePar))
      error <- paste(error,"-- there should be",(parSize[2]-1)*nbStates+ncol(angleDM),"initial angle parameters.")
    if(angleDist=="none" & length(anglePar)>0)
      error <- paste(error,"-- 'anglePar' should be NULL.")
    stop(error)
  }
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

  stepBounds <- bounds[1:ncol(stepDM),]
  if(length(which(stepPar<=stepBounds[,1] | stepPar>=stepBounds[,2]))>0)
    stop(paste("Check the step parameters bounds (the initial parameters should be",
               "strictly between the bounds of their parameter space)."))

  if(angleDist!="none") {
    # We can't really write distribution-agnostic code here, because the bounds
    # defined in parDef are not the actual bounds of the parameter space.
    if(is.null(angleMean)) {
      m <- anglePar[1:nbStates] # angle mean
      k <- anglePar[(nbStates+1):length(anglePar)] # angle concentration
      if(length(which(m<=(-pi) | m>pi))>0)
        stop("Check the angle parameters bounds. The angle mean should be in (-pi,pi].")
    } else {
      k <- anglePar # angle concentration
      if(length(which(angleMean<=-pi | angleMean>pi))>0)
        stop("The 'angleMean' should be in (-pi,pi].")
      if(length(angleMean)!=nbStates)
        stop("The argument 'angleMean' should be of length nbStates.")
    }
    if(angleDist!="wrpcauchy" & length(which(k<=0))>0)
      stop("Check the angle parameters bounds. The concentration should be strictly positive.")
    if(angleDist=="wrpcauchy" & ((length(which(k>=1))>0) | length(which(k<= -1))>0))
      stop("Check the angle parameters bounds. The concentration should be in (-1,1).")
  }
  else if(!is.null(angleMean))
    stop("'angleMean' shouldn't be specified if 'angleDist' is \"none\"")

  if(length(which(par0<bounds[,1] | par0>bounds[,2]))>0)
    stop(paste("Check the parameter bounds (the initial parameters should be",
               "strictly between the bounds of their parameter space)."))

  # check that verbose is in {0,1,2}
  if(!(verbose %in% c(0,1,2)))
    stop("verbose must be in {0,1,2}")

  # check that observations are within expected bounds
  if(length(which(data$step<0))>0)
    stop("The step lengths should be positive.")
  if(length(which(data$angle < -pi | data$angle > pi))>0)
    stop("The turning angles should be between -pi and pi.")

  # check that stationary==FALSE if there are covariates
  if(nbCovs>0 & stationary==TRUE)
    stop("stationary can't be set to TRUE if there are covariates.")

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

  estAngleMean <- (is.null(angleMean) & angleDist!="none")

  if(!is.numeric(p$bounds)){
    bounds <- gsub("step","",p$bounds,fixed=TRUE)
    bounds <- gsub("angle","",bounds,fixed=TRUE)
    bounds <- gsub("omega","",bounds,fixed=TRUE)
    bounds <- gsub("dry","",bounds,fixed=TRUE) 
    bounds <- gsub("dive","",bounds,fixed=TRUE)    
    bounds <- gsub("ice","",bounds,fixed=TRUE)
    bounds <- gsub("land","",bounds,fixed=TRUE)
  }
  
  # build the vector of initial working parameters
  wpar <- n2w(par0,stepDist,bounds,beta0,delta0,nbStates,estAngleMean,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,cons,logitcons)
      
  ##################
  ## Optimization ##
  ##################
  # this function is used to muffle the warning "NA/Inf replaced by maximum positive value" in nlm
  h <- function(w) {
    if(any(grepl("NA/Inf replaced by maximum positive value",w)))
      invokeRestart("muffleWarning")
  }

  if(fit) {
    # check additional parameters for nlm
    gradtol <- ifelse(is.null(nlmPar$gradtol),1e-6,nlmPar$gradtol)
    typsize = rep(1, length(wpar))
    defStepmax <- max(1000 * sqrt(sum((wpar/typsize)^2)),1000)
    stepmax <- ifelse(is.null(nlmPar$stepmax),defStepmax,nlmPar$stepmax)
    steptol <- ifelse(is.null(nlmPar$steptol),1e-6,nlmPar$steptol)
    iterlim <- ifelse(is.null(nlmPar$iterlim),1000,nlmPar$iterlim)

    # call to optimizer nlm
    withCallingHandlers(mod <- nlm(nLogLike,wpar,nbStates,formula,bounds,parSize,data,stepDist,
                                   angleDist,omegaDist,dryDist,diveDist,iceDist,landDist,
                                   angleMean,zeroInflation,stationary,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,p$boundInd,logitcons,
                                   print.level=verbose,gradtol=gradtol,
                                   stepmax=stepmax,steptol=steptol,
                                   iterlim=iterlim,hessian=TRUE),
                        warning=h) # filter warnings using function h

    # convert the parameters back to their natural scale
    mle <- w2n(mod$estimate,stepDist,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,p$boundInd,logitcons)
  }
  else {
    mod <- NA
    mle <- w2n(wpar,stepDist,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary,cons,stepDM,angleDM,omegaDM,dryDM,diveDM,iceDM,landDM,p$boundInd,logitcons)
  }

  ####################
  ## Prepare output ##
  ####################
  # include angle mean if it wasn't estimated
  if(!is.null(angleMean) & angleDist!="none")
    mle$anglePar <- rbind(angleMean,mle$anglePar)

  # name columns and rows of MLEs
  rownames(mle$stepPar) <- p$parNames[1:nrow(mle$stepPar)]

  if(is.null(stateNames)){
    for(i in 1:nbStates)
      stateNames[i] <- paste("state",i)
  }
  colnames(mle$stepPar) <- stateNames

  if(angleDist!="none") {
    rownames(mle$anglePar) <- c("mean","concentration")
    colnames(mle$anglePar) <- stateNames
  }
  
  if(omegaDist!="none") {
    rownames(mle$omegaPar) <- c("shape1","shape2")
    colnames(mle$omegaPar) <- stateNames
  }
  
  if(dryDist!="none") {
    rownames(mle$dryPar) <- c("shape1","shape2")
    colnames(mle$dryPar) <- stateNames
  }    

  if(diveDist!="none") {  
    rownames(mle$divePar) <- "lambda"
    colnames(mle$divePar) <- stateNames
  }
  
  if(iceDist!="none") {
    rownames(mle$icePar) <- c("shape1","shape2")
    colnames(mle$icePar) <- stateNames
  }    
  
  if(landDist!="none") {
    rownames(mle$landPar) <- c("shape1","shape2")
    colnames(mle$landPar) <- stateNames
  }    

  if(!is.null(mle$beta)) {
    rownames(mle$beta) <- c("intercept",attr(terms(formula),"term.labels"))
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
  conditions <- list(stepDist=stepDist,angleDist=angleDist,omegaDist=omegaDist,dryDist=dryDist,diveDist=diveDist,iceDist=iceDist,landDist=landDist,zeroInflation=zeroInflation,
                     estAngleMean=estAngleMean,stationary=stationary,formula=formula,cons=cons,stepDM=stepDM,angleDM=angleDM,omegaDM=omegaDM,dryDM=dryDM,diveDM=diveDM,iceDM=iceDM,landDM=landDM,cons=cons,logitcons=logitcons)

  mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,bounds=p$bounds,rawCovs=rawCovs,stateNames=stateNames)
  
  CI_real<-CI_real(momentuHMM(mh))
  CI_beta<-CI_beta(momentuHMM(mh))
  
  mh <- list(data=data,mle=mle,CI_real=CI_real,CI_beta=CI_beta,mod=mod,conditions=conditions,bounds=p$bounds,rawCovs=rawCovs,stateNames=stateNames)
  
  return(momentuHMM(mh))
}
