
#' Simulation tool
#'
#' Simulates movement data from an HMM.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are: gamma, weibull, lnorm, exp. Default: gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are: vm, wrpcauchy. Set to \code{"none"} if the angle distribution should
#' not be estimated. Default: vm.
#' @param stepPar Parameters of the step length distribution.
#' @param anglePar Parameters of the turning angle distribution.
#' @param beta Matrix of regression parameters for the transition probabilities (more information
#' in "Details").
#' @param covs Covariate values to include in the model, as a dataframe. Default: \code{NULL}.
#' Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL}, and specifying \code{nbCovs>0}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified of
#' \code{covs} is specified.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default: \code{FALSE}. If \code{TRUE}, values for the zero-mass parameters should be
#' included in \code{stepPar}.
#' @param obsPerAnimal Either the number of the number of observations per animal (if single value),
#' or the bounds of the number of observations per animal (if vector of two values). In the latter case,
#' the numbers of obervations generated for each animal are uniformously picked from this interval.
#' Default: \code{c(500,1500)}.
#' @param model A momentuHMM object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for nbAnimals,
#' obsPerAnimal, covs (if covariate values different from those in the data should be specified),
#' and states.
#' @param states \code{TRUE} if the simulated states should be returned, \code{FALSE} otherwise (default).
#'
#' @return An object moveData, i.e. a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{step}{The step lengths}
#' \item{angle}{The turning angles (if any)}
#' \item{x}{Either easting or longitude}
#' \item{y}{Either norting or latitude}
#' \item{...}{Covariates (if any)}
#'
#' @details \itemize{
#' \item The matrix \code{beta} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix - filled in
#' row-wise).
#' In a covariate-free model (default), \code{beta} has one row, for the intercept.
#'
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal', then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' }
#'
#' @examples
#' # 1. Pass a fitted model to simulate from
#' # (m is a momentuHMM object - as returned by fitHMM - automatically loaded with the package)
#' # We keep the default nbAnimals=1.
#' m <- example$m
#' obsPerAnimal=c(50,100)
#' data <- simData(model=m,obsPerAnimal=obsPerAnimal)
#'
#' # 2. Pass the parameters of the model to simulate from
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' anglePar <- c(pi,0,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                anglePar=anglePar,nbCovs=2,zeroInflation=TRUE,obsPerAnimal=obsPerAnimal)
#'
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(pi,0,0.5,0.7) # mean1, mean2, k1, k2
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                anglePar=anglePar,obsPerAnimal=obsPerAnimal)
#'
#' # step length only and zero-inflation
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' stepDist <- "gamma"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist="none",stepPar=stepPar,
#'                nbCovs=2,zeroInflation=TRUE,obsPerAnimal=obsPerAnimal)
#'
#' # include covariates
#' # (note that it is useless to specify "nbCovs", which respectively determined
#' # by the number of columns of "cov")
#' cov <- data.frame(temp=rnorm(500,20,5))
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(pi,0,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                 anglePar=anglePar,covs=cov)
#'
#' @export
#' @importFrom stats rnorm runif

simData <- function(nbAnimals=1,nbStates=2,dist,
                    Par,beta=NULL,
                    covs=NULL,nbCovs=0,zeroInflation=NULL,obsPerAnimal=c(500,1500),
                    DM=NULL,cons=NULL,userBounds=NULL,logitcons=NULL,stateNames=NULL,
                    model=NULL,states=FALSE)
{
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    # extract simulation parameters from model
    nbStates <- length(model$stateNames)
    dist<-model$conditions$dist
    distnames<-names(dist)
    userBounds <- model$bounds
    stateNames<-model$stateNames
  
    Par <- model$mle[distnames]
    for(i in distnames[which(dist %in% c("wrpcauchy","vm"))]){
      Par[[i]]<-Par[[i]][2,,drop=FALSE]
    }
    beta <- model$mle$beta
    
    DM <- model$conditions$DM
    cons <- model$conditions$cons
    logitcons <- model$conditions$logitcons
    zeroInflation <- model$conditions$zeroInflation

    if(is.null(covs)) {
      covsCol <- seq(1,ncol(model$data))[-match(c("ID","x","y",distnames),names(model$data))]
      covs <- model.matrix(model$conditions$formula,model$data)

      if(length(covsCol)>1) {
        # remove intercept column, which is not expected in 'covs'
        names <- colnames(covs)
        covs <- data.frame(covs[,-1]) # data.frame structure is lost when only one column
        colnames(covs) <- names[-1]
      } else
        covs <- NULL
    }
    # else, allow user to enter new values for covariates

  } else {
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
    distnames<-names(dist)
    if(length(setdiff(distnames,names(Par))) | (length(dist)!=length(Par))) stop("Length and names of 'dist' and 'Par' must match")
    Par<-Par[distnames]
  }
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
  
  #####################
  ## Check arguments ##
  #####################
  eval(parse(text=paste0("dist$",distnames,"=match.arg(dist$",distnames,",c('gamma','weibull','exp','beta','pois','wrpcauchy','vm'))")))
  Fun <- lapply(dist,function(x) paste("r",x,sep=""))

  
  if(nbAnimals<1)
    stop("nbAnimals should be at least 1.")
  if(nbStates<1)
    stop("nbStates should be at least 1.")

  if(is.null(DM)){
    DM <- vector('list',length(distnames))
    names(DM) <- distnames
  } else if(!is.list(DM) | is.null(names(DM))) stop("'DM' must be a named list")
  eval(parse(text=paste0("if(is.null(",DM[distnames],")) DM$",distnames," <- diag(",ifelse(dist=="exp" | dist=="pois" | dist=="wrpcauchy",1,2),"*nbStates)")))
  DM<-DM[distnames]
  if(any(unlist(lapply(Par,length))!=unlist(lapply(DM,ncol))))
    stop("Dimension mismatch between Par and DM for: ",paste(names(which(unlist(lapply(Par,length))!=unlist(lapply(DM,ncol)))),collapse=", "))
  
  if(is.null(cons)){
    cons <- vector('list',length(distnames))
    names(cons) <- distnames
  } else if(!is.list(cons) | is.null(names(cons))) stop("'cons' must be a named list")
  eval(parse(text=paste0("if(is.null(",cons[distnames],")) cons$",distnames," <- rep(1,ncol(DM$",distnames,"))")))
  cons<-cons[distnames]
  if(any(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and cons for: ",paste(names(which(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  
  if(is.null(logitcons)){
    logitcons <- vector('list',length(distnames))
    names(logitcons) <- distnames
  } else if(!is.list(logitcons) | is.null(names(logitcons))) stop("'logitcons' must be a named list")
  eval(parse(text=paste0("if(is.null(",logitcons[distnames],")) logitcons$",distnames," <- rep(0,ncol(DM$",distnames,"))")))
  for(i in which(is.na(match(dist,"wrpcauchy")))){
    logitcons[[distnames[i]]]<-rep(0,ncol(DM[[distnames[i]]]))
  }
  logitcons<-logitcons[distnames]
  if(any(unlist(lapply(logitcons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and logitcons for: ",paste(names(which(unlist(lapply(logitcons,length))!=unlist(lapply(Par,length)))),collapse=", "))

  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)
  
  if(is.null(zeroInflation)){
    zeroInflation <- vector('list',length(distnames))
    names(zeroInflation) <- distnames
    for(i in distnames){
      zeroInflation[[i]]<-FALSE
    }
  } else {
    if(!is.list(zeroInflation) | is.null(names(zeroInflation))) stop("'zeroInflation' must be a named list")
    for(i in distnames){
      if(dist[[i]]=="wrpcauchy" | dist[[i]]=="vm" | dist[[i]]=="pois"){
        zeroInflation[[i]]<-FALSE
      }
    }
  }
  eval(parse(text=paste0("if(is.null(",zeroInflation[distnames],")) zeroInflation$",distnames," <- FALSE")))
  if(!all(unlist(lapply(zeroInflation,is.logical)))) stop("zeroInflation must be a list of logical objects")

  p <- parDef(dist,nbStates,FALSE,zeroInflation,userBounds,DM)
  par0 <- unlist(Par) 
  
  bounds <- p$bounds
  for(i in distnames){
    if(!is.numeric(bounds[[i]])){
      bounds[[i]]<-matrix(sapply(bounds[[i]],function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(p$bounds[[i]])))
    }
  }
  parSize <- p$parSize
  if(sum((parSize>0)*unlist(lapply(DM,ncol)))!=length(par0)) {
    error <- "Wrong number of initial parameters"
    stop(error)
  }
  if(!is.null(beta)) {
    if(ncol(beta)!=nbStates*(nbStates-1) | nrow(beta)!=nbCovs+1) {
      error <- paste("beta has wrong dimensions: it should have",nbCovs+1,"rows and",
                     nbStates*(nbStates-1),"columns.")
      stop(error)
    }
  }
  
  for(i in distnames){
    if(length(which(Par[[i]]<bounds[[i]][,1] | Par[[i]]>bounds[[i]][,2]))>0)
      stop(paste0("Check the parameter bounds for '",i,"' (the initial parameters should be ",
                  "strictly between the bounds of their parameter space)."))
  }

  if(length(which(obsPerAnimal<1))>0)
    stop("obsPerAnimal should have positive values.")

  if(!is.null(covs) & nbCovs>0) {
    if(ncol(covs)!=nbCovs)
      warning("covs and nbCovs argument conflicting - nbCovs was set to ncol(covs)")
  }

  if(!is.null(covs)) {
    if(!is.data.frame(covs))
      stop("'covs' should be a data.frame")
  }

  if(!is.null(covs)) {
    nbCovs <- ncol(covs)

    # account for missing values of the covariates
    if(length(which(is.na(covs)))>0)
      warning(paste("There are",length(which(is.na(covs))),
                    "missing covariate values.",
                    "Each will be replaced by the closest available value."))
    for(i in 1:nbCovs) {
      if(length(which(is.na(covs[,i])))>0) { # if covariate i has missing values
        if(is.na(covs[1,i])) { # if the first value of the covariate is missing
          k <- 1
          while(is.na(covs[k,i])) k <- k+1
          for(j in k:2) covs[j-1,i] <- covs[j,i]
        }
        for(j in 2:nrow(trackData))
          if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
      }
    }
  }

  if(length(obsPerAnimal)==1)
    obsPerAnimal <- rep(obsPerAnimal,2)
  else if(length(obsPerAnimal)!=2)
    stop("obsPerAnimal should be of length 1 or 2.")

  #######################################
  ## Prepare parameters for simulation ##
  #######################################
  # define number of observations for each animal
  allNbObs <- rep(NA,nbAnimals)
  for(zoo in 1:nbAnimals) {
    if(obsPerAnimal[1]!=obsPerAnimal[2])
      allNbObs[zoo] <- sample(obsPerAnimal[1]:obsPerAnimal[2],size=1)
    else
      allNbObs[zoo] <- obsPerAnimal[1]
  }

  # extend covs if not enough covariate values
  if(!is.null(covs)) {
    covnames <- colnames(covs)
    while(sum(allNbObs)>nrow(covs))
      covs <- rbind(covs,covs)
    # shrink covs if too many covariate values
    covs <- data.frame(covs[1:sum(allNbObs),])
    colnames(covs) <- covnames
  }

  # generate regression parameters for transition probabilities
  if(is.null(beta))
    beta <- matrix(rnorm(nbStates*(nbStates-1)*(nbCovs+1)),nrow=nbCovs+1)
  else if(nrow(beta)!=nbCovs+1 | ncol(beta)!=nbStates*(nbStates-1)) {
    if(nbStates>1)
      stop(paste("beta should have ",nbCovs+1," rows and ",nbStates*(nbStates-1)," columns.",sep=""))
    else
      stop("beta should be NULL")
  }

  # initial state distribution
  delta <- rep(1,nbStates)/nbStates
  
  if(!all(unlist(lapply(p$bounds,is.numeric)))){
    bounds <- p$bounds
    for(i in distnames){
      if(!is.numeric(bounds[[i]])){
        bounds[[i]] <- gsub(i,"",bounds,fixed=TRUE)
      }
    }
  }
  
  # format parameters
  wpar <- n2w(lapply(Par,function(x) c(t(x))),bounds,beta,delta,nbStates,FALSE,DM,cons,logitcons) 
  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs,FALSE,stationary=FALSE,cons,DM,p$boundInd,logitcons)

  zeroMass<-vector('list',length(dist))
  names(zeroMass)<-distnames
  for(i in distnames){
    if(zeroInflation[[i]]) {
      zeroMass[[i]] <- par[[i]][nrow(par[[i]]),]
      par[[i]] <- par[[i]][-(nrow(par[[i]])),]
    }
    else {
      zeroMass[[i]] <- rep(0,nbStates)
    }
  }
  for(i in distnames[which(dist %in% c("wrpcauchy","vm"))]){
    par[[i]] <- rbind(rep(0,nbStates),par[[i]])
  }

  trackData <- NULL
  allCovs <- NULL
  allStates <- NULL

  # build the data frame to be returned
  eval(parse(text=paste0("data <- data.frame(ID=character(),x=numeric(),y=numeric(),",paste0(distnames,"=numeric()",collapse=","),")")))

  ###########################
  ## Loop over the animals ##
  ###########################
  for (zoo in 1:nbAnimals) {

    # number of observations for animal zoo
    nbObs <- allNbObs[zoo]

    ###############################
    ## Simulate covariate values ##
    ###############################
    if(nbCovs>0) {
      if(is.null(covs)) {
        subCovs <- data.frame(cov1=rnorm(nbObs))
        if(nbCovs>1) {
          for(j in 2:nbCovs) {
            c <- data.frame(rnorm(nbObs))
            colnames(c) <- paste("cov",j,sep="")
            subCovs <- cbind(subCovs,c)
          }
        }
      } else {
        # select covariate values which concern the current animal
        if(zoo<2)
          ind1 <- 1
        else
          ind1 <- sum(allNbObs[1:(zoo-1)])+1
        ind2 <- sum(allNbObs[1:zoo])
        subCovs <- data.frame(covs[ind1:ind2,])
        if(!is.null(covs))
          colnames(subCovs) <- colnames(covs) # keep covariates names from input
      }
      allCovs <- rbind(allCovs,subCovs)
    }

    ###############################
    ## Simulate state sequence Z ##
    ###############################
    if(nbStates>1) {
      Z <- rep(NA,nbObs)
      Z[1] <- sample(1:nbStates,size=1,prob=delta)
      for (k in 2:nbObs) {
        gamma <- diag(nbStates)

        g <- beta[1,]
        if(nbCovs==1) g <- g + beta[2,]*subCovs[k,1]
        if(nbCovs>1) {
          for(j in 1:nbCovs)
            g <- g + beta[j+1,]*subCovs[k,j]
        }

        gamma[!gamma] <- exp(g)
        gamma <- t(gamma)
        gamma <- gamma/apply(gamma,1,sum)
        Z[k] <- sample(1:nbStates,size=1,prob=gamma[Z[k-1],])
      }
      allStates <- c(allStates,Z)
    } else
      Z <- rep(1,nbObs)

    X <- matrix(0,nrow=nbObs,ncol=2)
    X[1,] <- c(0,0) # initial position of animal

    phi <- 0
    eval(parse(text=paste0(distnames," <- rep(NA,nbObs)")))
    
    if(all(c("step","angle") %in% distnames)){
      distnames<-c("step","angle",distnames[!(distnames %in% c("step","angle"))])
    }
    
    ############################
    ## Simulate movement path ##
    ############################
    for (k in 1:(nbObs-1)){
      # prepare lists of arguments for step and angle distributions
      for(i in distnames){
        genDist <- dist[[i]]
        genArgs <- list(1)  # first argument = 1 (one random draw)
        for(j in 1:nrow(par[[i]]))
          genArgs[[j+1]] <- par[[i]][j,Z[k]]


      if(genDist=="gamma") {
        shape <- genArgs[[2]]^2/genArgs[[3]]^2
        scale <- genArgs[[3]]^2/genArgs[[2]]
        genArgs[[2]] <- shape
        genArgs[[3]] <- 1/scale # rgamma expects rate=1/scale
      }

      if(runif(1)>zeroMass[[i]][Z[k]])
        eval(parse(text=paste0(i,"[k] <- do.call(Fun[[i]],genArgs)")))
      else
        eval(parse(text=paste0(i,"[k] <- 0")))

      if(i=="angle" & dist[[i]] %in% c("wrpcauchy","vm") & ("step" %in% distnames))
        if(dist[["step"]] %in% c("gamma","weibull","exp")) {
          if(step[k]>0){
            eval(parse(text=paste0(i,"[k] <- do.call(Fun[[i]],genArgs)")))
            eval(parse(text=paste0("if(",i,"[k] >  pi) ",i,"[k] <- ",i,"[k]-2*pi")))
            eval(parse(text=paste0("if(",i,"[k] < -pi) ",i,"[k] <- ",i,"[k]+2*pi")))
            eval(parse(text=paste0("phi <- phi + ",i,"[k]")))
          } else if(step[k]==0) {
            eval(parse(text=paste0(i,"[k] <- NA")))  # angle = NA if step = 0
          }
          m <- step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
          X[k+1,] <- X[k,] + m
        }
      }
    }

    eval(parse(text=paste0("d <- data.frame(ID=rep(zoo,nbObs),",paste0(distnames,"=",distnames,collapse=","),")")))
    if("angle" %in% distnames) 
      if(dist[["angle"]] %in% c("wrpcauchy","vm") & ("step" %in% distnames))
        if(dist[["step"]] %in% c("gamma","weibull","exp")){
            angle[1] <- NA # the first angle value is arbitrary
            eval(parse(text=paste0("d <- data.frame(ID=rep(zoo,nbObs),",paste0(distnames,"=",distnames,collapse=","),",x=X[,1],y=X[,2])")))
        }
    data <- rbind(data,d)
  }

  # if covs provided as argument
  if(!is.null(covs) & is.null(allCovs))
    allCovs <- covs

  if(nbCovs>0)
    data <- cbind(data,allCovs)

  # include states sequence in the data
  if(states)
    data <- cbind(data,states=allStates)
  return(moveData(data))
}
