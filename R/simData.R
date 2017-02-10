
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
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified if
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
#' @return An object momentuHMMData, i.e. a dataframe of:
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
#' @importFrom stats rnorm runif step terms.formula
#' @importFrom raster cellFromXY getValues
#' @importFrom moveHMM simData
#' @importFrom CircStats rvm

simData <- function(nbAnimals=1,nbStates=2,dist,
                    Par,beta=NULL,
                    formula=NULL,
                    covs=NULL,nbCovs=0,
                    spatialCovs=NULL,
                    zeroInflation=NULL,
                    circularAngleMean=NULL,
                    centers=NULL,
                    obsPerAnimal=c(500,1500),
                    DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL,stateNames=NULL,
                    model=NULL,states=FALSE,
                    lambda=NULL,
                    errorEllipse=NULL)
{
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    # extract simulation parameters from model
    nbStates <- length(model$stateNames)
    dist<-model$conditions$dist
    distnames<-names(dist)
    userBounds <- model$conditions$bounds
    stateNames<-model$stateNames
    estAngleMean<-model$conditions$estAngleMean
    DM <- model$conditions$DM
    cons <- model$conditions$cons
    workcons <- model$conditions$workcons
    zeroInflation <- model$conditions$zeroInflation
    formula <- model$conditions$formula
  
    Par <- model$mle[distnames]
    parindex <- c(0,cumsum(unlist(lapply(model$conditions$fullDM,ncol)))[-length(model$conditions$fullDM)])
    names(parindex) <- distnames
    for(i in distnames){
      if(!is.null(DM[[i]]) & model$conditions$DMind[[i]]){
        Par[[i]] <- model$mod$estimate[parindex[[i]]+1:ncol(model$conditions$fullDM[[i]])]
        names(Par[[i]])<-colnames(model$conditions$fullDM[[i]])
      }
    }
    for(i in distnames[which(dist %in% angledists)]){
      if(!estAngleMean[[i]]){
        estAngleMean[[i]]<-TRUE
        userBounds[[i]]<-rbind(matrix(rep(c(-pi,pi),nbStates),nbStates,2,byrow=TRUE),userBounds[[i]])
        cons[[i]] <- c(rep(1,nbStates),cons[[i]])
        workcons[[i]] <- c(rep(0,nbStates),workcons[[i]])
        if(!is.null(DM[[i]])){
          Par[[i]] <- c(rep(0,nbStates),Par[[i]])
          if(is.list(DM[[i]])){
            DM[[i]]$mean<- ~1
          } else {
            tmpDM <- matrix(0,nrow(DM[[i]])+nbStates,ncol(DM[[i]])+nbStates)
            tmpDM[nbStates+1:nrow(DM[[i]]),nbStates+1:ncol(DM[[i]])] <- DM[[i]]
            diag(tmpDM)[1:nbStates] <- 1
            DM[[i]] <- tmpDM
          }
        }
      }
    }
    beta <- model$mle$beta
    delta <- model$mle$delta
    Par<-lapply(Par,function(x) c(t(x)))
    
    if(states) model$data$states <- NULL

    if(is.null(covs)) {
      if(!is.null(spatialCovs)) spatialcovnames <- names(spatialCovs)
      else spatialcovnames <- NULL
      covsCol <- seq(1,ncol(model$data))[-match(c("ID","x","y",distnames,spatialcovnames),names(model$data),nomatch=0)]
      #covs <- model.matrix(model$conditions$formula,model$data)

      if(length(covsCol)) {
        covs <- model$data[,covsCol,drop=FALSE]
      }
    }
    # else, allow user to enter new values for covariates

  } else {
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
    distnames<-names(dist)
    if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
    Par <- Par[distnames]
    delta <- NULL
    
    mHind <- (is.null(DM) & is.null(userBounds) & is.null(spatialCovs) & ("step" %in% names(dist)) & is.null(lambda) & is.null(errorEllipse)) # indicator for moveHMM::simData
    if(all(names(dist) %in% c("step","angle")) & mHind){
      zi <- FALSE
      if(!is.null(zeroInflation$step)) zi <- zeroInflation$step
      if(is.null(dist$angle)) dist$angle<-"none"
      data <- moveHMM::simData(nbAnimals, nbStates, dist$step, dist$angle, Par$step, Par$angle, beta, covs, nbCovs, zi, obsPerAnimal, model, states)
      attr(data,"class") <- "data.frame"
      data$ID <- as.factor(data$ID)
      return(momentuHMMData(data))
    }
  }
  
  Fun <- lapply(dist,function(x) paste("r",x,sep=""))

  if(nbAnimals<1)
    stop("nbAnimals should be at least 1.")
  if(nbStates<1)
    stop("nbStates should be at least 1.")
  
  if(is.null(zeroInflation)){
    zeroInflation <- vector('list',length(distnames))
    names(zeroInflation) <- distnames
    for(i in distnames){
      zeroInflation[[i]]<-FALSE
    }
  } else {
    if(!is.list(zeroInflation) | is.null(names(zeroInflation))) stop("'zeroInflation' must be a named list")
    for(i in distnames){
      if(is.null(zeroInflation[[i]])) zeroInflation[[i]] <- FALSE
    }
  }

  if(!all(unlist(lapply(zeroInflation,is.logical)))) stop("zeroInflation must be a list of logical objects")
  for(i in distnames){
    if((dist[[i]] %in% angledists | dist[[i]]=="pois") & zeroInflation[[i]])
      stop(dist[[i]]," distribution cannot be zero inflated")
  }
  
  estAngleMean <- vector('list',length(distnames))
  names(estAngleMean) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% angledists) estAngleMean[[i]]<-TRUE
    else estAngleMean[[i]]<-FALSE
  }
  
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,DM,userBounds,cons,workcons,stateNames)
  p <- inputs$p
  parSize <- p$parSize
  bounds <- p$bounds

  spatialcovnames<-NULL
  if(!is.null(spatialCovs)){
    if(!all(c("step","angle") %in% distnames)) stop("spatialCovs can only be included when 'step' and 'angle' distributions are specified") 
    else if(!(dist[["angle"]] %in% angledists) | !(dist[["step"]] %in% stepdists)) stop("spatialCovs can only be included when valid 'step' and 'angle' distributions are specified") 
    nbSpatialCovs<-length(names(spatialCovs))
    for(j in 1:nbSpatialCovs){
      if(class(spatialCovs[[j]])!="RasterLayer") stop("spatialCovs must be of class 'RasterLayer'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs")
    }
    spatialcovnames<-names(spatialCovs)
  } else nbSpatialCovs <- 0

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
        for(j in 2:nrow(covs))
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
  cumNbObs <- c(0,cumsum(allNbObs))

  # extend covs if not enough covariate values
  if(!is.null(covs)) {
    covnames <- colnames(covs)
    while(sum(allNbObs)>nrow(covs))
      covs <- rbind(covs,covs)
    # shrink covs if too many covariate values
    covs <- data.frame(covs[1:sum(allNbObs),])
    colnames(covs) <- covnames
  }
  
  ###############################
  ## Simulate covariate values ##
  ###############################
  allCovs <- NULL
  if(nbCovs>0) {
    if(is.null(covs)) {
      allCovs <- data.frame(cov1=rnorm(sum(allNbObs)))
      if(nbCovs>1) {
        for(j in 2:nbCovs) {
          c <- data.frame(rnorm(sum(allNbObs)))
          colnames(c) <- paste("cov",j,sep="")
          allCovs <- cbind(allCovs,c)
        }
      }
    } else {
      allCovs <- covs
    }
  }
  
  if(anyDuplicated(colnames(allCovs))) stop("covariates must have unique names")
  if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  if(any(colnames(allCovs) %in% spatialcovnames)) stop("spatialCovs name(s) cannot match other covariate name(s)")
  
  if(!is.null(centers)){
    if(any(dim(centers)!=c(nbStates,2))) stop("centers must be of dimension ",nbStates,"x",2)
    centerInd <- which(!apply(centers,1,function(x) any(is.na(x))))
    if(length(centerInd)){
      centerNames<-paste0("center",".",rep(c("dist","angle"),length(centerInd)),rep(centerInd,each=2))
      centerCovs <- data.frame(matrix(NA,nrow=sum(allNbObs),ncol=length(centerInd)*2,dimnames=list(NULL,centerNames)))
    }  
  }
  
  allNbCovs <- nbCovs+nbSpatialCovs
  if(is.null(formula)) {
    if(allNbCovs) formula <- formula(paste0("~",paste0(c(colnames(allCovs),spatialcovnames),collapse="+")))
    else formula <- formula(~1)
  }
  if(is.null(beta))
    beta <- matrix(rnorm(nbStates*(nbStates-1)*(length(attr(terms.formula(formula),"term.labels"))+1)),nrow=length(attr(terms.formula(formula),"term.labels"))+1)
  else {
    if(ncol(beta)!=nbStates*(nbStates-1) | nrow(beta)!=length(attr(terms.formula(formula),"term.labels"))+1) {
      error <- paste("beta has wrong dimensions: it should have",length(attr(terms.formula(formula),"term.labels"))+1,"rows and",
                     nbStates*(nbStates-1),"columns.")
      stop(error)
    }
  }

  # initial state distribution
  if(is.null(delta)) delta <- rep(1,nbStates)/nbStates

  zeroMass<-vector('list',length(dist))
  names(zeroMass)<-distnames

  allStates <- NULL
  allSpatialcovs<-NULL
  
  #make sure 'step' preceeds 'angle'
  if(all(c("step","angle") %in% distnames)){
    distnames<-c("step","angle",distnames[!(distnames %in% c("step","angle"))])
  }
  
  # build the data frame to be returned
  data<-data.frame(ID=factor())
  for(i in distnames){
    data[[i]]<-numeric()
  }
  if("angle" %in% distnames) 
    if(dist[["angle"]] %in% angledists & ("step" %in% distnames))
      if(dist[["step"]] %in% stepdists){
        data$x<-numeric()
        data$y<-numeric()
      }
  
  ###########################
  ## Loop over the animals ##
  ###########################
  for (zoo in 1:nbAnimals) {

    # number of observations for animal zoo
    nbObs <- allNbObs[zoo]
    d <- data.frame(ID=factor(rep(zoo,nbObs)))
    #subPar <- lapply(par[distnames],function(x) x[,cumNbObs[zoo]+1:nbObs])
    
    ###############################
    ## Simulate covariate values ##
    ###############################
    subCovs<-as.data.frame(matrix(NA,nrow=nbObs,ncol=nbCovs))
    if(nbCovs>0) {
      # select covariate values which concern the current animal
      if(zoo<2)
        ind1 <- 1
      else
        ind1 <- sum(allNbObs[1:(zoo-1)])+1
        ind2 <- sum(allNbObs[1:zoo])
        subCovs <- data.frame(allCovs[ind1:ind2,,drop=FALSE])
        if(!is.null(covs))
          colnames(subCovs) <- colnames(covs) # keep covariates names from input
    }
    if(length(centerInd)) subCovs <- cbind(subCovs,centerCovs[cumNbObs[zoo]+1:nbObs,])
    
    subSpatialcovs<-as.data.frame(matrix(NA,nrow=nbObs,ncol=nbSpatialCovs))
    colnames(subSpatialcovs)<-spatialcovnames

    X <- matrix(0,nrow=nbObs,ncol=2)
    X[1,] <- c(0,0) # initial position of animal

    phi <- 0
    
    ############################
    ## Simulate movement path ##
    ############################
 
    genData <- genArgs <- vector('list',length(distnames))
    names(genData) <- names(genArgs) <- distnames
    
    for(i in distnames){
      genData[[i]] <- rep(NA,nbObs)
      genArgs[[i]] <- list(1)  # first argument = 1 (one random draw)
    }
    
    gamma <- diag(nbStates)
    
    if(!nbSpatialCovs & !length(centerInd)) {
      DMcov <- model.matrix(formula,subCovs)
      gFull <-  DMcov %*% beta
      
      # format parameters
      DMinputs<-getDM(subCovs,inputs$DM,dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,inputs$circularAngleMean)
      fullDM <- DMinputs$fullDM
      DMind <- DMinputs$DMind
      wpar <- n2w(Par,bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
      fullsubPar <- w2n(wpar,bounds,parSize,nbStates,length(attr(terms.formula(formula),"term.labels")),inputs$estAngleMean,inputs$circularAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nbObs,dist,p$Bndind)
      g <- gFull[1,,drop=FALSE]
    } else {
      
      if(nbSpatialCovs){
        for(j in 1:nbSpatialCovs){
          getCell<-raster::cellFromXY(spatialCovs[[j]],c(X[1,1],X[1,2]))
          if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
          subSpatialcovs[1,j]<-spatialCovs[[j]][getCell]
        }
      }
      
      if(length(centerInd)){
        for(j in 1:length(centerInd)){
          subCovs[1,centerNames[(j-1)*length(centerInd)+1:2]]<-distAngle(X[1,],X[1,],centers[centerInd[j],])
        }
      }
      g <- model.matrix(formula,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% beta
    }
    gamma[!gamma] <- exp(g)
    gamma <- t(gamma)
    gamma <- gamma/apply(gamma,1,sum)
    
    if(nbStates>1) {
      Z <- rep(NA,nbObs)
      Z[1] <- sample(1:nbStates,size=1,prob=delta%*%gamma)
    } else
      Z <- rep(1,nbObs)
    
    for (k in 1:(nbObs-1)){
      
      if(nbSpatialCovs |  length(centerInd)){
        # format parameters
        DMinputs<-getDM(cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE]),inputs$DM,dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par,bounds,beta,delta,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
        subPar <- w2n(wpar,bounds,parSize,nbStates,length(attr(terms.formula(formula),"term.labels")),inputs$estAngleMean,inputs$circularAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,dist,p$Bndind)
      } else {
        subPar <- lapply(fullsubPar[distnames],function(x) x[,k,drop=FALSE])#fullsubPar[,k,drop=FALSE]
      }
      
      for(i in distnames){
        
        if(zeroInflation[[i]]) {
          zeroMass[[i]] <- subPar[[i]][parSize[[i]]*nbStates-(nbStates-1):0]
          subPar[[i]] <- subPar[[i]][-(parSize[[i]]*nbStates-(nbStates-1):0)]
        }
        else {
          zeroMass[[i]] <- rep(0,nbStates)
        }

        for(j in 1:(parSize[[i]]-zeroInflation[[i]]))
          genArgs[[i]][[j+1]] <- subPar[[i]][(j-1)*nbStates+Z[k]]
        
        if(dist[[i]] %in% angledists){
          
          genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
          if(genData[[i]][k] >  pi) genData[[i]][k] <- genData[[i]][k]-2*pi
          if(genData[[i]][k] < -pi) genData[[i]][k] <- genData[[i]][k]+2*pi

          if(i=="angle" & ("step" %in% distnames)){
            if(dist[["step"]] %in% stepdists) {
              if(genData$step[k]>0){
                phi <- phi + genData[[i]][k]
              } else if(genData$step[k]==0) {
                genData[[i]][k] <- NA # angle = NA if step = 0
                #if(length(centerInd)) subCovs[k,centerNames[seq(2,2*length(centerInd),2)]] <- NA
              }
              m <- genData$step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
              X[k+1,] <- X[k,] + m
            }
          }
        } else {
          
          if(dist[[i]]=="gamma") {
            shape <- genArgs[[i]][[2]]^2/genArgs[[i]][[3]]^2
            scale <- genArgs[[i]][[3]]^2/genArgs[[i]][[2]]
            genArgs[[i]][[2]] <- shape
            genArgs[[i]][[3]] <- 1/scale # rgamma expects rate=1/scale
          }
    
          if(runif(1)>zeroMass[[i]][Z[k]])
            genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
          else
            genData[[i]][k] <- 0
        }
        
        d[[i]] <- genData[[i]]
        
      }
      # get next state
      gamma <- diag(nbStates)
      if(nbSpatialCovs | length(centerInd)){
        if(nbSpatialCovs){
          for(j in 1:nbSpatialCovs){
            getCell<-raster::cellFromXY(spatialCovs[[j]],c(X[k+1,1],X[k+1,2]))
            if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
            subSpatialcovs[k+1,j]<-spatialCovs[[j]][getCell]
          }
        }
        if(length(centerInd)){
          for(j in 1:length(centerInd)){
            subCovs[k+1,centerNames[(j-1)*length(centerInd)+1:2]]<-distAngle(X[k,],X[k+1,],centers[centerInd[j],])
          }
        }
        g <- model.matrix(formula,cbind(subCovs[k+1,,drop=FALSE],subSpatialcovs[k+1,,drop=FALSE])) %*% beta
      } else {
        g <- gFull[k+1,,drop=FALSE]
      }
      gamma[!gamma] <- exp(g)
      gamma <- t(gamma)
      gamma <- gamma/apply(gamma,1,sum)
      Z[k+1] <- sample(1:nbStates,size=1,prob=gamma[Z[k],])  
    }
    allStates <- c(allStates,Z)
    if(nbSpatialCovs>0) {
      allSpatialcovs <- rbind(allSpatialcovs,subSpatialcovs)
    }
    
    if("angle" %in% distnames) 
      if(dist[["angle"]] %in% angledists & ("step" %in% distnames))
        if(dist[["step"]] %in% stepdists){
            d$angle[1] <- NA # the first angle value is arbitrary
            #if(length(centerInd)) subCovs[1,centerNames[seq(2,2*length(centerInd),2)]] <- NA
            d$x=X[,1]
            d$y=X[,2]
        }
    
    if(length(centerInd)) centerCovs[cumNbObs[zoo]+1:nbObs,] <- subCovs[,centerNames]
    data <- rbind(data,d)
  }
  
  if(nbSpatialCovs>0) colnames(allSpatialcovs)<-spatialcovnames

  if(nbCovs>0)
    data <- cbind(data,allCovs)
  
  if(nbSpatialCovs>0)
    data <- cbind(data,allSpatialcovs)
  
  if(length(centerInd))
    data <- cbind(data,centerCovs)
  
  # include states sequence in the data
  if(states)
    data <- cbind(data,states=allStates)
  
  if(!is.null(lambda) | !is.null(errorEllipse)){
    if(!is.null(errorEllipse)){
      if(!is.list(errorEllipse) | any(!(c("M","m","r") %in% names(errorEllipse)))) stop("errorEllipse must be a list of scalars named 'M', 'm', and 'r'.")
      if(any(unlist(lapply(errorEllipse[c("M","m","r")],length))>1)) stop('errorEllipse must consist of positive scalars')
      if(any(unlist(lapply(errorEllipse[c("M","m")],function(x) x<0)))) stop("errorEllipse$M and errorEllipse$m must be >=0")
      if(errorEllipse$r < 0 | errorEllipse$r > 180) stop('errorEllipse$r must be in [0,180]')
    }
    if(!is.null(lambda))
      if(lambda<=0) stop('lambda must be >0')
    
    # account for location measurement error and/or temporal irregularity
    return(simObsData(data,dist,lambda,errorEllipse))
    
  } else {
    
    return(momentuHMMData(data))
    
  }
}
