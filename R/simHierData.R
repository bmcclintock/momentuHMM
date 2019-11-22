
#' @rdname simData
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. See details.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (only hierarchical-level effects, with no covariate effects).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy. See details.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' @param nbHierCovs A hierarchical data structure \code{\link[data.tree]{Node}} for the number of covariates ('nbCovs') to simulate for each level of the hierarchy (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1.1' and 'cov1.2' for \code{nbHierCovs$level1$nbCovs=2}) and can be included in \code{hierFormula} and/or \code{DM}.
#' @param obsPerLevel A hierarchical data structure \code{\link[data.tree]{Node}} indicating the number of observations for each level of the hierarchy ('obs'). For each level, the 'obs' field can either be the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated per level for each animal are uniformously picked from this interval. Alternatively, \code{obsPerLevel} can be specified as
#' a list of length \code{nbAnimals} with each element providing the hierarchical data structure for the number of observations for each level of the hierarchy for each animal, where the 'obs' field can either be the number of observations (if single value) or the bounds of the number of observations (if vector of two values) for each individual.
#'
#' @details \itemize{
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal' (or 'obsPerLevel' for \code{simHierData}), then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' 
#' \item For \code{simData}, when covariates are not included in \code{formulaDelta} (i.e. \code{formulaDelta=NULL}), then \code{delta} is specified as a vector of length \code{nbStates} that 
#' sums to 1.  When covariates are included in \code{formulaDelta}, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#' 
#' \item For \code{simHierData}, \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. 
#' }
#' 
#' @references
#' 
#' Leos-Barajas, V., Gangloff, E.J., Adam, T., Langrock, R., van Beest, F.M., Nabe-Nielsen, J. and Morales, J.M. 2017. 
#' Multi-scale modeling of animal movement and general behavior data using hidden Markov models with hierarchical structures. 
#' Journal of Agricultural, Biological and Environmental Statistics, 22 (3), 232-248.
#'
#' @export
#' @importFrom stats rnorm runif rmultinom step terms.formula
#' @importFrom raster cellFromXY getValues
#' @importFrom CircStats rvm
#' @importFrom Brobdingnag as.brob sum
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.tree Node Get Aggregate isLeaf Clone

simHierData <- function(nbAnimals=1,hierStates,hierDist,
                    Par,hierBeta=NULL,hierDelta=NULL,
                    hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,
                    covs=NULL,nbHierCovs=NULL,
                    spatialCovs=NULL,
                    zeroInflation=NULL,
                    oneInflation=NULL,
                    circularAngleMean=NULL,
                    centers=NULL,
                    centroids=NULL,
                    angleCovs=NULL,
                    obsPerLevel,
                    initialPosition=c(0,0),
                    DM=NULL,userBounds=NULL,workBounds=NULL,mvnCoords=NULL,
                    model=NULL,states=FALSE,
                    retrySims=0,
                    lambda=NULL,
                    errorEllipse=NULL)
{
  
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    
    if(!inherits(model,"momentuHierHMM") & !inherits(model,"hierarchical")) stop("model must be a 'momentuHierHMM' and/or 'hierarchical' object")
    
    if(is.miHMM(model)){
      model <- model$miSum
    }
    
    model <- delta_bc(model)
    
    # extract simulation parameters from model
    nbStates <- length(model$stateNames)
    dist<-model$conditions$dist
    distnames<-names(dist)
    
    if(is.miSum(model)){
      model <- formatmiSum(model)
      if(!is.null(model$mle$beta)) model$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(model$mle$beta),2,byrow=TRUE)
      if(!is.null(model$Par$beta$pi$est)) model$conditions$workBounds$pi<-matrix(c(-Inf,Inf),length(model$Par$beta$pi$est),2,byrow=TRUE)
      if(!is.null(model$Par$beta$delta$est)) model$conditions$workBounds$delta<-matrix(c(-Inf,Inf),length(model$Par$beta$delta$est),2,byrow=TRUE)
      if(!is.null(model$mle$g0)) model$conditions$workBounds$g0<-matrix(c(-Inf,Inf),length(model$mle$g0),2,byrow=TRUE)
      if(!is.null(model$mle$theta)) model$conditions$workBounds$theta<-matrix(c(-Inf,Inf),length(model$mle$theta),2,byrow=TRUE)
    }
    
    hierStates <- model$conditions$hierStates
    hierDist <- model$conditions$hierDist
    hierFormula <- model$conditions$hierFormula
    hierFormulaDelta <- model$conditions$hierFormulaDelta
    hierBeta <- model$conditions$hierBeta
    hierDelta <- model$conditions$hierDelta
    userBounds <- model$conditions$bounds
    workBounds <- model$conditions$workBounds
    mvnCoords <- model$conditions$mvnCoords
    stateNames<-model$stateNames
    estAngleMean<-model$conditions$estAngleMean
    circularAngleMean<-model$conditions$circularAngleMean
    DM <- model$conditions$DM
    cons <- model$conditions$cons
    workcons <- model$conditions$workcons
    betaRef <- model$conditions$betaRef
    zeroInflation <- model$conditions$zeroInflation
    oneInflation <- model$conditions$oneInflation
    formula <- model$conditions$formula
    formulaDelta <- formDelta <- model$condition$formulaDelta
    mixtures <- model$conditions$mixtures
    if(is.null(model$condition$formulaPi)){
      formulaPi <- formPi <- ~1
    } else formulaPi <- formPi <- model$condition$formulaPi
    
    Par <- model$mle[distnames]
    parCount<- lapply(model$conditions$fullDM,ncol)
    for(i in distnames[!unlist(lapply(circularAngleMean,isFALSE))]){
      parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(model$conditions$fullDM[[i]])))))
    }
    parindex <- c(0,cumsum(unlist(parCount))[-length(model$conditions$fullDM)])
    names(parindex) <- distnames
    for(i in distnames){
      if(!is.null(DM[[i]])){
        Par[[i]] <- model$mod$estimate[parindex[[i]]+1:parCount[[i]]]
        if(!isFALSE(circularAngleMean[[i]])){
          names(Par[[i]]) <- unique(gsub("cos","",gsub("sin","",colnames(model$conditions$fullDM[[i]]))))
        } else names(Par[[i]])<-colnames(model$conditions$fullDM[[i]])
        #cons[[i]]<-rep(1,length(cons[[i]]))
        #workcons[[i]]<-rep(0,length(workcons[[i]]))
      }
    }
    for(i in distnames[which(dist %in% angledists)]){
      if(!estAngleMean[[i]]){
        estAngleMean[[i]]<-TRUE
        userBounds[[i]]<-rbind(matrix(rep(c(-pi,pi),nbStates),nbStates,2,byrow=TRUE),userBounds[[i]])
        workBounds[[i]]<-rbind(matrix(rep(c(-Inf,Inf),nbStates),nbStates,2,byrow=TRUE),workBounds[[i]])
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
    
    if(!is.null(model$conditions$recharge)){
      g0 <- nw2w(model$mle$g0,model$conditions$workBounds$g0)
      theta <- nw2w(model$mle$theta,model$conditions$workBounds$theta)
    } else g0 <- theta <- NULL
    beta <- nw2w(model$mle$beta,model$conditions$workBounds$beta)
    nbCovsDelta <- ncol(model$covsDelta)-1
    foo <- length(model$mod$estimate)-length(g0)-length(theta)-(nbCovsDelta+1)*(nbStates-1)*mixtures+1
    delta <- matrix(model$mod$estimate[foo:(length(model$mod$estimate)-length(g0)-length(theta))],nrow=(nbCovsDelta+1)*mixtures) 
    if(mixtures>1) {
      nbCovsPi <- ncol(model$covsPi)-1
      foo <- length(model$mod$estimate)-length(g0)-length(theta)-(nbCovsDelta+1)*(nbStates-1)*mixtures-(nbCovsPi+1)*(mixtures-1)+1:((nbCovsPi+1)*(mixtures-1))
      pie <- matrix(model$mod$estimate[foo],nrow=nbCovsPi+1,ncol=mixtures-1)
      #pie <- model$mle$pi
      #workBounds$pi <- NULL
    } else {
      pie <- NULL
      nbCovsPi <- 0
    }
    beta <- list(beta=beta,pi=pie,g0=g0,theta=theta)
    
    workBounds$beta <- model$conditions$hierBeta
    workBounds$delta <- model$conditions$hierDelta
    
    Par<-lapply(Par,function(x) c(t(x)))
    
    if(states) model$data$states <- NULL
    
    if(is.null(covs)) {
      p<-parDef(lapply(dist,function(x) gsub("Consensus","",x)),nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds)
      covNames<-character()
      for(i in distnames){
        covNames<-c(covNames,getCovNames(model,p,i)$DMterms)
      }
      if(!is.null(model$rawCovs)){
        covNames <- c(colnames(model$rawCovs),covNames)
      }
      covNames <- c(covNames,colnames(model$covsPi)[which(colnames(model$covsPi)!="(Intercept)")],colnames(model$covsDelta)[which(colnames(model$covsDelta)!="(Intercept)")])
      covsCol <- unique(covNames)
      factorterms<-names(model$data)[unlist(lapply(model$data,is.factor))]
      factorcovs<-paste0(rep(factorterms,times=unlist(lapply(model$data[factorterms],nlevels))),unlist(lapply(model$data[factorterms],levels)))
      
      if(length(covsCol)){
        for(jj in 1:length(covsCol)){
          cov<-covsCol[jj]
          form<-formula(paste("~",cov))
          varform<-all.vars(form)
          if(any(varform %in% factorcovs) && !all(varform %in% factorterms)){
            factorvar<-factorcovs %in% (varform[!(varform %in% factorterms)])
            covsCol[jj]<-rep(factorterms,times=unlist(lapply(model$data[factorterms],nlevels)))[which(factorvar)]
          } 
        }
      }
      covsCol<-unique(covsCol)
      covsCol <- covsCol[!(covsCol %in% c("ID","level"))]
      
      if(length(covsCol)) covs <- model$data[covsCol]
      
      if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% rwdists){
        covs[[paste0(mvnCoords,".x_tm1")]] <- NULL
        covs[[paste0(mvnCoords,".y_tm1")]] <- NULL
        if(dist[[mvnCoords]]=="rw_mvnorm3") covs[[paste0(mvnCoords,".z_tm1")]] <- NULL
        if(!ncol(covs)) covs <- NULL
      }
    }
    # else, allow user to enter new values for covariates
    
  } else {
    
    inputHierHMM <- formatHierHMM(data=NULL,hierStates,hierDist,hierBeta,hierDelta,hierFormula,hierFormulaDelta,mixtures,workBounds)
    
    nbStates <- inputHierHMM$nbStates
    dist <- inputHierHMM$dist
    beta <- inputHierHMM$beta
    delta <- inputHierHMM$delta
    formula <- inputHierHMM$formula
    formulaDelta <- inputHierHMM$formulaDelta
    betaRef <- inputHierHMM$betaRef
    stateNames <- inputHierHMM$stateNames
    
    #if(!is.null(workBounds$beta)) stop("'workBounds$beta' cannot be specified; use 'hierBeta' instead")
    #if(!is.null(workBounds$delta)) stop("'workBounds$delta' cannot be specified; use 'hierDelta' instead")
    
    cons <- workcons <- NULL
    
    distnames <- names(dist)
    
    if(is.null(formulaPi)){
      formPi <- ~1
    } else formPi <- formulaPi
    formDelta <- formulaDelta
  }
  
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
  if(is.null(oneInflation)){
    oneInflation <- vector('list',length(distnames))
    names(oneInflation) <- distnames
    for(i in distnames){
      oneInflation[[i]]<-FALSE
    }
  } else {
    if(!is.list(oneInflation) | is.null(names(oneInflation))) stop("'oneInflation' must be a named list")
    for(i in distnames){
      if(is.null(oneInflation[[i]])) oneInflation[[i]] <- FALSE
    }
  }
  
  if(!all(unlist(lapply(zeroInflation,is.logical)))) stop("zeroInflation must be a list of logical objects")
  if(!all(unlist(lapply(oneInflation,is.logical)))) stop("oneInflation must be a list of logical objects")
  for(i in distnames){
    if(!(dist[[i]] %in% zeroInflationdists) & zeroInflation[[i]])
      stop(dist[[i]]," distribution cannot be zero inflated")
    if(!(dist[[i]] %in% oneInflationdists) & oneInflation[[i]])
      stop(dist[[i]]," distribution cannot be one inflated")
  }
  
  if(!is.null(mvnCoords)){
    if(length(mvnCoords)>1 | !is.character(mvnCoords)) stop("mvnCoords must be a character string")
    if(!(mvnCoords %in% distnames)) stop("mvnCoords not found. Permitted values are: ",paste0(distnames,collapse=", "))
    if(!(dist[[mvnCoords]] %in% mvndists)) stop("mvnCoords must correspond to a multivariate normal data stream")
    if(any(c("step","angle") %in% distnames)) stop("step and angle distributions cannot be specified when ",mvnCoords," ~ ",dist[[mvnCoords]])
    if(mvnCoords %in% c("x","y","z")) stop("'x', 'y', and 'z' are reserved and cannot be used for mvnCoords data stream names")
    if(sum(unlist(dist) %in% rwdists)>1) stop("sorry, simHierData currently only supports a single multivariate normal random walk distribution (and it must correspond to location data)")
  } else if(any(unlist(dist) %in% rwdists)) stop("sorry, simHierData currently only supports random walk distributions for multivariate location data identified through the mvnCoords argument")
  if(any(unlist(dist)=="rw_norm")) stop("sorry, 'rw_norm' is not currently supported by simHierData")
  
  estAngleMean <- vector('list',length(distnames))
  names(estAngleMean) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% angledists) estAngleMean[[i]]<-TRUE
    else estAngleMean[[i]]<-FALSE
  }
  
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,cons,workcons,stateNames,checkInflation = TRUE)
  p <- inputs$p
  parSize <- p$parSize
  bounds <- p$bounds
  
  Fun <- lapply(inputs$dist,function(x) paste("r",x,sep=""))
  
  spatialcovnames<-NULL
  if(!is.null(spatialCovs)){
    if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
    spatialcovnames<-names(spatialCovs)
    if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
    nbSpatialCovs<-length(spatialcovnames)
    if(is.null(mvnCoords)){
      if(!("step" %in% distnames)) stop("spatialCovs can only be included when 'step' distribution is specified") 
      else if(!(inputs$dist[["step"]] %in% stepdists)) stop("spatialCovs can only be included when valid 'step' distributions are specified") 
    }
    for(j in 1:nbSpatialCovs){
      if(!inherits(spatialCovs[[j]],c("RasterLayer","RasterBrick","RasterStack"))) stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer', 'RasterStack', or 'RasterBrick'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs$",spatialcovnames[j])
      if(inherits(spatialCovs[[j]],c("RasterBrick","RasterStack"))){
        if(is.null(raster::getZ(spatialCovs[[j]]))) stop("spatialCovs$",spatialcovnames[j]," is a raster stack or brick that must have set z values (see ?raster::setZ)")
        else if(!(names(attributes(spatialCovs[[j]])$z) %in% names(covs))) {
          if(!is.null(model)) covs[[names(attributes(spatialCovs[[j]])$z)]] <- model$data[[names(attributes(spatialCovs[[j]])$z)]]
          else stop("spatialCovs$",spatialcovnames[j]," z value '",names(attributes(spatialCovs[[j]])$z),"' not found in covs")
        }
        zname <- names(attributes(spatialCovs[[j]])$z)
        zvalues <- raster::getZ(spatialCovs[[j]])
        if(!all(unique(covs[[zname]]) %in% zvalues)) stop("data$",zname," includes z-values with no matching raster layer in spatialCovs$",spatialcovnames[j])
      }
    }
  } else nbSpatialCovs <- 0
  
  if(is.list(obsPerLevel)){
    if(length(obsPerLevel)!=nbAnimals) stop("obsPerLevel must be a list of length ",nbAnimals," where each element is a hierarchical Node")
    for(i in 1:length(obsPerLevel)){
      if(!inherits(obsPerLevel[[i]],"Node")) stop("obsPerLevel for individual ",i," must be a hierarchical Node")
      if(!("obs" %in% obsPerLevel[[i]]$fieldsAll)) stop("'obsPerLevel' for individual ",i," must include a 'obs' field")
      if(!all(obsPerLevel[[i]]$Get("name",filterFun=function(x) x$level==2) %in% hierDist$Get("name",filterFun=function(x) x$level==2))) stop("'obsPerLevel' level types for individual ",i," can only include ",paste0(hierDist$Get("name",filterFun=function(x) x$level==2),collapse = ", "))
      #if(!all(obsPerLevel[[i]]$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",1:obsPerLevel[[i]]$count))) stop("'obsPerLevel' level types for individual ",i," can only include ",paste0(paste0("level",1:obsPerLevel[[i]]$count),collapse = ", "))
      tmpObs <- data.tree::Clone(obsPerLevel[[i]])
      for(j in obsPerLevel[[i]]$Get("name",filterFun=function(x) x$level==2)){
        if(length(which(obsPerLevel[[i]][[j]]$obs<1))>0)
          stop("obsPerLevel 'obs' field should have positive values.")
        if(length(obsPerLevel[[i]][[j]]$obs)==1)
          tmpObs[[j]]$obs <- rep(obsPerLevel[[i]][[j]]$obs,2)
        else if(length(tmpObs[[j]]$obs)!=2)
          stop("obsPerLevel 'obs' field should be of length 1 or 2.")
      }
      obsPerLevel[[i]] <- tmpObs
    }
  } else {
    if(!inherits(obsPerLevel,"Node")) stop("obsPerLevel must be a hierarchical Node")
    if(!("obs" %in% obsPerLevel$fieldsAll)) stop("'obsPerLevel' must include a 'obs' field")
    if(!all(obsPerLevel$Get("name",filterFun=function(x) x$level==2) %in% hierDist$Get("name",filterFun=function(x) x$level==2))) stop("'obsPerLevel' level types can only include ",paste0(hierDist$Get("name",filterFun=function(x) x$level==2),collapse = ", "))
    #if(!all(obsPerLevel$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",1:obsPerLevel$count))) stop("'obsPerLevel' level types can only include ",paste0(paste0("level",1:obsPerLevel[[i]]$count),collapse = ", "))
    tmpObs <- data.tree::Clone(obsPerLevel)
    for(j in obsPerLevel$Get("name",filterFun=function(x) x$level==2)){
      if(length(which(obsPerLevel[[j]]$obs<1))>0)
        stop("obsPerLevel 'obs' field should have positive values.")
      if(length(obsPerLevel[[j]]$obs)==1)
        tmpObs[[j]]$obs <- rep(obsPerLevel[[j]]$obs,2)
      else if(length(tmpObs[[j]]$obs)!=2)
        stop("obsPerLevel 'obs' field should be of length 1 or 2.")
    }
    obsPerLevel<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      obsPerLevel[[i]]<-tmpObs
    }
  }
  
  if(is.list(initialPosition)){
    if(length(initialPosition)!=nbAnimals) stop("initialPosition must be a list of length ",nbAnimals)
    for(i in 1:nbAnimals){
      if(is.null(mvnCoords) || dist[[mvnCoords]] %in% c("mvnorm2","rw_mvnorm2")){
        if(length(initialPosition[[i]])!=2 | !is.numeric(initialPosition[[i]]) | any(!is.finite(initialPosition[[i]]))) stop("each element of initialPosition must be a finite numeric vector of length 2")
      } else if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")){
        if(length(initialPosition[[i]])!=3 | !is.numeric(initialPosition[[i]]) | any(!is.finite(initialPosition[[i]]))) stop("each element of initialPosition must be a finite numeric vector of length 3")
      }
    }
  } else {
    if(is.null(mvnCoords) || dist[[mvnCoords]] %in% c("mvnorm2","rw_mvnorm2")){
      if(length(initialPosition)!=2 | !is.numeric(initialPosition) | any(!is.finite(initialPosition))) stop("initialPosition must be a finite numeric vector of length 2")
    } else if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")){
      if(all(initialPosition==0)) initialPosition <- c(0,0,0)
      if(length(initialPosition)!=3 | !is.numeric(initialPosition) | any(!is.finite(initialPosition))) stop("initialPosition must be a finite numeric vector of length 3")
    }
    tmpPos<-initialPosition
    initialPosition<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      initialPosition[[i]]<-tmpPos
    }
  }
  
  if(is.null(nbHierCovs)){
    wnbHierCovs <- data.tree::Node$new(hierDist$Get("name",filterFun=isRoot))
    wnbHierCovs$AddChild(hierDist$Get("name",filterFun=function(x) x$level==2)[1],nbCovs=0)
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)[-1]){
      #wnbHierCovs$AddChild(paste0(j,"i"),nbCovs=0)
      wnbHierCovs$AddChild(j,nbCovs=0)
    }
  } else {
    if(!inherits(nbHierCovs,"Node")) stop("'nbHierCovs' must be of class Node; see ?data.tree::Node")
    if(!("nbCovs" %in% nbHierCovs$fieldsAll)) stop("'nbHierCovs' must include a 'nbCovs' field")
    if(!data.tree::AreNamesUnique(nbHierCovs)) stop("node names in 'nbHierCovs' must be unique")
    if(nbHierCovs$height!=2) stop("'nbHierCovs' hierarchy must contain 1 level")
    
    wnbHierCovs <- data.tree::Clone(nbHierCovs)
    
    #if(!all(hierDist$Get("name",filterFun=function(x) x$level==2)==wnbHierCovs$Get("name",filterFun=function(x) x$level==2)[seq(1,wnbHierCovs$count,2)])) stop("'hierDist' and 'nbHierCovs' are not consistent; check number of nodes, node names, and node order")
    for(j in hierDist$Get("name",filterFun=function(x) x$level==2)){
      if(is.null(wnbHierCovs[[j]])) {
        wnbHierCovs$AddChild(j,nbCovs=0)
      } else if(is.null(wnbHierCovs[[j]]$nbCovs)){
        wnbHierCovs[[j]]$nbCovs <- 0
      } else if(!is.numeric(wnbHierCovs[[j]]$nbCovs)) stop("'nbHierCovs$",j,"$nbCovs' must be numeric")
      
      #if(j!=hierDist$Get("name",filterFun=function(x) x$level==2)[1]){
      #  if(is.null(wnbHierCovs[[paste0(j,"i")]])) {
      #    wnbHierCovs$AddChild(paste0(j,"i"),nbCovs=0)
      #  } else if(is.null(wnbHierCovs[[paste0(j,"i")]]$nbCovs)){
      #    wnbHierCovs[[paste0(j,"i")]]$nbCovs <- 0
      #  } else if(!is.numeric(wnbHierCovs[[paste0(j,"i")]]$nbCovs)) stop("'nbHierCovs$",paste0(j,"i"),"$nbCovs' must be numeric")
      #}
    }
  }
  
  if(!all(sort(wnbHierCovs$Get("name",filterFun=function(x) x$level==2))==sort(hierDist$Get("name",filterFun=function(x) x$level==2)))) 
    stop("'nbHierCovs' level types can only include ",paste(hierDist$Get("name",filterFun=function(x) x$level==2),collapse=", "))
  
  nbCovs <- data.tree::Aggregate(wnbHierCovs,"nbCovs",sum,filterFun=isLeaf)
  
  if(is.null(model)){
    if(!is.null(covs) & nbCovs>0) {
      if(ncol(covs)!=nbCovs)
        warning("covs and nbCovs argument conflicting - nbCovs was set to ncol(covs)")
    }
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
  
  #######################################
  ## Prepare parameters for simulation ##
  #######################################
  # define number of observations for each animal
  allNbObs <- rep(NA,nbAnimals)
  levelObs <- list()
  level <- lLevels <- list()
  for(zoo in 1:nbAnimals) {
    levelObs[[zoo]] <- list()
    for(j in obsPerLevel[[zoo]]$Get("name",filterFun=function(x) x$level==2)){
      if(obsPerLevel[[zoo]][[j]]$obs[1]!=obsPerLevel[[zoo]][[j]]$obs[2])
        levelObs[[zoo]][[j]] <- sample(obsPerLevel[[zoo]][[j]]$obs[1]:obsPerLevel[[zoo]][[j]]$obs[2],size=1)
      else
        levelObs[[zoo]][[j]] <- obsPerLevel[[zoo]][[j]]$obs[1]
    }
    lInd <- unname(gsub("level","",sort(obsPerLevel[[zoo]]$Get("name",filterFun=function(x) x$level==2),decreasing=TRUE)[1]))
    level[[zoo]] <- c(paste0(lInd,"i"),rep(lInd,levelObs[[zoo]][[paste0("level",lInd)]]))
    lLevels[[zoo]] <- c(paste0(lInd,"i"),lInd)
    for(j in sort(obsPerLevel[[zoo]]$Get("name",filterFun=function(x) x$level==2),decreasing=TRUE)[-1]){
      level[[zoo]] <- c(gsub("level","",j),level[[zoo]])
      lLevels[[zoo]] <- c(gsub("level","",j),lLevels[[zoo]])
      level[[zoo]] <- rep(level[[zoo]],levelObs[[zoo]][[j]])
      if(j!="level1") {
        level[[zoo]] <- c(paste0(gsub("level","",j),"i"),level[[zoo]])
        lLevels[[zoo]] <- c(paste0(gsub("level","",j),"i"),lLevels[[zoo]])
      }
    }
    allNbObs[zoo] <- length(level[[zoo]]) #prod(unlist(levelObs[[zoo]]))
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
    rownames(covs) <- 1:sum(allNbObs)
  }
  
  ###############################
  ## Simulate covariate values ##
  ###############################
  allCovs <- NULL
  if(nbCovs>0) {
    if(is.null(covs)) {
      for(i in hierDist$Get("name",filterFun=function(x) x$level==2)){
        lnbCovs <- wnbHierCovs[[i]]$nbCovs
        if(lnbCovs>0){
          for(j in 1:lnbCovs){
            if(is.null(allCovs)){
              allCovs <- data.frame(rep(NA,sum(allNbObs)))
              lInd <- which(unlist(level)==gsub("level","",i))
              allCovs[lInd,] <- rnorm(length(lInd))
              colnames(allCovs) <- paste("cov",gsub("level","",i),".",j,sep="")
            } else {
              c <- data.frame(rep(NA,sum(allNbObs)))
              lInd <- which(unlist(level)==gsub("level","",i))
              c[lInd,] <- rnorm(length(lInd))
              colnames(c) <- paste("cov",gsub("level","",i),".",j,sep="")
              allCovs <- cbind(allCovs,c)
            }
          }
        }
      }
      # account for missing values of the covariates
      for(i in 1:nbCovs) {
        if(length(which(is.na(allCovs[,i])))>0) { # if covariate i has missing values
          if(is.na(allCovs[1,i])) { # if the first value of the covariate is missing
            k <- 1
            while(is.na(allCovs[k,i])) k <- k+1
            for(j in k:2) allCovs[j-1,i] <- allCovs[j,i]
          }
          for(j in 2:nrow(allCovs))
            if(is.na(allCovs[j,i])) allCovs[j,i] <- allCovs[j-1,i]
        }
      }
    } else {
      allCovs <- covs
    }
  }
  
  if(anyDuplicated(colnames(allCovs))) stop("covariates must have unique names")
  if(anyDuplicated(spatialcovnames)) stop("spatialCovs must have unique names")
  if(!is.null(model) & nbSpatialCovs>0){
    spInd <- which(!(colnames(allCovs) %in% spatialcovnames))
    if(length(spInd)) {
      allCovs <- allCovs[,spInd,drop=FALSE]
      nbCovs <- ncol(allCovs)
    } else {
      allCovs <- NULL
      nbCovs <- 0
    }
  } else if(any(colnames(allCovs) %in% spatialcovnames)) stop("spatialCovs name(s) cannot match other covariate name(s)")
  
  if(!all(angleCovs %in% c(colnames(allCovs),spatialcovnames))){
    stop("angleCovs ",paste0("'",angleCovs[!(angleCovs %in% c(colnames(allCovs),spatialcovnames))],"'",collapse=", ")," not found in covs or spatialCovs")
  }
  
  centerInd<-NULL
  if(!is.null(centers)){
    if(!is.matrix(centers)) stop("centers must be a matrix")
    if(dim(centers)[2]!=2) stop("centers must be a matrix consisting of 2 columns (i.e., x- and y-coordinates)")
    centerInd <- which(!apply(centers,1,function(x) any(is.na(x))))
    if(length(centerInd)){
      if(is.null(rownames(centers))) centerNames<-paste0("center",rep(centerInd,each=2),".",rep(c("dist","angle"),length(centerInd)))
      else centerNames <- paste0(rep(rownames(centers),each=2),".",rep(c("dist","angle"),length(centerInd)))
      centerCovs <- data.frame(matrix(NA,nrow=sum(allNbObs),ncol=length(centerInd)*2,dimnames=list(NULL,centerNames)))
    }  
  } else centerNames <- NULL
  
  centroidInd<-NULL
  if(!is.null(centroids)){
    if(!is.list(centroids)) stop("centroids must be a named list")
    centroidNames <- character()
    for(j in 1:length(centroids)){
      if(!is.data.frame(centroids[[j]])) stop("each element of centroids must be a data frame")
      if(dim(centroids[[j]])[1]<max(allNbObs) | dim(centroids[[j]])[2]!=2) stop("each element of centroids must be a data frame consisting of at least ",max(allNbObs)," rows (i.e., the maximum number of observations per animal) and 2 columns (i.e., x- and y-coordinates)")
      if(!all(c("x","y") %in% colnames(centroids[[j]]))) stop("centroids columns must be named 'x' (x-coordinate) and 'y' (y-coordinate)")
      #centroidInd <- which(!apply(centroids[[j]],1,function(x) any(is.na(x))))
      #if(length(centroidInd)){
      if(any(is.na(centroids[[j]]))) stop("centroids cannot contain missing values")
      if(is.null(names(centroids[j]))) centroidNames <- c(centroidNames,paste0("centroid",rep(j,each=2),".",c("dist","angle")))
      else centroidNames <- c(centroidNames,paste0(rep(names(centroids[j]),each=2),".",c("dist","angle")))
    }
    centroidCovs <- data.frame(matrix(NA,nrow=sum(allNbObs),ncol=length(centroidNames),dimnames=list(NULL,centroidNames)))
    centroidInd <- length(centroidNames)/2
    #}  
  } else centroidNames <- NULL
  
  if(!is.null(model) & length(centerInd)){
    cInd <- which(!(colnames(allCovs) %in% centerNames))
    if(length(cInd)) {
      allCovs <- allCovs[,cInd,drop=FALSE]
      nbCovs <- ncol(allCovs)
    } else {
      allCovs <- NULL
      nbCovs <- 0
    }
  } else if(any(colnames(allCovs) %in% centerNames)) stop("centers name(s) cannot match other covariate name(s)")
  
  if(!is.null(model) & length(centroidInd)){
    cInd <- which(!(colnames(allCovs) %in% centroidNames))
    if(length(cInd)) {
      allCovs <- allCovs[,cInd,drop=FALSE]
      nbCovs <- ncol(allCovs)
    } else {
      allCovs <- NULL
      nbCovs <- 0
    }
  } else if(any(colnames(allCovs) %in% centroidNames)) stop("centroids name(s) cannot match other covariate name(s)")
  
  allNbCovs <- nbCovs+nbSpatialCovs
  
  zeroMass<-oneMass<-vector('list',length(inputs$dist))
  names(zeroMass)<-names(oneMass)<-distnames
  
  allStates <- NULL
  allSpatialcovs<-NULL
  
  #make sure 'step' preceeds 'angle'
  if(all(c("step","angle") %in% distnames)){
    distnames<-c(distnames[!(distnames %in% c("step","angle"))],"step","angle")
  }
  
  # build the data frame to be returned
  data<-data.frame(ID=factor())
  for(i in distnames){
    if(dist[[i]] %in% mvndists){
      data[[paste0(i,".x")]] <- numeric()
      data[[paste0(i,".y")]] <- numeric()
      if(dist[[i]] %in% c("mvnorm3","rw_mvnorm3")) data[[paste0(i,".z")]] <- numeric()
    } else data[[i]]<-numeric()
  }
  if("angle" %in% distnames){ 
    if(inputs$dist[["angle"]] %in% angledists & ("step" %in% distnames))
      if(inputs$dist[["step"]] %in% stepdists){
        if(hierDist$Get("parent",filterFun=isLeaf)$step$name != hierDist$Get("parent",filterFun=isLeaf)$angle$name) stop("step and angle must be in the same level of the hierarchy")
        data$x<-numeric()
        data$y<-numeric()
        coordLevel <- gsub("level","",hierDist$Get("parent",filterFun=isLeaf)$step$name)
      }
  } else if("step" %in% distnames){
    if(inputs$dist[["step"]] %in% stepdists){
      data$x<-numeric()
      data$y<-numeric()
      coordLevel <- gsub("level","",hierDist$Get("parent",filterFun=isLeaf)$step$name)
    }    
  } else if(!is.null(mvnCoords)){
    data[[paste0(mvnCoords,".x")]]<-numeric()
    data[[paste0(mvnCoords,".y")]]<-numeric()
    if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) data[[paste0(mvnCoords,".z")]]<-numeric()
    coordLevel <- gsub("level","",hierDist$Get("parent",filterFun=isLeaf)[[mvnCoords]]$name)
  } else {
    if(nbSpatialCovs | length(centerInd) | length(centroidInd) | length(angleCovs)) stop("spatialCovs, angleCovs, centers, and/or centroids cannot be specified without valid step length and turning angle distributions")
    coordLevel <- NULL
  }
  
  if(!is.null(coordLevel)) distCoordLevel <- names(hierDist[[paste0("level",coordLevel)]]$children)
  
  rwInd <- any(unlist(dist) %in% rwdists)
  
  #if(is.null(formula)) {
  #  if(allNbCovs) formula <- formula(paste0("~",paste0(c(colnames(allCovs),spatialcovnames),collapse="+")))
  #  else formula <- formula(~1)
  #}
  
  if(length(all.vars(formula)))
    if(!all(all.vars(formula) %in% c("ID","level",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formula' covariate(s) not found")
  if(length(all.vars(formPi)))
    if(!all(all.vars(formPi) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formulaPi' covariate(s) not found")
  if(length(all.vars(formDelta)))
    if(!all(all.vars(formDelta) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formulaDelta' covariate(s) not found")
  if(("ID" %in% all.vars(formula) | "ID" %in% all.vars(formPi) | "ID" %in% all.vars(formDelta)) & nbAnimals<2) stop("ID cannot be a covariate when nbAnimals=1")
  
  tmpCovs <- data.frame(ID=factor(1,levels=1:nbAnimals))
  if(!is.null(allCovs))
    tmpCovs <- cbind(tmpCovs,allCovs[1,,drop=FALSE])
  if(nbSpatialCovs){
    for(j in 1:nbSpatialCovs){
      for(i in 1:length(initialPosition)){
        if(is.na(raster::cellFromXY(spatialCovs[[j]],initialPosition[[i]]))) stop("initialPosition for individual ",i," is not within the spatial extent of the ",spatialcovnames[j]," raster")
      }
      getCell<-raster::cellFromXY(spatialCovs[[j]],initialPosition[[1]])
      spCov <- spatialCovs[[j]][getCell]
      if(inherits(spatialCovs[[j]],c("RasterStack","RasterBrick"))){
        zname <- names(attributes(spatialCovs[[j]])$z)
        zvalues <- raster::getZ(spatialCovs[[j]])
        spCov <- spCov[1,which(zvalues==tmpCovs[[zname]][1])]
      }
      tmpCovs[[spatialcovnames[j]]]<-spCov
    }
  }
  if(length(centerInd)){
    for(j in 1:length(centerInd)){
      tmpDistAngle <- distAngle(initialPosition[[1]],initialPosition[[1]],centers[centerInd[j],])
      tmpCovs[[centerNames[(j-1)*2+1]]]<- tmpDistAngle[1]
      tmpCovs[[centerNames[(j-1)*2+2]]]<- tmpDistAngle[2]
    }
  }
  if(length(centroidInd)){
    for(j in 1:centroidInd){
      tmpDistAngle <- distAngle(initialPosition[[1]],initialPosition[[1]],as.numeric(centroids[[j]][1,]))
      tmpCovs[[centroidNames[(j-1)*2+1]]]<- tmpDistAngle[1]
      tmpCovs[[centroidNames[(j-1)*2+2]]]<- tmpDistAngle[2]
    }
  }
  
  tmpCovs$level <- factor(level[[1]][1],levels=lLevels[[1]])
  
  tmpCovs <- tmpCovs[rep(1,length(lLevels[[1]])),,drop=FALSE]
  tmpCovs$level <- factor(lLevels[[1]],levels=lLevels[[1]])
  class(tmpCovs) <- append("hierarchical",class(tmpCovs))
  
  inputHierHMM <- formatHierHMM(data=tmpCovs,hierStates,hierDist,hierBeta,hierDelta,hierFormula,hierFormulaDelta,mixtures,workBounds,checkData=FALSE)
  nbStates <- inputHierHMM$nbStates
  dist <- inputHierHMM$dist
  formula <- inputHierHMM$formula
  formulaDelta <- inputHierHMM$formulaDelta
  betaRef <- inputHierHMM$betaRef
  stateNames <- inputHierHMM$stateNames
  
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  # build design matrix for recharge model
  if(!is.null(recharge)){
    reForm <- formatRecharge(nbStates,formula,data=tmpCovs)
    formulaStates <- reForm$formulaStates
    formterms <- newForm$formterms
    tmpCovs[colnames(reForm$newdata)] <- reForm$newdata
    recharge <- reForm$recharge
    hierRecharge <- reForm$hierRecharge
    recLevelNames <- names(hierRecharge)
    rechargeNames <- paste0("recharge",gsub("level","",recLevelNames))
    g0covs <- reForm$g0covs
    nbG0covs <- ncol(g0covs)-1
    recovs <- reForm$recovs
    nbRecovs <- ncol(recovs)-1
  } else {
    nbRecovs <- 0
    nbG0covs <- 0
    g0covs <- NULL
    recovs <- NULL
  }
  
  if(is.null(model)){
    beta <- inputHierHMM$beta
    delta <- inputHierHMM$delta
  }
  
  wworkBounds <- inputHierHMM$workBounds
  
  if(is.null(model) & !is.list(beta)){
      beta0 <- list(beta=beta)
  } else beta0 <- beta
  
  nbBetaCovs <- ncol(model.matrix(newformula,tmpCovs))
  
  if(is.null(beta0$beta)){
    beta0$beta <- matrix(rnorm(nbStates*(nbStates-1)*nbBetaCovs*mixtures)[inputHierHMM$betaCons],nrow=nbBetaCovs*mixtures)
    beta0$beta[which(!is.na(inputHierHMM$fixPar$beta))] <- inputHierHMM$fixPar$beta[which(!is.na(inputHierHMM$fixPar$beta))]
  }
  if(nbRecovs){
    if(is.null(beta0$g0) | is.null(beta0$theta)) stop("beta$g0 and beta$theta must be specified for recharge model")
    if(length(beta0$g0)!=(nbG0covs+1) | any(!is.numeric(beta0$g0))) stop("beta$g0 must be a numeric vector of length ",nbG0covs+1)
    if(length(beta0$theta)!=(nbRecovs+1) | any(!is.numeric(beta0$theta))) stop("beta$theta must be a numeric vector of length ",nbRecovs+1)
  }
  if(ncol(beta0$beta)!=nbStates*(nbStates-1) | (nrow(beta0$beta)/mixtures)!=nbBetaCovs) {
    error <- paste("beta has wrong dimensions: it should have",nbBetaCovs*mixtures,"rows and",
                   nbStates*(nbStates-1),"columns.")
    stop(error)
  }
  if(nbStates>1){
    for(state in 1:(nbStates*(nbStates-1))){
      noBeta<-which(match(colnames(model.matrix(newformula,tmpCovs)),colnames(model.matrix(formulaStates[[state]],tmpCovs)),nomatch=0)==0)
      if(length(noBeta)) beta0$beta[noBeta,state] <- 0
    }
  }
  
  covsPi <- model.matrix(formPi,tmpCovs)
  nbCovsPi <- ncol(covsPi)-1
  if(!nbCovsPi & is.null(formulaPi)){
    if(is.null(beta0$pi)){
      beta0$pi <- matrix(1/mixtures,(nbCovsPi+1),mixtures)
    } else {
      beta0$pi <- matrix(beta0$pi,(nbCovsPi+1),mixtures)
    }
    if(length(beta0$pi) != (nbCovsPi+1)*mixtures)
      stop(paste("beta$pi has the wrong length: it should have",mixtures,"elements."))
    beta0$pi <- matrix(log(beta0$pi[-1]/beta0$pi[1]),nbCovsPi+1,mixtures-1)
  } else {
    if(is.null(beta0$pi)) beta0$pi <- matrix(0,nrow=(nbCovsPi+1),ncol=mixtures-1)
    if(is.null(dim(beta0$pi)) || (ncol(beta0$pi)!=mixtures-1 | nrow(beta0$pi)!=(nbCovsPi+1)))
      stop(paste("beta$pi has wrong dimensions: it should have",(nbCovsPi+1),"rows and",
                 mixtures-1,"columns."))
  }
  
  covsDelta <- model.matrix(formDelta,tmpCovs)
  nbCovsDelta <- ncol(covsDelta)-1
  # initial state distribution
  if(is.null(delta)) {
    delta0 <- matrix(0,(nbCovsDelta+1)*mixtures,nbStates-1)
    delta0[which(!is.na(inputHierHMM$fixPar$delta))] <- inputHierHMM$fixPar$delta[which(!is.na(inputHierHMM$fixPar$delta))]
  } else delta0 <- delta
  #if(!nbCovsDelta & is.null(formulaDelta)){
  #  if(mixtures==1){
  #    if(length(delta0) != (nbCovsDelta+1)*nbStates)
  #      stop(paste("delta has the wrong length: it should have",nbStates*mixtures,"elements."))
  #    deltaB <- matrix(log(delta0[-1]/delta0[1]),1)
  #  } else {
  #    if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates | nrow(delta0)!=mixtures))
  #      stop(paste("delta has wrong dimensions: it should have",mixtures,"rows and",
  #                 nbStates,"columns."))
  #    deltaB <- apply(delta0,1,function(x) log(x[-1]/x[1]))
  #  }
    
  #} else {
  if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates-1 | nrow(delta0)!=(nbCovsDelta+1)*mixtures))
    stop(paste("delta has wrong dimensions: it should have",(nbCovsDelta+1)*mixtures,"rows and",
               nbStates-1,"columns."))
  deltaB <- delta0
  #}
  
  parCount<- lapply(Par[distnames],length)
  parindex <- c(0,cumsum(unlist(parCount))[-length(distnames)])
  names(parindex) <- distnames
  
  wworkBounds <- getWorkBounds(wworkBounds,distnames,unlist(Par[distnames]),parindex,parCount,inputs$DM,beta0,deltaB)
  
  wnbeta <- w2wn(beta0$beta,wworkBounds$beta)
  wnpi <- w2wn(beta0$pi,wworkBounds$pi)
  if(!is.null(recharge)){
    wng0 <- w2wn(beta0$g0,wworkBounds$g0)
    wntheta <- w2wn(beta0$theta,wworkBounds$theta)
  }
  
  mix <- rep(1,nbAnimals)
  
  
  printMessage(nbStates,dist,p,DM,formula,formDelta,formPi,mixtures,"Simulating",hierarchical=TRUE)
  
  if(!nbSpatialCovs | !retrySims){
    
    ###########################
    ## Loop over the animals ##
    ###########################
    for (zoo in 1:nbAnimals) {
      
      # number of observations for animal zoo
      nbObs <- allNbObs[zoo]
      d <- data.frame(ID=factor(rep(zoo,nbObs)),level=factor(level[[zoo]],levels=lLevels[[zoo]]))
      
      ###############################
      ## Simulate covariate values ##
      ###############################
      subCovs<-data.frame(ID=rep(factor(zoo,levels=1:nbAnimals),nbObs),level=factor(level[[zoo]],levels=lLevels[[zoo]]))
      if(nbCovs>0) {
        # select covariate values which concern the current animal
        if(zoo<2)
          ind1 <- 1
        else
          ind1 <- sum(allNbObs[1:(zoo-1)])+1
        ind2 <- sum(allNbObs[1:zoo])
        subCovs <- cbind(subCovs,data.frame(allCovs[ind1:ind2,,drop=FALSE]))
      }
      if(length(centerInd)) subCovs <- cbind(subCovs,centerCovs[cumNbObs[zoo]+1:nbObs,])
      if(length(centroidInd)) subCovs <- cbind(subCovs,centroidCovs[cumNbObs[zoo]+1:nbObs,])
      
      subSpatialcovs<-as.data.frame(matrix(NA,nrow=nbObs,ncol=nbSpatialCovs))
      colnames(subSpatialcovs)<-spatialcovnames
      subAnglecovs <- as.data.frame(matrix(NA,nrow=nbObs,ncol=length(angleCovs)))
      colnames(subAnglecovs) <- angleCovs
      
      X <- matrix(NA,nrow=nbObs,ncol=2)
      if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) X <- matrix(0,nrow=nbObs,ncol=3)
      X[1,] <- initialPosition[[zoo]] # initial position of animal
      
      phi <- 0
      
      ############################
      ## Simulate movement path ##
      ############################
      
      genData <- genArgs <- vector('list',length(distnames))
      names(genData) <- names(genArgs) <- distnames
      
      for(i in distnames){
        if(inputs$dist[[i]] %in% mvndists){
          genData[[i]] <- matrix(NA,nbObs,ncol(X))
        } else genData[[i]] <- rep(NA,nbObs)
        genArgs[[i]] <- list(1)  # first argument = 1 (one random draw)
      }
      
      gamma <- matrix(0,nbStates,nbStates)
      gamma[cbind(1:nbStates,betaRef)] <- 1
      gamma <- t(gamma)
      
      class(subCovs) <- append("hierarchical",class(subCovs))
      
      if(!nbSpatialCovs & !length(centerInd) & !length(centroidInd) & !length(angleCovs) & !rwInd) {
        if(!is.null(recharge)){
          for(j in rechargeNames){
            subCovs[[j]] <- formatRecharge(nbStates,formula,data=subCovs,par=list(g0=wng0,theta=wntheta))$newdata[[j]]
          }
        }
        DMcov <- model.matrix(newformula,subCovs)
        
        # format parameters
        DMinputs<-getDM(subCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind,inputs$dist)
        if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")
        
        ncmean <- get_ncmean(distnames,fullDM,inputs$circularAngleMean,nbStates)
        nc <- ncmean$nc
        meanind <- ncmean$meanind
        
        covsDelta <- model.matrix(formDelta,subCovs[1,,drop=FALSE])
        covsPi <- model.matrix(formPi,subCovs[1,,drop=FALSE])
        fullsubPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nbObs,inputs$dist,p$Bndind,nc,meanind,covsDelta,wworkBounds,covsPi)
        
        pie <- fullsubPar$pi
        
        # assign individual to mixture
        if(mixtures>1) mix[zoo] <- sample.int(mixtures,1,prob=pie)
        
        gFull <-  DMcov %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
        g <- gFull[1,,drop=FALSE]
        delta0 <- fullsubPar$delta[mix[zoo],]
      } else {
        
        if(nbSpatialCovs){
          for(j in 1:nbSpatialCovs){
            getCell<-raster::cellFromXY(spatialCovs[[j]],c(X[1,1],X[1,2]))
            if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
            spCov <- spatialCovs[[j]][getCell]
            if(inherits(spatialCovs[[j]],c("RasterStack","RasterBrick"))){
              zname <- names(attributes(spatialCovs[[j]])$z)
              zvalues <- raster::getZ(spatialCovs[[j]])
              spCov <- spCov[1,which(zvalues==subCovs[1,zname])]
            }
            subSpatialcovs[1,j]<-spCov
            if(spatialcovnames[j] %in% angleCovs) {
              subAnglecovs[1,spatialcovnames[j]] <- subSpatialcovs[1,j]
              subSpatialcovs[1,j] <- 0  # set to zero because can't have NA covariates
            }
          }
        }
        
        for(j in angleCovs[which(angleCovs %in% names(subCovs))]){
          subAnglecovs[1,j] <- subCovs[1,j]
          subCovs[1,j] <- 0 # set to zero because can't have NA covariates
        }
        
        if(length(centerInd)){
          for(j in 1:length(centerInd)){
            subCovs[1,centerNames[(j-1)*2+1:2]]<-distAngle(X[1,],X[1,],centers[centerInd[j],])
          }
        }
        
        if(length(centroidInd)){
          for(j in 1:centroidInd){
            subCovs[1,centroidNames[(j-1)*2+1:2]]<-distAngle(X[1,],X[1,],as.numeric(centroids[[j]][1,]))
          }
        }
        if(!is.null(recharge)){
          tmpSubSpatCovs <- subCovs[1,,drop=FALSE]
          tmpSubSpatCovs[colnames(subSpatialcovs[1,,drop=FALSE])] <- subSpatialcovs[1,,drop=FALSE]
          tmpSubSpatCovs <- tmpSubSpatCovs[rep(1,length(lLevels[[1]])),,drop=FALSE]
          tmpSubSpatCovs$level <- factor(lLevels[[1]],levels=lLevels[[1]])
          for(j in rechargeNames){
            subCovs[1,j] <- formatRecharge(nbStates,formula,data=tmpSubSpatCovs,par=list(g0=wng0,theta=wntheta))$newdata[[j]][1]
          }
          #g0 <- model.matrix(recharge$g0,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wng0
          #subCovs[1,"recharge"] <- g0 #+ model.matrix(recharge$theta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wntheta
        }
        
        for(i in distnames){
          if(dist[[i]] %in% rwdists){
            if(dist[[i]] %in% c("rw_mvnorm2")){
              subCovs[1,paste0(i,".x_tm1")] <- X[1,1]
              subCovs[1,paste0(i,".y_tm1")] <- X[1,2]
            } else if(dist[[i]] %in% c("rw_mvnorm3")){
              subCovs[1,paste0(i,".x_tm1")] <- X[1,1]
              subCovs[1,paste0(i,".y_tm1")] <- X[1,2]
              subCovs[1,paste0(i,".z_tm1")] <- X[1,3]
            }
          }
        }
        
        # get max crw lag
        maxlag <- getDM(cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]),inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean,wlag=TRUE)$lag
        
        covsPi <- model.matrix(formPi,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        pie <- mlogit(wnpi,covsPi,nbCovsPi,1,mixtures)
        
        # assign individual to mixture
        if(mixtures>1) mix[zoo] <- sample.int(mixtures,1,prob=pie)
        
        g <- model.matrix(newformula,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
        covsDelta <- model.matrix(formDelta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        
        delta0 <- mlogit(deltaB[(mix[zoo]-1)*(nbCovsDelta+1)+1:(nbCovsDelta+1),,drop=FALSE],covsDelta,nbCovsDelta,1,nbStates)
        #delta0 <- c(rep(0,nbCovsDelta+1),deltaB[(mix[zoo]-1)*(nbCovsDelta+1)+1:(nbCovsDelta+1),])
        #deltaXB <- covsDelta %*% matrix(delta0,nrow=nbCovsDelta+1)
        #expdelta <- exp(deltaXB)
        #delta0 <- expdelta/rowSums(expdelta)
        #for(i in which(!is.finite(rowSums(delta0)))){
        #  tmp <- exp(Brobdingnag::as.brob(deltaXB[i,]))
        #  delta0[i,] <- as.numeric(tmp/Brobdingnag::sum(tmp))
        #}
      }
      gamma[!gamma] <- exp(g)
      gamma <- t(gamma)
      gamma <- gamma/apply(gamma,1,sum)
      
      if(nbStates>1) {
        Z <- rep(NA,nbObs)
        Z[1] <- sample(1:nbStates,size=1,prob=delta0%*%gamma)
      } else
        Z <- rep(1,nbObs)
      
      #for(kk in sort(hierDist$Get("name",filterFun=function(x) x$level==2))){
        
        #kIndi <- which(subCovs$level==paste0(gsub("level","",kk),"i"))
        
        #if(length(kIndi)){
        #  if(nbSpatialCovs |  length(centerInd) | length(centroidInd) | length(angleCovs) | rwInd){
        #    subCovs[kIndi,which(colnames(subCovs)!="level")] <- subCovs[kIndi-1,which(colnames(subCovs)!="level"),drop=FALSE]
        #    subSpatialcovs[kIndi,] <- subSpatialcovs[kIndi-1,,drop=FALSE]
        #  }
        #  X[kIndi,] <- X[kIndi-1,]
        #}
        
        #kInd <- which(subCovs$level==gsub("level","",kk))
      
        if(!is.null(coordLevel)) coordNA <- which(level[[zoo]]==coordLevel)[sum(level[[zoo]]==coordLevel)]
      
        for (k in 1:nbObs){
          
          #if(level[[zoo]][k] %in% lLevels[[zoo]][seq(2,length(lLevels[[zoo]]),2)]){
          #  if(nbSpatialCovs | length(centerInd) | length(centroidInd) | length(angleCovs) | rwInd){
          #    subCovs[k,which(colnames(subCovs)!="level")] <- subCovs[k-1,which(colnames(subCovs)!="level"),drop=FALSE]
          #    subSpatialcovs[k,] <- subSpatialcovs[k-1,,drop=FALSE]
          #    g <- model.matrix(newformula,cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
          #  } else {
          #    g <- gFull[k,,drop=FALSE]
          #  }
            
          #  X[k,] <- X[k-1,]
            
          #  # get initial state
          #  gamma <- matrix(0,nbStates,nbStates)
          #  gamma[cbind(1:nbStates,betaRef)] <- 1
          #  gamma <- t(gamma)
          #  gamma[!gamma] <- exp(g)
          #  gamma <- t(gamma)
          #  gamma <- gamma/apply(gamma,1,sum)
          #  Z[k] <- sample(1:nbStates,size=1,prob=gamma[Z[k-1],])  
          #}
            
          if(!is.null(names(as.list(hierDist)[[paste0("level",level[[zoo]][k])]]))){
            
            if(nbSpatialCovs |  length(centerInd) | length(centroidInd) | length(angleCovs) | rwInd){
              # format parameters
              DMinputs<-getDM(cbind(subCovs[k-maxlag:0,,drop=FALSE],subSpatialcovs[k-maxlag:0,,drop=FALSE]),inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean,wlag=TRUE)
              fullDM <- DMinputs$fullDM
              DMind <- DMinputs$DMind
              wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind,inputs$dist)
              if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")
              
              nc <- meanind <- vector('list',length(distnames))
              names(nc) <- names(meanind) <- distnames
              for(i in distnames){
                nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
                if(!isFALSE(inputs$circularAngleMean[[i]])) {
                  meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
                  # deal with angular covariates that are exactly zero
                  if(length(meanind[[i]])){
                    angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
                    sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
                    nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
                    nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
                  }
                }
              }
              subPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,inputs$dist,p$Bndind,nc,meanind,covsDelta,wworkBounds,covsPi)
            } else {
              subPar <- lapply(fullsubPar[distnames],function(x) x[,k,drop=FALSE])#fullsubPar[,k,drop=FALSE]
            }
            
            for(i in distnames[which(distnames %in% names(as.list(hierDist)[[paste0("level",level[[zoo]][k])]]))]){
              
              zeroMass[[i]] <- rep(0,nbStates)
              oneMass[[i]] <- rep(0,nbStates)
              if(zeroInflation[[i]] | oneInflation[[i]]) {
                if(zeroInflation[[i]]) zeroMass[[i]] <- subPar[[i]][parSize[[i]]*nbStates-nbStates*oneInflation[[i]]-(nbStates-1):0]
                if(oneInflation[[i]])  oneMass[[i]] <- subPar[[i]][parSize[[i]]*nbStates-(nbStates-1):0]
                subPar[[i]] <- subPar[[i]][-(parSize[[i]]*nbStates-(nbStates*oneInflation[[i]]-nbStates*zeroInflation[[i]]-1):0)]
              }
              
              if(inputs$dist[[i]] %in% mvndists){
                if(inputs$dist[[i]]=="mvnorm2" || inputs$dist[[i]]=="rw_mvnorm2"){
                  genArgs[[i]][[2]] <- c(subPar[[i]][Z[k]],
                                         subPar[[i]][nbStates+Z[k]])
                  genArgs[[i]][[3]] <- matrix(c(subPar[[i]][nbStates*2+Z[k]], #x
                                                subPar[[i]][nbStates*3+Z[k]], #xy
                                                subPar[[i]][nbStates*3+Z[k]], #xy
                                                subPar[[i]][nbStates*4+Z[k]]) #y
                                              ,2,2)
                } else if(inputs$dist[[i]]=="mvnorm3" || inputs$dist[[i]]=="rw_mvnorm3"){
                  genArgs[[i]][[2]] <- c(subPar[[i]][Z[k]],
                                         subPar[[i]][nbStates+Z[k]],
                                         subPar[[i]][2*nbStates+Z[k]])
                  genArgs[[i]][[3]] <- matrix(c(subPar[[i]][nbStates*3+Z[k]], #x
                                                subPar[[i]][nbStates*4+Z[k]], #xy
                                                subPar[[i]][nbStates*5+Z[k]], #xz
                                                subPar[[i]][nbStates*4+Z[k]], #xy
                                                subPar[[i]][nbStates*6+Z[k]], #y
                                                subPar[[i]][nbStates*7+Z[k]], #yz
                                                subPar[[i]][nbStates*5+Z[k]], #xz
                                                subPar[[i]][nbStates*7+Z[k]], #yz
                                                subPar[[i]][nbStates*8+Z[k]]) #z          
                                              ,3,3)
                } 
              } else if(inputs$dist[[i]]=="cat"){
                genArgs[[i]][[2]] <- subPar[[i]][seq(Z[k],(parSize[[i]]+1)*nbStates,nbStates)]
              } else {
                for(j in 1:(parSize[[i]]-zeroInflation[[i]]-oneInflation[[i]]))
                  genArgs[[i]][[j+1]] <- subPar[[i]][(j-1)*nbStates+Z[k]]
              }
              
              if(inputs$dist[[i]] %in% angledists){
                
                genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
                if(genData[[i]][k] >  pi) genData[[i]][k] <- genData[[i]][k]-2*pi
                if(genData[[i]][k] < -pi) genData[[i]][k] <- genData[[i]][k]+2*pi
                
                if(i=="angle" & ("step" %in% distnames)){
                  if(inputs$dist[["step"]] %in% stepdists) {
                    if(k != coordNA){
                      if(genData$step[k]>0){
                        phi <- phi + genData[[i]][k]
                      } #else if(genData$step[k]==0) {
                      #genData[[i]][k] <- NA # angle = NA if step = 0
                      #}
                      m <- genData$step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
                      X[k+1,] <- X[k,] + m
                    } else {
                      for(ii in distCoordLevel){
                        d[[ii]][k] <- genData[[ii]][k] <- NA
                      }
                      #X[k+1,] <- X[k,]
                    }
                  }
                }
              } else {
                
                if(inputs$dist[[i]]=="gamma") {
                  shape <- genArgs[[i]][[2]]^2/genArgs[[i]][[3]]^2
                  scale <- genArgs[[i]][[3]]^2/genArgs[[i]][[2]]
                  genArgs[[i]][[2]] <- shape
                  genArgs[[i]][[3]] <- 1/scale # rgamma expects rate=1/scale
                }
                
                probs <- c(1.-zeroMass[[i]][Z[k]]-oneMass[[i]][Z[k]],zeroMass[[i]][Z[k]],oneMass[[i]][Z[k]])
                rU <- which(rmultinom(1,1,prob=probs)==1)
                if(rU==1){
                  if(inputs$dist[[i]] %in% mvndists){
                    genData[[i]][k,] <- do.call(Fun[[i]],genArgs[[i]])
                  } else genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
                } else if(rU==2) {
                  genData[[i]][k] <- 0
                } else {
                  genData[[i]][k] <- 1
                }
              }
              
              if(!is.null(mvnCoords) && i==mvnCoords){
                if(k < nbObs) X[k+1,] <- genData[[i]][k,]
                d[[i]] <- X
              } else d[[i]] <- genData[[i]]
              
            }
          }
          
          if(k < nbObs){
            
            if(all(is.na(X[k+1,]))){
              X[k+1,] <- X[k,]
            }

            if(nbSpatialCovs | length(centerInd) | length(centroidInd) | length(angleCovs) | rwInd){
              if(nbSpatialCovs){
                for(j in 1:nbSpatialCovs){
                  getCell<-raster::cellFromXY(spatialCovs[[j]],c(X[k+1,1],X[k+1,2]))
                  if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
                  spCov <- spatialCovs[[j]][getCell]
                  if(inherits(spatialCovs[[j]],c("RasterStack","RasterBrick"))){
                    zname <- names(attributes(spatialCovs[[j]])$z)
                    zvalues <- raster::getZ(spatialCovs[[j]])
                    spCov <- spCov[1,which(zvalues==subCovs[k+1,zname])]
                  }
                  subSpatialcovs[k+1,j]<-spCov
                  if(spatialcovnames[j] %in% angleCovs) {
                    subAnglecovs[k+1,spatialcovnames[j]] <- subSpatialcovs[k+1,j]
                    subSpatialcovs[k+1,j] <- circAngles(subAnglecovs[c(k,k+1),spatialcovnames[j]],data.frame(x=X[c(k,k+1),1],y=X[c(k,k+1),2]))[2] 
                  }
                }
              }
              
              for(j in angleCovs[which(angleCovs %in% names(subCovs))]){
                subAnglecovs[k+1,j] <- subCovs[k+1,j]
                subCovs[k+1,j] <- circAngles(subAnglecovs[c(k,k+1),j],data.frame(x=X[c(k,k+1),1],y=X[c(k,k+1),2]))[2] 
              }
              
              if(length(centerInd)){
                for(j in 1:length(centerInd)){
                  subCovs[k+1,centerNames[(j-1)*2+1:2]]<-distAngle(X[k,],X[k+1,],centers[centerInd[j],])
                }
              }
              if(length(centroidInd)){
                for(j in 1:centroidInd){
                  subCovs[k+1,centroidNames[(j-1)*2+1:2]]<-distAngle(X[k,],X[k+1,],as.numeric(centroids[[j]][k+1,]))
                }
              }
              if(!is.null(recharge)){
                colInd <- lapply(recLevelNames,function(x) which(grepl(paste0("I((level == \"",gsub("level","",x),"\")"),colnames(recovs),fixed=TRUE)))
                
                for(j in 1:length(rechargeNames)){
                  tmpSubCovs <- cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE])
                  if(subCovs[k+1,"level"]==gsub("level","",recLevelNames[j]) & match(subCovs[k,"level"],levels(subCovs$level))>=match(subCovs[k+1,"level"],levels(subCovs$level))){
                    tmpSubCovs$level <- subCovs[k+1,"level"]
                    subCovs[k+1,rechargeNames[j]] <- subCovs[k,rechargeNames[j]] + model.matrix(recharge$theta,tmpSubCovs)[,colInd[[j]],drop=FALSE] %*% wntheta[colInd[[j]]]
                  } else if(match(subCovs[k,"level"],levels(subCovs$level))>match(subCovs[k+1,"level"],levels(subCovs$level))){
                    #tmpSubCovs$level <- subCovs[k,"level"]
                    subCovs[k+1,rechargeNames[j]] <- subCovs[k,rechargeNames[j]] + model.matrix(recharge$theta,tmpSubCovs)[,colInd[[j]],drop=FALSE] %*% wntheta[colInd[[j]]]
                  } else {
                    subCovs[k+1,rechargeNames[j]] <- subCovs[k,rechargeNames[j]]
                  }
                }
                #subCovs[k+1,"recharge"] <- subCovs[k,"recharge"] + model.matrix(recharge$theta,cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE])) %*% wntheta
              }
              for(i in distnames){
                if(dist[[i]] %in% rwdists){
                  if(dist[[i]] %in% c("rw_mvnorm2")) subCovs[k+1,paste0(i,c(".x_tm1",".y_tm1"))] <- X[k+1,]
                  else if(dist[[i]] %in% c("rw_mvnorm3")) subCovs[k+1,paste0(i,c(".x_tm1",".y_tm1",".z_tm1"))] <- X[k+1,]
                }
              }
              g <- model.matrix(newformula,cbind(subCovs[k+1,,drop=FALSE],subSpatialcovs[k+1,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
            } else {
              g <- gFull[k+1,,drop=FALSE]
            }
            # get next state
            gamma <- matrix(0,nbStates,nbStates)
            gamma[cbind(1:nbStates,betaRef)] <- 1
            gamma <- t(gamma)
            gamma[!gamma] <- exp(g)
            gamma <- t(gamma)
            gamma <- gamma/apply(gamma,1,sum)
            Z[k+1] <- sample(1:nbStates,size=1,prob=gamma[Z[k],])  
          }
        }
      #}
      allStates <- c(allStates,Z)
      if(nbSpatialCovs>0) {
        allSpatialcovs <- rbind(allSpatialcovs,subSpatialcovs)
      }
      
      if("angle" %in% distnames){ 
        if(inputs$dist[["angle"]] %in% angledists & ("step" %in% distnames))
          if(inputs$dist[["step"]] %in% stepdists){
            d$angle[1] <- NA # the first angle value is arbitrary
            step0 <- which(d$step==0)
            d$angle[c(step0,step0+1)] <- NA
            #if(length(centerInd)) subCovs[1,centerNames[seq(2,2*length(centerInd),2)]] <- NA
            d$x=X[,1]
            d$y=X[,2]
            coordNames <- c("x","y")
          }
      } else if("step" %in% distnames){
        if(inputs$dist[["step"]] %in% stepdists){
          d$step[nrow(d)] <- NA
          tmpstep <- d$step
          tmpstep[which(d$level!=gsub("level","",coordLevel))] <- 0
          d$x=c(initialPosition[[zoo]][1],initialPosition[[zoo]][1]+cumsum(tmpstep)[-nrow(d)])
          d$y=rep(initialPosition[[zoo]][2],nrow(d))
          coordNames <- c("x","y")
        }    
      } else if(!is.null(mvnCoords)){
        d[[paste0(mvnCoords,".x")]] <- d[[mvnCoords]][,1]
        d[[paste0(mvnCoords,".y")]] <- d[[mvnCoords]][,2]
        if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) d[[paste0(mvnCoords,".z")]] <- d[[mvnCoords]][,3]
        d[[mvnCoords]] <- NULL
        coordNames <- paste0(mvnCoords,c(".x",".y"))
      } else {
        coordNames <- NULL
      }
      for(j in angleCovs[which(angleCovs %in% names(subCovs))])
        allCovs[cumNbObs[zoo]+1:nbObs,j] <- subCovs[,j]
      if(length(centerInd)) centerCovs[cumNbObs[zoo]+1:nbObs,] <- subCovs[,centerNames]
      if(length(centroidInd)) centroidCovs[cumNbObs[zoo]+1:nbObs,] <- subCovs[,centroidNames]
      data <- rbind(data,d)
    }
    
    if(nbCovs>0)
      data <- cbind(data,allCovs)
    
    if(nbSpatialCovs>0){
      colnames(allSpatialcovs)<-spatialcovnames
      for(j in spatialcovnames){
        if(any(raster::is.factor(spatialCovs[[j]]))){
          allSpatialcovs[[j]] <- factor(allSpatialcovs[[j]],levels=unique(unlist(raster::levels(spatialCovs[[j]]))))
        }
      }
      data <- cbind(data,allSpatialcovs)
    }
    
    if(length(centerInd)){
      data <- cbind(data,centerCovs)
      for(j in which(grepl(".angle",names(data)))){
        if(names(data[j]) %in% centerNames)
          class(data[[j]]) <- c(class(data[[j]]), "angle")
      }
    }
    
    if(length(centroidInd)){
      data <- cbind(data,centroidCovs)
      for(j in which(grepl(".angle",names(data)))){
        if(names(data[j]) %in% centroidNames)
          class(data[[j]]) <- c(class(data[[j]]), "angle")
      }
    }
    
    # include states sequence in the data
    if(states)
      data <- cbind(data,states=allStates)
    
    for(i in distnames){
      if(inputs$dist[[i]] %in% angledists)
        class(data[[i]]) <- c(class(data[[i]]), "angle")
    }
    
    for(i in angleCovs){
      class(data[[i]]) <- c(class(data[[i]]), "angle")
    }
    
    if(!is.null(coordNames)){
      attr(data,'coords') <- coordNames
      attr(data,"coordLevel") <- coordLevel
      data[which(data$level!=coordLevel),coordNames] <- NA
      #tmpNames <- colnames(data)[-which(colnames(data)=="ID")]
      #prepDat <- prepData(data,coordNames = paste0(mvnCoords,c(".x",".y")))
      #data$step <- prepDat$step
      #data$angle <- prepDat$angle
      #data$x <- prepDat$x
      #data$y <- prepDat$y
      #data <- data[,c("ID","step","angle",tmpNames,"x","y")]
    }
    
    # account for observation error (if any)
    out<-simObsData(momentuHierHMMData(data),lambda,errorEllipse,coordLevel)
    
    message("DONE")
    return(out)
  } else {
    simCount <- 0
    cat("Attempting to simulate tracks within spatial extent(s) of raster layers(s). Press 'esc' to force exit from 'simHierData'\n",sep="")
    while(simCount < retrySims){
      cat("\r    Attempt ",simCount+1," of ",retrySims,"...",sep="")
      tmp<-suppressMessages(tryCatch(simHierData(nbAnimals,hierStates,hierDist,
                                             Par,hierBeta,hierDelta,
                                             hierFormula,hierFormulaDelta,mixtures,formulaPi,
                                             covs,nbHierCovs,
                                             spatialCovs,
                                             zeroInflation,
                                             oneInflation,
                                             circularAngleMean,
                                             centers,
                                             centroids,
                                             angleCovs,
                                             obsPerLevel,
                                             initialPosition,
                                             DM,userBounds,workBounds,mvnCoords,
                                             model,states,
                                             retrySims=0,
                                             lambda,
                                             errorEllipse),error=function(e) e))
      if(inherits(tmp,"error")){
        if(grepl("Try expanding the extent of the raster",tmp)) simCount <- simCount+1
        else stop(tmp)
      } else {
        simCount <- retrySims
        cat("DONE\n")
        return(tmp)
      }
    }
    cat("FAILED\n")
    stop(tmp)
  }
}
