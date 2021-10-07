
#' @rdname simCTDS
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
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal' (or 'obsPerLevel' for \code{simHierCTDS}), then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' 
#' \item For \code{simCTDS}, when covariates are not included in \code{formulaDelta} (i.e. \code{formulaDelta=NULL}), then \code{delta} is specified as a vector of length \code{nbStates} that 
#' sums to 1.  When covariates are included in \code{formulaDelta}, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#' 
#' \item For \code{simHierCTDS}, \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. 
#' }
#' @export
#' @importFrom stats rnorm runif rmultinom step terms.formula
# @importFrom raster cellFromXY getValues
#' @importFrom CircStats rvm
#' @importFrom Brobdingnag as.brob sum
#' @importFrom mvtnorm rmvnorm
# #' @importFrom data.tree Node Get Aggregate isLeaf Clone
#' @importFrom doParallel registerDoParallel stopImplicitCluster

simHierCTDS <- function(nbAnimals=1,hierStates,hierDist,
                         Par,hierBeta=NULL,hierDelta=NULL,
                         hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,
                         covs=NULL,nbHierCovs=NULL,
                         rast,
                         spatialCovs=NULL,
                         spatialCovs.grad=NULL,
                         directions=4,
                         normalize.gradients=FALSE,
                         grad.point.decreasing=FALSE,
                         zero.idx = integer(),
                         obsPerLevel,
                         initialPosition=c(0,0),
                         DM=NULL,userBounds=NULL,workBounds=NULL,#mvnCoords=NULL,
                         model=NULL,
                         matchModelObs = TRUE,
                         states=FALSE,
                         #retrySims=0,
                         lambda=1,
                         errorEllipse=NULL,
                         ncores=1,
                         export=NULL)
{
  
  installDataTree()
  
  if (!requireNamespace("ctmcmove", quietly = TRUE)) {
    stop("Package \"ctmcmove\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package \"raster\" needed for spatial covariates. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("extraDistr", quietly = TRUE)) {
    stop("Package \"extraDistr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if(!is.null(model)){
    if(is.miHMM(model)){
      model <- model$miSum
    } 
    if(!inherits(model,"momentuHierHMM") & !inherits(model,"hierarchical")) stop("model must be a 'momentuHierHMM' or 'hierarchical' object")
    if(!inherits(model,"CTHMM")) stop("model must be of class 'CTHMM'")
    if(!inherits(model,"ctds")) stop("model must be of class 'ctds'")
    attributes(model)$class <- attributes(model)$class[which(!attributes(model)$class %in% "CTHMM")]
    #mvnCoords <- model$conditions$mvnCoords
    directions <- attr(model$data,"directions")
    normalize.gradients <- attr(model$data,"normalize.gradients")
    grad.point.decreasing <- attr(model$data,"grad.point.decreasing")
    if(matchModelObs){
      #rwInd <- any(unlist(lapply(model$conditions$dist,function(x) x %in% rwdists)))
      nbAnimals <- length(unique(model$data$ID))
      if(!missing(obsPerLevel)) warning("'obsPerLevel' is ignored when 'matchModelObs' is TRUE")
      obsPerLevel <- NULL
      lambda <- list()
      Time.name <- model$conditions$Time.name
      Time.unit <- model$conditions$Time.unit
      for(zoo in 1:nbAnimals){
        iTime <- model$data[[Time.name]][which(model$data$ID==unique(model$data$ID)[zoo])]
        lambda[[zoo]] <- iTime #as.numeric(iTime-min(iTime),units=Time.unit)
        #if(rwInd) lambda[[zoo]] <- c(lambda[[zoo]],tail(lambda[[zoo]],1)+model$data$dt[tail(which(model$data$ID==unique(model$data$ID)[zoo]),1)])
        if(inherits(lambda[[zoo]] ,"POSIXt")) {
          attr(lambda[[zoo]],"units") <- Time.unit
          lambda[[zoo]] <- c(lambda[[zoo]],tail(lambda[[zoo]],1)+lubridate::duration(tail(model$data$dt[which(model$data$ID==unique(model$data$ID)[zoo])],1),units=Time.unit))
        } else {
          lambda[[zoo]] <- c(lambda[[zoo]],tail(lambda[[zoo]],1)+tail(model$data$dt[which(model$data$ID==unique(model$data$ID)[zoo])],1))
        }
      }
    } else {
      if(missing(obsPerLevel)) stop('argument "obsPerLevel" is missing, with no default')
      if(is.null(obsPerLevel)) stop("'obsPerLevel' cannot be NULL when 'matchModelObs' is FALSE")
    }
  } else {
    dist <- formatHierHMM(NULL,hierStates=hierStates,hierDist=hierDist)$dist
    if(!("z" %in% names(dist)) || dist$z!="ctds") stop("'ctds' data stream named 'z' must be included")
    if(sum(unlist(dist)=="ctds")>1) stop("Only one 'ctds' data stream is allowed")
    for(i in names(dist)){
      if(!dist[[i]] %in% CTHMMdists) stop("Sorry, currently simCTDS only supports the following distributions: ",paste0(CTHMMdists,sep=", "))
    }
    if(!(directions %in% c(4,8))) stop("'directions' must be 4 or 8")
    if(!is.null(lambda)){
      if(length(lambda)>1 || lambda<=0) stop('lambda must be a scalar and >0')
    } else stop("lambda cannot be NULL")
    Time.name = "time"
  }
  
  if(missing(rast) || !inherits(rast,"RasterLayer")) stop("'rast' must be provided as a RasterLayer")
  raster::values(rast) <- 1
  
  if(is.null(spatialCovs) & is.null(spatialCovs.grad)) stop("spatialCovs and/or spatialCovs.grad rasters must be provided")
  
  if(!is.null(spatialCovs)){
    for(i in names(spatialCovs)){
      raster::compareRaster(rast,spatialCovs[[i]])
      if(raster::nlayers(spatialCovs[[i]])>1){
        zname <- names(attributes(spatialCovs[[i]])$z)
        zvalues <- raster::getZ(spatialCovs[[i]])
        spatialCovs[[i]] <- raster::extend(spatialCovs[[i]],1,value=0)
        spatialCovs[[i]] <- raster::setZ(spatialCovs[[i]],zvalues,zname)
      } else spatialCovs[[i]] <- raster::extend(spatialCovs[[i]],1,value=0)
      attr(spatialCovs[[i]],"nograd") <- TRUE
    }
  }
  
  if(!is.null(spatialCovs.grad)){
    for(i in names(spatialCovs.grad)){
      raster::compareRaster(rast,spatialCovs.grad[[i]])
    }
  }
  
  rast <- raster::extend(rast,1)
  
  if(!is.null(spatialCovs.grad)){
    if(is.null(spatialCovs)) spatialCovs <- list()
    else if(any(names(spatialCovs) %in% names(spatialCovs.grad))) stop("spatialCovs and spatialCovs.grad names must be unique")
    if (length(zero.idx) > 0) {
      notzero.rast <- spatialCovs.grad[[names(spatialCovs.grad)[1]]]
      raster::values(notzero.rast) <- NA
      notzero.idx = (1:raster::ncell(notzero.rast))[-zero.idx]
      raster::values(notzero.rast)[notzero.idx] <- 1
      notzero.rast <- raster::extend(notzero.rast,1,value=1)
      zero.idx <- which(is.na(raster::values(notzero.rast)))
    }
    adj <- raster::adjacent(rast,which(!is.na(raster::values(rast))), pairs = TRUE, sorted = TRUE, id = TRUE, directions = directions)
    xy.cell = raster::xyFromCell(rast, adj[,2])
    xy.adj = raster::xyFromCell(rast, adj[,3])
    v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 1, sum))
    gradMat <- get.grad(directions=directions, start.cells =  adj[,1] ,v.adj = v.adj, spatialCovs.grad = spatialCovs.grad, normalize.gradients = normalize.gradients, grad.point.decreasing = grad.point.decreasing, gradMat = FALSE)
    #if(any(!is.finite(gradMat))){
    #  warning("some gradients were not finite and were set to zero")
    #  gradMat[which(!is.finite(gradMat))] <- 0
    #}
    if(length(zero.idx)) raster::values(rast)[zero.idx] <- 0
    for(i in names(spatialCovs.grad)){
      for(j in 1:directions){
        spatialCovs[[paste0(i,".",j)]] <- spatialCovs.grad[[i]]
        names(spatialCovs[[paste0(i,".",j)]]) <- paste0(names(spatialCovs[[paste0(i,".",j)]]),".",j)
        if(raster::nlayers(spatialCovs[[paste0(i,".",j)]])>1){
          for(l in 1:raster::nlayers(spatialCovs[[paste0(i,".",j)]])){
            raster::values(spatialCovs[[paste0(i,".",j)]][[l]]) <- gradMat[[i]][seq(j,raster::ncell(spatialCovs.grad[[i]])*directions,directions),l]
          }
          spatialCovs[[paste0(i,".",j)]] <- raster::extend(spatialCovs[[paste0(i,".",j)]],1,value=0)
          zname <- names(attributes(spatialCovs.grad[[i]])$z)
          zvalues <- raster::getZ(spatialCovs.grad[[i]])
          spatialCovs[[paste0(i,".",j)]] <- raster::setZ(spatialCovs[[paste0(i,".",j)]],zvalues,zname)
        } else {
          raster::values(spatialCovs[[paste0(i,".",j)]]) <- gradMat[[i]][seq(j,raster::ncell(spatialCovs.grad[[i]])*directions,directions)]
          spatialCovs[[paste0(i,".",j)]] <- raster::extend(spatialCovs[[paste0(i,".",j)]],1,value=0)
        }
        attr(spatialCovs[[paste0(i,".",j)]],"grad") <- TRUE
      }
    }
  }
  
  withCallingHandlers(out <- simHierData(nbAnimals,hierStates,hierDist,
                                     Par,hierBeta,hierDelta,
                                     hierFormula,hierFormulaDelta,mixtures,formulaPi,
                                     covs,nbHierCovs,
                                     spatialCovs,
                                     zeroInflation=NULL,
                                     oneInflation=NULL,
                                     circularAngleMean=NULL,
                                     centers=NULL,
                                     centroids=NULL,
                                     angleCovs=NULL,
                                     obsPerLevel,
                                     initialPosition,
                                     DM,userBounds,workBounds,mvnCoords=NULL,
                                     model,states,
                                     retrySims=0,
                                     lambda,
                                     errorEllipse,
                                     ncores,
                                     export=export,
                                     CT=TRUE,
                                     ctds=TRUE,
                                     rast=rast,
                                     directions=directions),warning=muffleCTDSwarning)
  
  out<- out[,c("ID","time",colnames(out)[which(!colnames(out) %in% c("ID","time"))])]
  class(out) <- unique(append(c("momentuHierHMMData","ctds"),class(out)))
  if(inherits(out$time,"POSIXt")) attr(out,"Time.unit") <- Time.unit
  attr(out,"directions") <- directions
  attr(out,"coords") <- c("x","y")
  attr(out,"ctdsData") <- "z"
  attr(out,"normalize.gradients") <- normalize.gradients
  attr(out,"grad.point.decreasing") <- grad.point.decreasing
  attr(out,"CT") <- TRUE
  attr(out,"Time.name") <- "time"
  return(out)
}
