#' Simulation tool
#'
#' Simulates data from a (multivariate) continuous-time discrete-space hidden Markov model. Note that only a single state and cell transition can occur between observations. Movement data are assumed to be in Cartesian coordinates (not longitude/latitude) and can be generated with or without observation error attributable to location measurement error.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param dist A named list indicating the probability distributions of the data streams.
#' @param Par A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be in the order expected by the pdfs of \code{dist}. 
#' 
#' If \code{DM} is not specified for a given data stream, then \code{Par} 
#' is on the natural (i.e., real) scale of the parameters. However, if \code{DM} is specified for a given data stream, then 
#' \code{Par} must be on the working (i.e., beta) scale of the parameters, and the length of \code{Par} must match the number 
#' of columns in the design matrix. See details below.
#' @param beta Matrix of regression parameters for the transition rates.
#' @param delta Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}. If \code{formulaDelta} includes a formula, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix, where k is the number of covariates and the columns correspond to states 2:\code{nbStates}. 
#' @param formula Regression formula for the transition rate covariates. Default: \code{~1} (no covariate effect). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition rates.
#' @param formulaDelta Regression formula for the initial distribution. Default: \code{NULL} (no covariate effects and \code{delta} is specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then \code{delta} must be specified on the working scale.
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. Default: \code{NULL} (no covariate effects; both \code{beta$pi} and \code{fixPar$pi} are specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{beta$pi} and \code{fixPar$pi} are specified on the working scale.
#' Note that only the covariate values corresponding to the first time step for each individual ID are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param covs Covariate values to include in the simulated data, as a dataframe. The names of any covariates specified by \code{covs} can
#' be included in \code{formula} and/or \code{DM}. Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL} (the default), and specifying \code{nbCovs>0}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1' and 'cov2' for \code{nbCovs=2}) and can be included in \code{formula} and/or \code{DM}.
#' @param rast A raster object or raster stack object that will define the discrete-space grid cells for the CTMC movement path. \code{spatialCovs} and \code{spatialCovs.grad} must have the same extent, number of rows and columns, projection, resolution, and origin as \code{rast}.
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date'). In the \code{\link{momentuHMMData}} object returned by \code{prepCTDS}, covariates for the current position (e.g.\ for use in \code{formula} or \code{DM}) are named with a \code{.cur} suffix (e.g. \code{cov1.cur}).
#' @param spatialCovs.grad List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates, where a directional gradient is to be calculated internally using \code{\link[ctmcmove]{rast.grad}}. Gradient-based covariates specified by \code{spatialCovs.grad} are extracted from the raster 
#' layer(s) based on the location data (and the z values for a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs.grad} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{data} must include column(s) of the corresponding z value(s) for each observation (e.g., 'Date').
#' @param directions Integer. Either 4 (indicating a "Rook's neighborhood" of 4 neighboring grid cells) or 8 (indicating a "King's neighborhood" of 8 neighboring grid cells).
#' @param normalize.gradients	Logical. Default is FALSE. If TRUE, then all gradient covariates for \code{spatialCovs.grad} are normalized by dividing by the length of the gradient vector at each point.
#' @param grad.point.decreasing	Logical. If TRUE, then the gradient covariates are positive in the direction of decreasing values of the covariate. If FALSE, then the gradient covariates are positive in the direction of increasing values of the covariate (like a true gradient).
#' @param zero.idx Integer vector of the indices of raster cells that are not passable and should be excluded. These are cells where movement should be impossible. Default is zero.idx=integer().
#' @param obsPerAnimal Either the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated for each animal are uniformously picked from this interval. Alternatively, \code{obsPerAnimal} can be specified as
#' a list of length \code{nbAnimals} with each element providing the number of observations (if single value) or the bounds (if vector of two values) for each individual.
#' Default: \code{c(500,1500)}.
#' @param initialPosition 2-vector providing the x- and y-coordinates of the initial position for all animals. Alternatively, \code{initialPosition} can be specified as
#' a list of length \code{nbAnimals} with each element a 2-vector providing the x- and y-coordinates of the initial position for each individual.
#' Default: \code{c(0,0)}.  If \code{mvnCoord} corresponds to a data stream with ``mvnorm3'' or ''rw_mvnorm3'' probability distributions, then \code{initialPosition} must be composed of 3-vector(s) for the x-, y-, and z-coordinates.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of regression formulas or a ``pseudo'' design matrix.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. For each matrix, the first column pertains to the lower bound and the second column the upper bound.
#' For data streams, each element of \code{workBounds} should be a k x 2 matrix with the same name of the corresponding element of 
#' \code{Par}, where k is the number of parameters. For transition rate parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``beta'', where k=\code{length(beta)}. For initial distribution parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``delta'', where k=\code{length(delta)}.
#' \code{workBounds} is ignored for any given data stream unless \code{DM} is also specified.
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the state transition rate matrix. Default: NULL, in which case
#' the diagonal elements of the transition rate matrix are the reference.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerAnimal}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
#' \code{covs}, and \code{spatialCovs}. It is not appropriate to simulate movement data from a \code{model} that was fitted to latitude/longitude data (because \code{simData} assumes Cartesian coordinates).
#' @param matchModelObs If \code{model} is provided, logical indicating whether to match \code{nbAnimals}, \code{obsPerAnimal}, and observation times to the fitted model data. If \code{TRUE}, then \code{nbAnimals}, \code{obsPerAnimal}, and \code{lambda} are ignored. Default: \code{TRUE}.
#' @param states \code{TRUE} if the simulated states should be returned, \code{FALSE} otherwise (default).
# #' @param times time values (numeric or POSIX) over which the CTHMM will be simulated for all animals. Alternatively, \code{times} can be specified as a list of length \code{nbAnimals} with each element the time values (numeric or POSIX) for each individual.
#' @param lambda Observation rate. \code{lambda} is the rate parameter of the exponential distribution for the waiting times between successive observations, i.e., 
#' \code{1/lambda} is the expected time between successive location observations. Note that only a single state transition can occur between observations. If \code{model} is specified and \code{model$data} time column is of class \code{\link[base]{date-time}} or \code{\link[base]{date}}, \code{lambda} has the same units as the \code{Time.unit} argument in \code{\link{fitCTHMM}}. Default: 1.
#' @param errorEllipse List providing the upper bound for the semi-major axis (\code{M}; on scale of x- and y-coordinates), semi-minor axis (\code{m}; 
#' on scale of x- and y-coordinates), and orientation (\code{r}; in degrees) of location error ellipses. If \code{NULL} (the default), no location 
#' measurement error is simulated. If \code{errorEllipse} is specified, then each observed location is subject to bivariate normal errors as described 
#' in McClintock et al. (2015), where the components of the error ellipse for each location are randomly drawn from \code{runif(1,min(errorEllipse$M),max(errorEllipse$M))}, 
#' \code{runif(1,min(errorEllipse$m),max(errorEllipse$m))}, and \code{runif(1,min(errorEllipse$r),max(errorEllipse$r))}. If only a single value is provided for any of the 
#' error ellipse elements, then the corresponding component is fixed to this value for each location. Only coordinate data streams are subject to location measurement error;
#' any other data streams are observed without error.
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @param export Character vector of the names of any additional objects or functions in the global environment that are used in \code{DM}, \code{formula}, \code{formulaDelta}, and/or \code{formulaPi}. Only necessary if \code{ncores>1} so that the needed items will be exported to the workers.
#' 
#' @return If the simulated data have no measurement error (i.e., \code{errorEllipse=NULL}), a \code{\link{momentuHMMData}} object, 
#' i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{time}{Numeric time of each observed observation}
#' \item{...}{Data streams as specified by \code{dist}}
#' \item{x}{Either easting or longitude (if data streams include valid non-negative distribution for 'step')}
#' \item{y}{Either norting or latitude (if data streams include valid non-negative distribution for 'step')}
#' \item{...}{Covariates (if any)}
#' 
#' If simulated location data include measurement error (i.e., \code{errorEllipse!=NULL}), a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{time}{Numeric time of each observed (and missing) observation}
#' \item{z}{Categorial index indicating cell movement (where \code{z=(directions+1)} indicates no movement)}
#' \item{x}{Either easting or longitude observed location}
#' \item{y}{Either norting or latitude observed location}
#' \item{tau}{Time difference between consecutive observations}
#' \item{...}{Data streams that are not derived from location (if applicable)}
#' \item{...}{Covariates at true (\code{mux},\code{muy}) locations (if any) and neighboring cell locations (with suffixes indicating neighbor, e.g., \code{cov.1}, \code{cov.2}, ..., \code{cov.directions})}
#' \item{mux}{Either easting or longitude true location}
#' \item{muy}{Either norting or latitude true location}
#' \item{error_semimajor_axis}{error ellipse semi-major axis (if applicable)}
#' \item{error_semiminor_axis}{error ellipse semi-minor axis (if applicable)}
#' \item{error_ellipse_orientation}{error ellipse orientation (if applicable)}
#' \item{ln.sd.x}{log of the square root of the x-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{ln.sd.y}{log of the square root of the y-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{error.corr}{correlation term of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' @export
#'

simCTDS <- function(nbAnimals=1,nbStates=2,dist,
                     Par,beta=NULL,delta=NULL,
                     formula=~1,formulaDelta=NULL,mixtures=1,formulaPi=NULL,
                     covs=NULL,nbCovs=0,
                     rast,
                     spatialCovs=NULL,
                     spatialCovs.grad=NULL,
                     directions=4,
                     normalize.gradients=FALSE,
                     grad.point.decreasing=FALSE,
                     zero.idx = integer(),
                     #zeroInflation=NULL,
                     #oneInflation=NULL,
                     #circularAngleMean=NULL,
                     #centers=NULL,
                     #centroids=NULL,
                     #angleCovs=NULL,
                     obsPerAnimal=c(500,1500),
                     initialPosition=c(0,0),
                     DM=NULL,userBounds=NULL,workBounds=NULL,betaRef=NULL,#mvnCoords=NULL,
                     stateNames=NULL,
                     model=NULL,
                     matchModelObs=TRUE,
                     states=FALSE,
                     #retrySims=0,
                     #times=NULL,
                     lambda=1,
                     errorEllipse=NULL,
                     ncores=1,
                     export=NULL)
{
  
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
    if(inherits(model,"momentuHierHMM") | inherits(model,"hierarchical")) stop("model can not be a 'momentuHierHMM' or 'hierarchical' object")
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
      obsPerAnimal <- as.list(table(model$data$ID)+1)#+ifelse(rwInd,1,0))
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
    }
  } else {
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
  
  withCallingHandlers(out <- simData(nbAnimals,nbStates,dist,
                 Par,beta,delta,
                 formula,formulaDelta,mixtures,formulaPi,
                 covs,nbCovs,
                 spatialCovs,
                 zeroInflation=NULL,
                 oneInflation=NULL,
                 circularAngleMean=NULL,
                 centers=NULL,
                 centroids=NULL,
                 angleCovs=NULL,
                 obsPerAnimal,
                 initialPosition,
                 DM,userBounds,workBounds,betaRef,mvnCoords=NULL,stateNames,
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
  class(out) <- unique(append(c("momentuHMMData","ctds"),class(out)))
  if(inherits(out$time,"POSIXt")) attr(out,"Time.unit") <- Time.unit
  attr(out,"directions") <- directions
  attr(out,"coords") <- c("x","y")
  attr(out,"ctdsData") <- "z"
  attr(out,"normalize.gradients") <- normalize.gradients
  attr(out,"grad.point.decreasing") <- grad.point.decreasing
  attr(out,"CT") <- TRUE
  attr(out,"Time.name") <- Time.name
  return(out)
  
}
