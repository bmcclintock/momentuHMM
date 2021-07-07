#' Simulation tool
#'
#' Simulates data from a (multivariate) continuous-time hidden Markov model. Note that only a single state transition can occur between observations. Movement data are assumed to be in Cartesian coordinates (not longitude/latitude) and can be generated with or without observation error attributable to location measurement error.
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
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on any simulated location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{covs} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' The names of the raster layer(s) can be included in 
#' \code{formula} and/or \code{DM}.  Note that \code{simCTHMM} usually takes longer to generate simulated data when \code{spatialCovs} is specified.
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
#' @param mvnCoords Character string indicating the name of location data that are to be simulated using a multivariate normal distribution. For example, if \code{mu="rw_mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are intended to be location data, then \code{mvnCoords="mu"} needs to be specified in order for these data to be treated as such.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerAnimal}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
#' \code{covs}, and \code{spatialCovs}. It is not appropriate to simulate movement data from a \code{model} that was fitted to latitude/longitude data (because \code{simData} assumes Cartesian coordinates).
#' @param matchModelObs If \code{model} is provided, logical indicating whether to match \code{nbAnimals}, \code{obsPerAnimal}, and observation times to the fitted model data. If \code{TRUE}, then \code{nbAnimals}, \code{obsPerAnimal}, and \code{lambda} are ignored. Default: \code{TRUE}.
#' @param states \code{TRUE} if the simulated states should be returned, \code{FALSE} otherwise (default).
#' @param retrySims Number of times to attempt to simulate data within the spatial extent of \code{spatialCovs}. If \code{retrySims=0} (the default), an
#' error is returned if the simulated tracks(s) move beyond the extent(s) of the raster layer(s). Instead of relying on \code{retrySims}, in many cases
#' it might be better to simply expand the extent of the raster layer(s) and/or adjust the movement parameters. 
#' Ignored if \code{spatialCovs=NULL}.
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
#' \item{x}{Either easting or longitude observed location}
#' \item{y}{Either norting or latitude observed location}
#' \item{...}{Data streams that are not derived from location (if applicable)}
#' \item{...}{Covariates at true (\code{mux},\code{muy}) locations (if any)}
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

simCTHMM <- function(nbAnimals=1,nbStates=2,dist,
                     Par,beta=NULL,delta=NULL,
                     formula=~1,formulaDelta=NULL,mixtures=1,formulaPi=NULL,
                     covs=NULL,nbCovs=0,
                     spatialCovs=NULL,
                     #zeroInflation=NULL,
                     #oneInflation=NULL,
                     #circularAngleMean=NULL,
                     #centers=NULL,
                     #centroids=NULL,
                     #angleCovs=NULL,
                     obsPerAnimal=c(500,1500),
                     initialPosition=c(0,0),
                     DM=NULL,userBounds=NULL,workBounds=NULL,betaRef=NULL,mvnCoords=NULL,stateNames=NULL,
                     model=NULL,
                     matchModelObs=TRUE,
                     states=FALSE,
                     retrySims=0,
                     #times=NULL,
                     lambda=1,
                     errorEllipse=NULL,
                     ncores=1)
{
  
  if(!is.null(model)){
    if(is.miHMM(model)){
      model <- model$miSum
    } 
    if(inherits(model,"momentuHierHMM") | inherits(model,"hierarchical")) stop("model can not be a 'momentuHierHMM' or 'hierarchical' object")
    if(!inherits(model,"CTHMM")) stop("model must be of class 'CTHMM'")
    attributes(model)$class <- attributes(model)$class[which(!attributes(model)$class %in% "CTHMM")]
    mvnCoords <- model$conditions$mvnCoords
    if(matchModelObs){
      rwInd <- any(unlist(lapply(model$conditions$dist,function(x) x %in% rwdists)))
      nbAnimals <- length(unique(model$data$ID))
      obsPerAnimal <- as.list(table(model$data$ID)+ifelse(rwInd,1,0))
      lambda <- list()
      for(zoo in 1:nbAnimals){
       lambda[[zoo]] <- model$data[[model$conditions$Time.name]][which(model$data$ID==unique(model$data$ID)[zoo])]
       if(rwInd) lambda[[zoo]] <- c(lambda[[zoo]],tail(lambda[[zoo]],1)+model$data$dt[tail(which(model$data$ID==unique(model$data$ID)[zoo]),1)])
       if(inherits(lambda[[zoo]] ,"POSIXt")) attr(lambda[[zoo]],"units") <- model$conditions$Time.unit
      }
    }
  } else {
    for(i in names(dist)){
      if(!dist[[i]] %in% CTHMMdists) stop("Sorry, currently simCTHMM only supports the following distributions: ",paste0(CTHMMdists,sep=", "))
    }
    if(!is.null(lambda)){
      if(length(lambda)>1 || lambda<=0) stop('lambda must be a scalar and >0')
    } else stop("lambda cannot be NULL")
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
                 DM,userBounds,workBounds,betaRef,mvnCoords,stateNames,
                 model,states,
                 retrySims,
                 lambda,
                 errorEllipse,
                 ncores,
                 CT=TRUE),warning=muffleCTwarning)
  
  out<- out[,c("ID","time",colnames(out)[which(!colnames(out) %in% c("ID","time"))])]
  if(!is.null(mvnCoords)){
    attr(out,'coords') <- paste0(mvnCoords,c(".x",".y"))
  }
  return(out)
  
}