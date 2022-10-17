#' Simulation tool
#' 
#' Simulates data from an approximate (multivariate) hierarchical continuous-time hidden Markov model. Note that the state active at observation time \emph{t} determines the (state-dependent) observation distribution from observation time \emph{t} to time \emph{t+1} (see Details) and time-varying covariates are assumed piece-wise constant between observations. Movement data are assumed to be in Cartesian coordinates (not longitude/latitude) and can be generated with or without observation error attributable to location measurement error.
# #' @rdname simCTHMM
#' @param nbAnimals Number of observed individuals to simulate.
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. See details.
#' @param Par A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be in the order expected by the pdfs of \code{dist}. 
#' 
#' If \code{DM} is not specified for a given data stream, then \code{Par} 
#' is on the natural (i.e., real) scale of the parameters. However, if \code{DM} is specified for a given data stream, then 
#' \code{Par} must be on the working (i.e., beta) scale of the parameters, and the length of \code{Par} must match the number 
#' of columns in the design matrix. See details below.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (only hierarchical-level effects, with no covariate effects).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy. See details.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. Default: \code{NULL} (no covariate effects; both \code{beta$pi} and \code{fixPar$pi} are specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{beta$pi} and \code{fixPar$pi} are specified on the working scale.
#' Note that only the covariate values corresponding to the first time step for each individual ID are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param covs Covariate values to include in the simulated data, as a dataframe. The names of any covariates specified by \code{covs} can
#' be included in \code{formula} and/or \code{DM}. Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL} (the default), and specifying \code{nbCovs>0}.
#' @param nbHierCovs A hierarchical data structure \code{\link[data.tree]{Node}} for the number of covariates ('nbCovs') to simulate for each level of the hierarchy (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1.1' and 'cov1.2' for \code{nbHierCovs$level1$nbCovs=2}) and can be included in \code{hierFormula} and/or \code{DM}.
#' @param spatialCovs List of \code{\link[raster]{raster}} objects for spatio-temporally referenced covariates. Covariates specified by \code{spatialCovs} are extracted from the raster 
#' layer(s) based on any simulated location data (and the z values for a raster \code{\link[raster]{stack}} 
#' or \code{\link[raster]{brick}}) for each time step.  If an element of \code{spatialCovs} is a raster \code{\link[raster]{stack}} or \code{\link[raster]{brick}}, 
#' then z values must be set using \code{raster::setZ} and \code{covs} must include column(s) of the corresponding z value(s) for each observation (e.g., 'time').
#' The names of the raster layer(s) can be included in 
#' \code{formula} and/or \code{DM}.  Note that \code{simHierCTHMM} usually takes longer to generate simulated data when \code{spatialCovs} is specified.
#' @param obsPerLevel A hierarchical data structure \code{\link[data.tree]{Node}} indicating the number of observations for each level of the hierarchy ('obs'). For each level, the 'obs' field can either be the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated per level for each animal are uniformously picked from this interval. Alternatively, \code{obsPerLevel} can be specified as
#' a list of length \code{nbAnimals} with each element providing the hierarchical data structure for the number of observations for each level of the hierarchy for each animal, where the 'obs' field can either be the number of observations (if single value) or the bounds of the number of observations (if vector of two values) for each individual.
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
#' @param mvnCoords Character string indicating the name of location data that are to be simulated using a multivariate normal distribution. For example, if \code{mu="rw_mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are intended to be location data, then \code{mvnCoords="mu"} needs to be specified in order for these data to be treated as such.
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerLevel}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
#' \code{covs}, and \code{spatialCovs}. It is not appropriate to simulate movement data from a \code{model} that was fitted to latitude/longitude data (because \code{simHierCTHMM} assumes Cartesian coordinates).
#' @param matchModelObs If \code{model} is provided, logical indicating whether to match \code{nbAnimals}, \code{obsPerLevel}, and observation times to the fitted model data. If \code{TRUE}, then \code{nbAnimals}, \code{obsPerLevel}, and \code{lambda} are ignored. Default: \code{TRUE}.
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
#' @param export Character vector of the names of any additional objects or functions in the global environment that are used in \code{DM}, \code{formula}, \code{formulaDelta}, and/or \code{formulaPi}. Only necessary if \code{ncores>1} so that the needed items will be exported to the workers.
#' @param gradient Logical indicating whether or not to calculate gradients of \code{spatialCovs} using bilinear interpolation (e.g. for inclusion in potential functions). Default: \code{FALSE}. If \code{TRUE}, the gradients are returned with ``\code{.x}'' (easting gradient) and ``\code{.y}'' (northing gradient) suffixes added to the names of \code{spatialCovs}. For example, if \code{cov1} is the name of a spatial covariate, then the returned \code{\link{momentuHMMData}} object will include the fields ``\code{cov1.x}'' and ``\code{cov1.y}''.
#'
#' @details \itemize{
#' \item Importantly, \code{simHierCTHMM} assumes the state active at observation time \emph{t} determines the (state-dependent) observation distribution from observation time \emph{t} to time \emph{t+1} (i.e. state switches only occur at the times of observations). 
#' The snapshot property applies to all data stream distributions (i.e. observations are "instantaneous") except for the (multivariate) normal random walk (\code{rw_norm}, \code{rw_mvnorm2}, \code{rw_mvnorm3}) and Poisson (\code{pois}) distributions. For these particular distributions, the observed data are not "instantaneous"; they depend on the time interval between observations \eqn{(\Delta_t)}.
#' It is critical that the frequency of observations (specified by \code{lambda}) is high relative to the serial correlation in the hidden state process (specified by \code{hierBeta}) in order for this discrete-time approximation to be reasonably accurate.
#' 
#' \item Model specification in \code{simHierCTHMM} is similar to \code{\link{simCTHMM}} except that instead of simply specifying the number of states (\code{nbStates}), distributions (\code{dist}), observations (\code{obsPerAnimal}), covariates (\code{nbCovs}), and a single t.p.m. formula (\code{formula}), the \code{hierStates} argument specifies the hierarchical nature of the states,
#' the \code{hierDist} argument specifies the hierarchical nature of the data streams, the \code{obsPerLevel} argument specifies the number of observations for each level of the hierarchy, the \code{nbHierCovs} argument specifies the number of covariates for each level of the hierarchy, and the \code{hierFormula} argument specifies a t.p.m. formula for each level of the hierarchy.
#' All of the hierarchical arguments in \code{simHierCTHMM} are specified as \code{\link[data.tree]{Node}} objects from the \code{\link[data.tree]{data.tree}} package.
#' 
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerLevel', then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' 
#' \item In \code{hierDelta}, 'delta' must be specified
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
# #' @importFrom data.tree Node Get Aggregate isLeaf Clone
#' @importFrom doParallel registerDoParallel stopImplicitCluster

simHierCTHMM <- function(nbAnimals=1,hierStates,hierDist,
                        Par,hierBeta=NULL,hierDelta=NULL,
                        hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,
                        covs=NULL,nbHierCovs=NULL,
                        spatialCovs=NULL,
                        obsPerLevel,
                        initialPosition=c(0,0),
                        DM=NULL,userBounds=NULL,workBounds=NULL,mvnCoords=NULL,
                        model=NULL,
                        matchModelObs = TRUE,
                        states=FALSE,
                        retrySims=0,
                        lambda=1,
                        errorEllipse=NULL,
                        ncores=1,
                        export=NULL,
                        gradient=FALSE)
{
  
  installDataTree()
  
  if(!is.null(model)){
    if(is.miHMM(model)){
      model <- model$miSum
    } 
    if(!(inherits(model,"momentuHierHMM") | inherits(model,"hierarchical"))) stop("model must be a 'momentuHierHMM' and/or 'hierarchical' object; use simCTHMM instead")
    if(!inherits(model,"CTHMM")) stop("model must be of class 'CTHMM'; use simData instead")
    attributes(model)$class <- attributes(model)$class[which(!attributes(model)$class %in% "CTHMM")]
    mvnCoords <- model$conditions$mvnCoords
    if(matchModelObs){
      rwInd <- any(unlist(lapply(model$conditions$dist,function(x) x %in% rwdists)))
      if(nbAnimals!=1) warning("'nbAnimals' is ignored when 'matchModelObs' is TRUE")
      nbAnimals <- length(unique(model$data$ID))
      if(!missing(obsPerLevel)) warning("'obsPerLevel' is ignored when 'matchModelObs' is TRUE")
      obsPerLevel <- NULL
      if(lambda!=1) warning("'lambda' is ignored when 'matchModelObs' is TRUE")
      lambda <- list()
      for(zoo in 1:nbAnimals){
        lambda[[zoo]] <- model$data[[model$conditions$Time.name]][which(model$data$ID==unique(model$data$ID)[zoo])]
        if(rwInd) lambda[[zoo]] <- c(lambda[[zoo]],tail(lambda[[zoo]],1)+model$data$dt[tail(which(model$data$ID==unique(model$data$ID)[zoo]),1)])
        if(inherits(lambda[[zoo]] ,"POSIXt")) attr(lambda[[zoo]],"units") <- model$conditions$Time.unit
      }
    } else {
      if(missing(obsPerLevel)) stop('argument "obsPerLevel" is missing, with no default')
      if(is.null(obsPerLevel)) stop("'obsPerLevel' cannot be NULL when 'matchModelObs' is FALSE")
    }
  } else {
    dist <- formatHierHMM(NULL,hierStates=hierStates,hierDist=hierDist)$dist
    for(i in names(dist)){
      if(!dist[[i]] %in% CTHMMdists) stop("Sorry, currently simHierCTHMM only supports the following distributions: ",paste0(CTHMMdists,sep=", "))
    }
    if(!is.null(lambda)){
      if(length(lambda)>1 || lambda<=0) stop('lambda must be a scalar and >0')
    } else stop("lambda cannot be NULL")
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
                                         DM,userBounds,workBounds,mvnCoords,
                                         model,states,
                                         retrySims,
                                         lambda,
                                         errorEllipse,
                                         ncores,
                                         export,
                                         gradient,
                                         CT=TRUE),warning=muffleCTwarning)
  
  coordLevel <- attr(out,"coordLevel")
  out<- out[,c("ID","time",colnames(out)[which(!colnames(out) %in% c("ID","time"))])]
  if(!is.null(mvnCoords)){
    attr(out,'coords') <- paste0(mvnCoords,c(".x",".y"))
    attr(out,"coordLevel") <- coordLevel
  }
  attr(out,"CT") <- TRUE
  attr(out,"Time.name") <- "time"
  attr(out,"gradient") <- ifelse(!is.null(model),isTRUE(attr(model$data,"gradient")),gradient)
  return(out)
}
