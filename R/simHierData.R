
#' Hierarchical simulation tool
#'
#' Simulates data from a (multivariate) hierarchical hidden Markov model. Movement data can be generated with or without observation error attributable to temporal irregularity or location measurement error.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states.  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams. Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. See details.
#' @param Par A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be in the order expected by the pdfs of \code{dist}, 
#' and any zero-mass and/or one-mass parameters should be the last (if both are present, then zero-mass parameters must preceed one-mass parameters). 
#' 
#' If \code{DM} is not specified for a given data stream, then \code{Par} 
#' is on the natural (i.e., real) scale of the parameters. However, if \code{DM} is specified for a given data stream, then 
#' \code{Par} must be on the working (i.e., beta) scale of the parameters, and the length of \code{Par} must match the number 
#' of columns in the design matrix. See details below.
#' @param beta Matrix of regression parameters for the transition probabilities (more information
#' in "Details").
#' @param delta Initial value for the initial distribution at the top level of the hierarchy. Default \code{NULL}: all top-level states have equal initial probability.
#' \code{delta} must be specified as a k x (\code{nbStates}-1) matrix, where k is the number of covariates and the columns correspond to states 2:\code{nbStates}. See details below.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy. Default: \code{NULL} (no covariate effects). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities.
#' @param formulaDelta Regression formula for the initial distribution at the top level of the hierarchy. Default: \code{~1} (both \code{delta} and \code{fixPar$delta} are specified on the working scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}).
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. Default: \code{NULL} (no covariate effects; both \code{beta0$pi} and \code{fixPar$pi} are specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{beta0$pi} and \code{fixPar$pi} are specified on the working scale.
#' Note that only the covariate values corresponding to the first time step for each individual ID are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param covs Covariate values to include in the simulated data, as a dataframe. The names of any covariates specified by \code{covs} can
#' be included in \code{formula} and/or \code{DM}. Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL} (the default), and specifying \code{nbCovs>0}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1' and 'cov2' for \code{nbCovs=2}) and can be included in \code{formula} and/or \code{DM}.
#' @param spatialCovs List of \code{\link[raster]{RasterLayer-class}} objects for spatially-referenced covariates. Covariates specified by \code{spatialCovs} are
#' extracted from the raster layer(s) based on the simulated location data for each time step (if applicable). The names of the raster layer(s) can be included in 
#' \code{formula} and/or \code{DM}.  Note that \code{simHierData} usually takes longer to generate simulated data when \code{spatialCovs} is specified.
#' @param zeroInflation A named list of logicals indicating whether the probability distributions of the data streams should be zero-inflated. If \code{zeroInflation} is \code{TRUE} 
#' for a given data stream, then values for the zero-mass parameters should be
#' included in the corresponding element of \code{Par}.
#' @param oneInflation A named list of logicals indicating whether the probability distributions of the data streams should be one-inflated. If \code{oneInflation} is \code{TRUE} 
#' for a given data stream, then values for the one-mass parameters should be
#' included in the corresponding element of \code{Par}.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  For example, 
#' \code{circularAngleMean=list(angle=TRUE)} indicates the angle mean is be estimated for 'angle' using circular-circular 
#' regression.  Whenever circular-circular regression is used for an angular data stream, a corresponding design matrix (\code{DM}) 
#' must be specified for the data stream, and the previous movement direction (i.e., a turning angle of zero) is automatically used 
#' as the reference angle (i.e., the intercept). Default is \code{NULL}, which assumes circular-linear regression is 
#' used for any angular distributions. Any \code{circularAngleMean} elements 
#' corresponding to data streams that do not have angular distributions are ignored.
#' \code{circularAngleMean} is also ignored for any 'vmConsensus' data streams (because the consensus model is a circular-circular regression model).
#' 
#' Alternatively, \code{circularAngleMean} can be specified as a numeric scalar, where the value specifies the coefficient for the reference angle (i.e., directional persistence) term in the circular-circular regression model. For example, setting \code{circularAngleMean} to \code{0} specifies a 
#' circular-circular regression model with no directional persistence term (thus specifying a biased random walk instead of a biased correlated random walk). Setting \code{circularAngleMean} to 1 is equivalent to setting it to TRUE, i.e., a circular-circular regression model with a coefficient of 1 for the directional persistence reference angle.
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on the simulated location data. These distance and angle 
#' covariates can be included in \code{formula} and \code{DM} using the row names of \code{centers}.  If no row names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'center1.dist', 'center1.angle', 'center2.dist', 'center2.angle'); otherwise the covariate names are derived from the row names
#' of \code{centers} as \code{paste0(rep(rownames(centers),each=2),c(".dist",".angle"))}. Note that the angle covariates for each activity center are calculated relative to 
#' the previous movement direction instead of standard directions relative to the x-axis; this is to allow turning angles to be simulated as a function of these covariates using circular-circular regression.
#' @param centroids List where each element is a data frame with rows providing the x-coordinates ('x') and y-coordinates ('y) for centroids (i.e., dynamic activity centers where the coordinates can change for each time step)
#' from which distance and angle covariates will be calculated based on the simulated location data. These distance and angle 
#' covariates can be included in \code{hierFormula} and \code{DM} using the names of \code{centroids}.  If no list names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'centroid1.dist', 'centroid1.angle', 'centroid2.dist', 'centroid2.angle'); otherwise the covariate names are derived from the list names
#' of \code{centroids} as \code{paste0(rep(names(centroids),each=2),c(".dist",".angle"))}. Note that the angle covariates for each centroid are calculated relative to 
#' the previous movement direction instead of standard directions relative to the x-axis; this is to allow turning angles to be simulated as a function of these covariates using circular-circular regression.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{covs} or \code{spatialCovs} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' using \code{\link{circAngles}}.
#' @param obsPerLevel A hierarchical data structure \code{\link[data.tree]{Node}} indicating the number of observations for each level of the hierarchy. For each level, the 'obs' field can either be the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated per level for each animal are uniformously picked from this interval. Alternatively, \code{obsPerLevel} can be specified as
#' a list of length \code{nbAnimals} with each element providing the hierarchical data structure for the number of observations for each level of the hierarchy for each animal, where the 'obs' field can either be the number of observations (if single value) or the bounds of the number of observations (if vector of two values) for each individual.
#' @param initialPosition 2-vector providing the x- and y-coordinates of the initial position for all animals. Alternatively, \code{initialPosition} can be specified as
#' a list of length \code{nbAnimals} with each element a 2-vector providing the x- and y-coordinates of the initial position for each individual.
#' Default: \code{c(0,0)}.  If \code{mvnCoord} corresponds to a data stream with ``mvnorm3'' or ''rw_mvnorm3'' probability distributions, then \code{initialPosition} must be composed of 3-vector(s) for the x-, y-, and z-coordinates.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of regression formulas or a ``pseudo'' design matrix.  For example, for a 2-state 
#' model using the gamma distribution for a data stream named 'step', \code{DM=list(step=list(mean=~cov1, sd=~1))} specifies the mean 
#' parameters as a function of the covariate 'cov1' for each state.  This model could equivalently be specified as a 4x6 ``pseudo'' design matrix using 
#' character strings for the covariate: 
#' \code{DM=list(step=matrix(c(1,0,0,0,'cov1',0,0,0,0,1,0,0,0,'cov1',0,0,0,0,1,0,0,0,0,1),4,6))}
#' where the 4 rows correspond to the state-dependent paramaters (mean_1,mean_2,sd_1,sd_2) and the 6 columns correspond to the regression 
#' coefficients. 
#' 
#' Design matrices specified using formulas allow standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}).  Special formula functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' \code{angleFormula(cov,strength,by)} for the angle mean of circular-circular regression models, and state-specific formulas (see details). Any formula terms that are not state-specific are included on the parameters for all \code{nbStates} states.
#' @param cons Deprecated. An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. While there could be other uses, primarily intended to constrain specific parameters to be positive. For example, 
#' \code{cons=list(step=c(1,2,1,1))} raises the second parameter to the second power. Default=NULL, which simply raises all parameters to 
#' the power of 1. \code{cons} is ignored for any given data stream unless \code{DM} is specified.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2,dimnames=list(c("mean_1",
#' "mean_2","concentration_1","concentration_2"))))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. For each matrix, the first column pertains to the lower bound and the second column the upper bound.
#' For data streams, each element of \code{workBounds} should be a k x 2 matrix with the same name of the corresponding element of 
#' \code{Par}, where k is the number of parameters. For transition probability parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``beta'', where k=\code{length(beta)}. For initial distribution parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``delta'', where k=\code{length(delta)}.
#' \code{workBounds} is ignored for any given data stream unless \code{DM} is also specified.
#' @param workcons Deprecated. An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. Warning: use of \code{workcons} is recommended only for advanced users implementing unusual parameter constraints 
#' through a combination of \code{DM}, \code{cons}, and \code{workcons}. \code{workcons} is ignored for any given data stream unless \code{DM} is specified.
#' @param mvnCoords Character string indicating the name of location data that are to be simulated using a multivariate normal distribution. For example, if \code{mu="rw_mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are intended to be location data, then \code{mvnCoords="mu"} needs to be specified in order for these data to be treated as such.
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerLevel}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
#' \code{covs}, \code{spatialCovs}, \code{centers}, and \code{centroids}.
#' @param states \code{TRUE} if the simulated states should be returned, \code{FALSE} otherwise (default).
#' @param retrySims Number of times to attempt to simulate data within the spatial extent of \code{spatialCovs}. If \code{retrySims=0} (the default), an
#' error is returned if the simulated tracks(s) move beyond the extent(s) of the raster layer(s). Instead of relying on \code{retrySims}, in many cases
#' it might be better to simply expand the extent of the raster layer(s) and/or adjust the step length and turning angle probability distributions. 
#' Ignored if \code{spatialCovs=NULL}.
#' @param lambda Observation rate for location data. If \code{NULL} (the default), location data are obtained at regular intervals. Otherwise 
#' \code{lambda} is the rate parameter of the exponential distribution for the waiting times between successive location observations, i.e., 
#' \code{1/lambda} is the expected time between successive location observations. Only the 'step' and 'angle' data streams are subject to temporal irregularity;
#' any other data streams are observed at temporally-regular intervals.  Ignored unless a valid distribution for the 'step' data stream is specified.
#' @param errorEllipse List providing the upper bound for the semi-major axis (\code{M}; on scale of x- and y-coordinates), semi-minor axis (\code{m}; 
#' on scale of x- and y-coordinates), and orientation (\code{r}; in degrees) of location error ellipses. If \code{NULL} (the default), no location 
#' measurement error is simulated. If \code{errorEllipse} is specified, then each observed location is subject to bivariate normal errors as described 
#' in McClintock et al. (2015), where the components of the error ellipse for each location are randomly drawn from \code{runif(1,min(errorEllipse$M),max(errorEllipse$M))}, 
#' \code{runif(1,min(errorEllipse$m),max(errorEllipse$m))}, and \code{runif(1,min(errorEllipse$r),max(errorEllipse$r))}. If only a single value is provided for any of the 
#' error ellipse elements, then the corresponding component is fixed to this value for each location. Only the 'step' and 'angle' data streams are subject to location measurement error;
#' any other data streams are observed without error.  Ignored unless a valid distribution for the 'step' data stream is specified.
#' 
#' @return If the simulated data are temporally regular (i.e., \code{lambda=NULL}) with no measurement error (i.e., \code{errorEllipse=NULL}), an object \code{\link{momentuHMMData}}, 
#' i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{...}{Data streams as specified by \code{dist}}
#' \item{x}{Either easting or longitude (if data streams include valid non-negative distribution for 'step')}
#' \item{y}{Either norting or latitude (if data streams include valid non-negative distribution for 'step')}
#' \item{...}{Covariates (if any)}
#' 
#' If simulated location data are temporally irregular (i.e., \code{lambda>0}) and/or include measurement error (i.e., \code{errorEllipse!=NULL}), a dataframe of:
#' \item{time}{Numeric time of each observed (and missing) observation}
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{x}{Either easting or longitude observed location}
#' \item{y}{Either norting or latitude observed location}
#' \item{...}{Data streams that are not derived from location (if applicable)}
#' \item{...}{Covariates at temporally-regular true (\code{mux},\code{muy}) locations (if any)}
#' \item{mux}{Either easting or longitude true location}
#' \item{muy}{Either norting or latitude true location}
#' \item{error_semimajor_axis}{error ellipse semi-major axis (if applicable)}
#' \item{error_semiminor_axis}{error ellipse semi-minor axis (if applicable)}
#' \item{error_ellipse_orientation}{error ellipse orientation (if applicable)}
#' \item{ln.sd.x}{log of the square root of the x-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{ln.sd.y}{log of the square root of the y-variance of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' \item{error.corr}{correlation term of bivariate normal error (if applicable; required for error ellipse models in \code{\link{crawlWrap}})}
#' 
#'
#' @details \itemize{
#' \item \code{simHierData} is very similar to \code{\link{simData}} except that instead of simply specifying the number of states (\code{nbStates}), distributions (\code{dist}), and a single t.p.m. formula (\code{formula}), the \code{hierStates} argument specifies the hierarchical nature of the states,
#' the \code{hierDist} argument specifies the hierarchical nature of the data streams, and the \code{hierFormula} argument specifies a t.p.m. formula for each level of the hierarchy.  All are specified as 
#' \code{\link[data.tree]{Node}} objects from the \code{\link[data.tree]{data.tree}} package.
#' 
#' \item \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#' 
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerLevel', then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' }
#' 
#' @seealso \code{\link{simData}}
#'                 
#' @references
#' 
#' Cornelissen, G. 2014. Cosinor-based rhythmometry. Theoretical Biology and Medical Modelling 11:16.
#' 
#' Leos-Barajas, V., Gangloff, E.J., Adam, T., Langrock, R., van Beest, F.M., Nabe-Nielsen, J. and Morales, J.M. 2017. 
#' Multi-scale modeling of animal movement and general behavior data using hidden Markov models with hierarchical structures. 
#' Journal of Agricultural, Biological and Environmental Statistics, 22 (3), 232-248.
#'
#' McClintock BT, London JM, Cameron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. 
#' Methods in Ecology and Evolution 6(3):266-277.
#' 
#' Rivest, LP, Duchesne, T, Nicosia, A, Fortin, D, 2016. A general angular regression model for the analysis of data on animal movement in ecology. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 65(3):445-463.
#'
#' @export
#' @importFrom stats rnorm runif rmultinom step terms.formula
#' @importFrom raster cellFromXY getValues
#' @importFrom CircStats rvm
#' @importFrom LaplacesDemon rbern
#' @importFrom Brobdingnag as.brob sum
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.tree Node Get

simHierData <- function(nbAnimals=1,hierStates,hierDist,
                    Par,beta=NULL,delta=NULL,
                    hierFormula=NULL,formulaDelta=~1,mixtures=1,formulaPi=NULL,
                    covs=NULL,nbCovs=0,
                    spatialCovs=NULL,
                    zeroInflation=NULL,
                    oneInflation=NULL,
                    circularAngleMean=NULL,
                    centers=NULL,
                    centroids=NULL,
                    angleCovs=NULL,
                    obsPerLevel,
                    initialPosition=c(0,0),
                    DM=NULL,cons=NULL,userBounds=NULL,workBounds=NULL,workcons=NULL,mvnCoords=NULL,
                    model=NULL,states=FALSE,
                    retrySims=0,
                    lambda=NULL,
                    errorEllipse=NULL)
{
  
  if(!is.null(cons) & is.null(model)) warning("cons argument is deprecated in momentuHMM >= 1.4.0. Please use workBounds instead.")
  if(!is.null(workcons) & is.null(model)) warning("workcons argument is deprecated in momentuHMM >= 1.4.0. Please use workBounds instead.")
  if(!is.null(workBounds) & (!is.null(cons) | !is.null(workcons)) & is.null(model)) stop("workBounds cannot be specified when using deprecated arguments cons or workcons; either workBounds or both cons and workcons must be NULL")
  
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    
    if(!inherits(model,"momentuHierHMM")) stop("model must be a 'momentuHierHMM' object")
    
    if(is.miHMM(model)){
      model <- model$miSum
    }
    
    model <- delta_bc(model)
    
    # extract simulation parameters from model
    nbStates <- length(model$stateNames)
    dist<-model$conditions$dist
    distnames<-names(dist)
    
    if(is.miSum(model)){
      model$mle <- lapply(model$Par$real,function(x) x$est)
      model$mle$beta <- model$Par$beta$beta$est
      model$mle$pi <- model$Par$real$pi$est
      model$mle$delta <- model$Par$real$delta$est
      model$mod <- list()
      if(!is.null(model$conditions$recharge)){
        nbRecovs <- ncol(model$g0covs) + ncol(model$reCovs)
        model$mle$g0 <- c(model$Par$beta$g0$est)
        names(model$mle$g0) <- colnames(model$Par$beta$g0$est)
        model$mle$theta <- c(model$Par$beta$theta$est)
        names(model$mle$theta) <- colnames(model$Par$beta$theta$est)
      } else nbRecovs <- 0
      model$mod$estimate <- expandPar(model$MIcombine$coefficients,model$conditions$optInd,unlist(model$conditions$fixPar),model$conditions$wparIndex,model$conditions$betaCons,nbStates,ncol(model$covsDelta)-1,model$conditions$stationary,nrow(model$Par$beta$beta$est)/model$conditions$mixtures-1,nbRecovs,model$conditions$mixtures,ncol(model$covsPi)-1)
      if(!is.null(model$mle$beta)) model$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(model$mle$beta),2,byrow=TRUE)
      if(!is.null(model$Par$beta$pi$est)) model$conditions$workBounds$pi<-matrix(c(-Inf,Inf),length(model$Par$beta$pi$est),2,byrow=TRUE)
      if(!is.null(model$Par$beta$delta$est)) model$conditions$workBounds$delta<-matrix(c(-Inf,Inf),length(model$Par$beta$delta$est),2,byrow=TRUE)
      if(!is.null(model$mle$g0)) model$conditions$workBounds$g0<-matrix(c(-Inf,Inf),length(model$mle$g0),2,byrow=TRUE)
      if(!is.null(model$mle$theta)) model$conditions$workBounds$theta<-matrix(c(-Inf,Inf),length(model$mle$theta),2,byrow=TRUE)
    }
    
    hierStates <- model$conditions$hierStates
    hierDist <- model$conditions$hierDist
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
          if(any(varform %in% factorcovs)){
            factorvar<-factorcovs %in% varform
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
    
    inputHierHMM <- formatHierHMM(data=NULL,hierStates,hierDist,hierFormula,formulaDelta,mixtures,workBounds,betaCons=NULL,fixPar=NULL)
    
    nbStates <- inputHierHMM$nbStates
    dist <- inputHierHMM$dist
    formula <- inputHierHMM$formula
    betaRef <- inputHierHMM$betaRef
    stateNames <- inputHierHMM$stateNames
    
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
      if(class(spatialCovs[[j]])!="RasterLayer") stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs")
    }
  } else nbSpatialCovs <- 0
  
  if(is.list(obsPerLevel)){
    if(length(obsPerLevel)!=nbAnimals) stop("obsPerLevel must be a list of length ",nbAnimals," where each element is a hierarchical Node")
    for(i in 1:length(obsPerLevel)){
      if(!inherits(obsPerLevel[[i]],"Node")) stop("obsPerLevel for individual ",i," must be a hierarchical Node")
      if(!("obs" %in% obsPerLevel[[i]]$fieldsAll)) stop("'obsPerLevel' for individual ",i," must include a 'obs' field")
      if(!all(obsPerLevel[[i]]$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",1:obsPerLevel[[i]]$count))) stop("'obsPerLevel' level types for individual ",i," can only include ",paste0(paste0("level",1:obsPerLevel[[i]]$count),collapse = ", "))
      for(j in obsPerLevel[[i]]$Get("name",filterFun=function(x) x$level==2)){
        if(length(which(obsPerLevel[[i]][[j]]$obs<1))>0)
          stop("obsPerLevel 'obs' field should have positive values.")
        if(length(obsPerLevel[[i]][[j]]$obs)==1)
          obsPerLevel[[i]][[j]]$obs <- rep(obsPerLevel[[i]][[j]]$obs,2)
        else if(length(obsPerLevel[[i]][[j]]$obs)!=2)
          stop("obsPerLevel 'obs' field should be of length 1 or 2.")
      }
    }
  } else {
    if(!inherits(obsPerLevel,"Node")) stop("obsPerLevel must be a hierarchical Node")
    if(!("obs" %in% obsPerLevel$fieldsAll)) stop("'obsPerLevel' must include a 'obs' field")
    if(!all(obsPerLevel$Get("name",filterFun=function(x) x$level==2) %in% paste0("level",1:obsPerLevel$count))) stop("'obsPerLevel' level types can only include ",paste0(paste0("level",1:obsPerLevel[[i]]$count),collapse = ", "))
    for(j in obsPerLevel$Get("name",filterFun=function(x) x$level==2)){
      if(length(which(obsPerLevel[[j]]$obs<1))>0)
        stop("obsPerLevel 'obs' field should have positive values.")
      if(length(obsPerLevel[[j]]$obs)==1)
        obsPerLevel[[j]]$obs <- rep(obsPerLevel[[j]]$obs,2)
      else if(length(obsPerLevel[[j]]$obs)!=2)
        stop("obsPerLevel 'obs' field should be of length 1 or 2.")
    }
    tmpObs<-obsPerLevel
    obsPerLevel<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      obsPerLevel[[i]]<-tmpObs
    }
  }
  
  if(is.list(initialPosition)){
    if(length(initialPosition)!=nbAnimals) stop("initialPosition must be a list of length ",nbAnimals)
    for(i in 1:nbAnimals){
      if(is.null(mvnCoords) || dist[[mvnCoords]] %in% c("mvnorm2","rw_mvnorm2")){
        if(length(initialPosition[[i]])!=2 | !is.numeric(initialPosition[[i]])) stop("each element of initialPosition must be a numeric vector of length 2")
      } else if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")){
        if(length(initialPosition[[i]])!=3 | !is.numeric(initialPosition[[i]])) stop("each element of initialPosition must be a numeric vector of length 3")
      }
    }
  } else {
    if(is.null(mvnCoords) || dist[[mvnCoords]] %in% c("mvnorm2","rw_mvnorm2")){
      if(length(initialPosition)!=2 | !is.numeric(initialPosition)) stop("initialPosition must be a numeric vector of length 2")
    } else if(!is.null(mvnCoords) && dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")){
      if(all(initialPosition==0)) initialPosition <- c(0,0,0)
      if(length(initialPosition)!=3 | !is.numeric(initialPosition)) stop("initialPosition must be a numeric vector of length 3")
    }
    tmpPos<-initialPosition
    initialPosition<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      initialPosition[[i]]<-tmpPos
    }
  }
  
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
        coordLevel <- hierDist$Get("parent",filterFun=isLeaf)$step$name
      }
  } else if("step" %in% distnames){
    if(inputs$dist[["step"]] %in% stepdists){
      data$x<-numeric()
      data$y<-numeric()
      coordLevel <- hierDist$Get("parent",filterFun=isLeaf)$step$name
    }    
  } else if(!is.null(mvnCoords)){
    data[[paste0(mvnCoords,".x")]]<-numeric()
    data[[paste0(mvnCoords,".y")]]<-numeric()
    if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) data[[paste0(mvnCoords,".z")]]<-numeric()
    coordLevel <- hierDist$Get("parent",filterFun=isLeaf)[[mvnCoords]]$name
  } else {
    if(nbSpatialCovs | length(centerInd) | length(centroidInd) | length(angleCovs)) stop("spatialCovs, angleCovs, centers, and/or centroids cannot be specified without valid step length and turning angle distributions")
    coordLevel <- NULL
  }
  
  if(!is.null(coordLevel)) distCoordLevel <- names(hierDist[[coordLevel]]$children)
  
  rwInd <- any(unlist(dist) %in% rwdists)
  
  #if(is.null(formula)) {
  #  if(allNbCovs) formula <- formula(paste0("~",paste0(c(colnames(allCovs),spatialcovnames),collapse="+")))
  #  else formula <- formula(~1)
  #}
  
  printMessage(nbStates,dist,p,DM,formula,formDelta,formPi,mixtures,"Simulating")
  
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
  
  newForm <- newFormulas(formula,nbStates)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  tmpCovs <- data.frame(ID=factor(1,levels=1:nbAnimals))
  if(!is.null(allCovs))
    tmpCovs <- cbind(tmpCovs,allCovs[1,,drop=FALSE])
  if(nbSpatialCovs){
    for(j in 1:nbSpatialCovs){
      getCell<-raster::cellFromXY(spatialCovs[[j]],initialPosition[[1]])
      if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
      tmpCovs[[spatialcovnames[j]]]<-spatialCovs[[j]][getCell]
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
  
  if(mixtures>1){
    if(!is.null(beta)){
      if(!is.list(beta)) stop("beta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
  }
  
  # build design matrix for recharge model
  if(!is.null(recharge)){
    g0covs <- model.matrix(recharge$g0,tmpCovs[1,,drop=FALSE])
    nbG0covs <- ncol(g0covs)-1
    recovs <- model.matrix(recharge$theta,tmpCovs)
    nbRecovs <- ncol(recovs)-1
    if(!nbRecovs) stop("invalid recharge model -- theta must include an intercept and at least 1 covariate")
    tmpcovs <- cbind(tmpCovs,rep(0,nrow(tmpCovs)))
    colnames(tmpcovs) <- c(colnames(tmpCovs),"recharge")
    tmpCovs <- tmpcovs
    if(!is.null(beta)){
      if(!is.list(beta)) stop("beta must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
    }
    for(j in 1:nbStates){
      formulaStates[[j]] <- as.formula(paste0(Reduce( paste, deparse(formulaStates[[j]]) ),"+recharge"))
    }
    formterms <- c(formterms,"recharge")
    newformula <- as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
    beta0 <- beta
  } else {
    nbRecovs <- 0
    nbG0covs <- 0
    g0covs <- NULL
    recovs <- NULL
    if(is.null(model) & !is.list(beta)){
      beta0 <- list(beta=beta)
    } else beta0 <- beta
  }
  
  tmpCovs$level <- factor(level[[1]][1],levels=lLevels[[1]])
  
  nbBetaCovs <- ncol(model.matrix(newformula,tmpCovs))
  
  inputHierHMM <- formatHierHMM(data=tmpCovs,hierStates,hierDist,hierFormula,formulaDelta,mixtures,workBounds,betaCons=NULL,fixPar=NULL,checkData=FALSE)
  
  if(is.null(beta0$beta)){
    beta0$beta <- matrix(rnorm(nbStates*(nbStates-1)*nbBetaCovs*mixtures)[inputHierHMM$betaCons],nrow=nbBetaCovs*mixtures)
    beta0$beta[which(!is.na(inputHierHMM$fixPar$beta))] <- inputHierHMM$fixPar$beta[which(!is.na(inputHierHMM$fixPar$beta))]
    workBounds <- inputHierHMM$workBounds
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
  
  if(is.null(workBounds)) {
    workBounds <- vector('list',length(distnames))
    names(workBounds) <- distnames
  }
  workBounds <- getWorkBounds(workBounds,distnames,unlist(Par[distnames]),parindex,parCount,inputs$DM,beta0,deltaB)
  
  wnbeta <- w2wn(beta0$beta,workBounds$beta)
  wnpi <- w2wn(beta0$pi,workBounds$pi)
  if(!is.null(recharge)){
    wng0 <- w2wn(beta0$g0,workBounds$g0)
    wntheta <- w2wn(beta0$theta,workBounds$theta)
  }
  
  mix <- rep(1,nbAnimals)
  
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
      
      if(!nbSpatialCovs & !length(centerInd) & !length(centroidInd) & !length(angleCovs) & !rwInd) {
        if(!is.null(recharge)){
          g0 <- model.matrix(recharge$g0,subCovs[1,,drop=FALSE]) %*% wng0
          subCovs[,"recharge"] <- cumsum(c(g0,model.matrix(recharge$theta,subCovs[-nrow(subCovs),])%*%wntheta))
        }
        DMcov <- model.matrix(newformula,subCovs)
        
        # format parameters
        DMinputs<-getDM(subCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
        
        ncmean <- get_ncmean(distnames,fullDM,inputs$circularAngleMean,nbStates)
        nc <- ncmean$nc
        meanind <- ncmean$meanind
        
        covsDelta <- model.matrix(formDelta,subCovs[1,,drop=FALSE])
        covsPi <- model.matrix(formPi,subCovs[1,,drop=FALSE])
        fullsubPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nbObs,inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
        
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
            subSpatialcovs[1,j]<-spatialCovs[[j]][getCell]
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
          g0 <- model.matrix(recharge$g0,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wng0
          subCovs[1,"recharge"] <- g0 #+ model.matrix(recharge$theta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wntheta
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
        
        covsPi <- model.matrix(formPi,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        pie <- mlogit(wnpi,covsPi,nbCovsPi,1,mixtures)
        
        # assign individual to mixture
        if(mixtures>1) mix[zoo] <- sample.int(mixtures,1,prob=pie)
        
        g <- model.matrix(newformula,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
        covsDelta <- model.matrix(formDelta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        
        delta0 <- mlogit(deltaB[(mix[zoo]-1)*(nbCovsDelta+1)+1:(nbCovsDelta+1),,drop=FALSE],covsDelta,nbCovsDelta,1,nbStates,mixtures)
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
      
        if(!is.null(coordLevel)) coordNA <- which(level[[zoo]]==gsub("level","",coordLevel))[sum(level[[zoo]]==gsub("level","",coordLevel))]
      
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
              DMinputs<-getDM(cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE]),inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
              fullDM <- DMinputs$fullDM
              DMind <- DMinputs$DMind
              wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
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
              subPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
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
                  subSpatialcovs[k+1,j]<-spatialCovs[[j]][getCell]
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
                subCovs[k+1,"recharge"] <- subCovs[k,"recharge"] + model.matrix(recharge$theta,cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE])) %*% wntheta
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
          }
      } else if("step" %in% distnames){
        if(inputs$dist[["step"]] %in% stepdists){
          d$x=c(0,cumsum(d$step)[-nrow(d)])
          d$y=X[,2]
        }    
      } else if(!is.null(mvnCoords)){
        d[[paste0(mvnCoords,".x")]] <- d[[mvnCoords]][,1]
        d[[paste0(mvnCoords,".y")]] <- d[[mvnCoords]][,2]
        if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) d[[paste0(mvnCoords,".z")]] <- d[[mvnCoords]][,3]
        d[[mvnCoords]] <- NULL
      }
      for(j in angleCovs[which(angleCovs %in% names(subCovs))])
        allCovs[cumNbObs[zoo]+1:nbObs,j] <- subCovs[,j]
      if(length(centerInd)) centerCovs[cumNbObs[zoo]+1:nbObs,] <- subCovs[,centerNames]
      if(length(centroidInd)) centroidCovs[cumNbObs[zoo]+1:nbObs,] <- subCovs[,centroidNames]
      data <- rbind(data,d)
    }
    
    if(nbSpatialCovs>0) colnames(allSpatialcovs)<-spatialcovnames
    
    if(nbCovs>0)
      data <- cbind(data,allCovs)
    
    if(nbSpatialCovs>0)
      data <- cbind(data,allSpatialcovs)
    
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
    
    if(!is.null(mvnCoords)){
      attr(data,'coords') <- paste0(mvnCoords,c(".x",".y"))
      #tmpNames <- colnames(data)[-which(colnames(data)=="ID")]
      #prepDat <- prepData(data,coordNames = paste0(mvnCoords,c(".x",".y")))
      #data$step <- prepDat$step
      #data$angle <- prepDat$angle
      #data$x <- prepDat$x
      #data$y <- prepDat$y
      #data <- data[,c("ID","step","angle",tmpNames,"x","y")]
    }
    
    # account for observation error (if any)
    out<-simObsHierData(momentuHierHMMData(data),lambda,errorEllipse,coordLevel)
    
    message("DONE")
    return(out)
  } else {
    simCount <- 0
    cat("Attempting to simulate tracks within spatial extent(s) of raster layers(s). Press 'esc' to force exit from 'simHierData'\n",sep="")
    while(simCount < retrySims){
      cat("\r    Attempt ",simCount+1," of ",retrySims,"...",sep="")
      tmp<-suppressMessages(tryCatch(simHierData(nbAnimals,hierStates,hierDist,
                                             Par,beta,delta,
                                             hierFormula,formulaDelta,mixtures,formulaPi,
                                             covs,nbCovs,
                                             spatialCovs,
                                             zeroInflation,
                                             oneInflation,
                                             circularAngleMean,
                                             centers,
                                             centroids,
                                             angleCovs,
                                             obsPerLevel,
                                             initialPosition,
                                             DM,cons,userBounds,workBounds,workcons,mvnCoords,
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
