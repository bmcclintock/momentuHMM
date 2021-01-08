
#' Simulation tool
#'
#' Simulates data from a (multivariate) hidden Markov model. Movement data are assumed to be in Cartesian coordinates (not longitude/latitude) and can be generated with or without observation error attributable to temporal irregularity or location measurement error.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'bern', 'beta', 'cat', 'exp', 'gamma', 'lnorm', 'logis', 'negbinom', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. For example,
#' \code{dist=list(step='gamma', angle='vm', dives='pois')} indicates 3 data streams ('step', 'angle', and 'dives')
#' and their respective probability distributions ('gamma', 'vm', and 'pois').
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
#' @param delta Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}. If \code{formulaDelta} includes a formula, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix, where k is the number of covariates and the columns correspond to states 2:\code{nbStates}. See details below.
#' @param formula Regression formula for the transition probability covariates. Default: \code{~1} (no covariate effect). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities.
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
#' \code{formula} and/or \code{DM}.  Note that \code{simData} usually takes longer to generate simulated data when \code{spatialCovs} is specified.
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
#' @param centroids List where each element is a data frame consisting of at least \code{max(unlist(obsPerAnimal))} rows that provides the x-coordinates ('x') and y-coordinates ('y) for centroids (i.e., dynamic activity centers where the coordinates can change for each time step)
#' from which distance and angle covariates will be calculated based on the simulated location data. These distance and angle 
#' covariates can be included in \code{formula} and \code{DM} using the names of \code{centroids}.  If no list names are provided, then generic names are generated 
#' for the distance and angle covariates (e.g., 'centroid1.dist', 'centroid1.angle', 'centroid2.dist', 'centroid2.angle'); otherwise the covariate names are derived from the list names
#' of \code{centroids} as \code{paste0(rep(names(centroids),each=2),c(".dist",".angle"))}. Note that the angle covariates for each centroid are calculated relative to 
#' the previous movement direction instead of standard directions relative to the x-axis; this is to allow turning angles to be simulated as a function of these covariates using circular-circular regression.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{covs} or \code{spatialCovs} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' using \code{\link{circAngles}}.
#' @param obsPerAnimal Either the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated for each animal are uniformously picked from this interval. Alternatively, \code{obsPerAnimal} can be specified as
#' a list of length \code{nbAnimals} with each element providing the number of observations (if single value) or the bounds (if vector of two values) for each individual.
#' Default: \code{c(500,1500)}.
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
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2,dimnames=list(c("mean_1",
#' "mean_2","concentration_1","concentration_2"))))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. For each matrix, the first column pertains to the lower bound and the second column the upper bound.
#' For data streams, each element of \code{workBounds} should be a k x 2 matrix with the same name of the corresponding element of 
#' \code{Par}, where k is the number of parameters. For transition probability parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``beta'', where k=\code{length(beta)}. For initial distribution parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``delta'', where k=\code{length(delta)}.
#' \code{workBounds} is ignored for any given data stream unless \code{DM} is also specified.
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the t.p.m. multinomial logit link. Default: NULL, in which case
#' the diagonal elements of the t.p.m. are the reference. See \code{\link{fitHMM}}.
#' @param mvnCoords Character string indicating the name of location data that are to be simulated using a multivariate normal distribution. For example, if \code{mu="rw_mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are intended to be location data, then \code{mvnCoords="mu"} needs to be specified in order for these data to be treated as such.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param model A \code{\link{momentuHMM}}, \code{\link{momentuHierHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerAnimal}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
#' \code{covs}, \code{spatialCovs}, \code{centers}, and \code{centroids}. It is not appropriate to simulate movement data from a \code{model} that was fitted to latitude/longitude data (because \code{simData} assumes Cartesian coordinates).
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
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' 
#' @return If the simulated data are temporally regular (i.e., \code{lambda=NULL}) with no measurement error (i.e., \code{errorEllipse=NULL}), an object \code{\link{momentuHMMData}} (or \code{\link{momentuHierHMMData}}), 
#' i.e., a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{...}{Data streams as specified by \code{dist} (or \code{hierDist})}
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
#' \item \code{simHierData} is very similar to \code{\link{simData}} except that instead of simply specifying the number of states (\code{nbStates}), distributions (\code{dist}), observations (\code{obsPerAnimal}), covariates (\code{nbCovs}), and a single t.p.m. formula (\code{formula}), the \code{hierStates} argument specifies the hierarchical nature of the states,
#' the \code{hierDist} argument specifies the hierarchical nature of the data streams, the \code{obsPerLevel} argument specifies the number of observations for each level of the hierarchy, the \code{nbHierCovs} argument specifies the number of covariates for each level of the hierarchy, and the \code{hierFormula} argument specifies a t.p.m. formula for each level of the hierarchy.
#' All of the hierarhcial arguments in \code{simHierData} are specified as \code{\link[data.tree]{Node}} objects from the \code{\link[data.tree]{data.tree}} package.

#' \item x- and y-coordinate location data are generated only if valid 'step' and 'angle' data streams are specified.  Vaild distributions for 'step' include 
#' 'gamma', 'weibull', 'exp', and 'lnorm'.  Valid distributions for 'angle' include 'vm' and 'wrpcauchy'.  If only a valid 'step' data stream is specified, then only x-coordinates
#' are generated.
#' 
#' \item If \code{DM} is specified for a particular data stream, then the initial values are specified on 
#' the working (i.e., beta) scale of the parameters. The working scale of each parameter is determined by the link function used.
#' The function \code{\link{getParDM}} is intended to help with obtaining initial values on the working scale when specifying a design matrix and other 
#' parameter constraints. 
#' 
#' \item Simulated data that are temporally regular (i.e., \code{lambda=NULL}) and without location measurement error (i.e., \code{errorEllipse=NULL}) are returned
#' as a \code{\link{momentuHMMData}} (or \code{\link{momentuHierHMMData}}) object suitable for analysis using \code{\link{fitHMM}}.
#' 
#' \item Simulated location data that are temporally-irregular (i.e., \code{lambda>0}) and/or with location measurement error (i.e., \code{errorEllipse!=NULL}) are returned
#' as a data frame suitable for analysis using \code{\link{crawlWrap}}.
#' 
#' \item The matrix \code{beta} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 \code{formula} covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix -
#' filled in row-wise).
#' In a covariate-free model (default), \code{beta} has one row, for the intercept.
#' 
#' \item State-specific formulas can be specified in \code{DM} using special formula functions. These special functions can take
#' the names \code{paste0("state",1:nbStates)} (where the integer indicates the state-specific formula).  For example, 
#' \code{DM=list(step=list(mean=~cov1+state1(cov2),sd=~cov2+state2(cov1)))} includes \code{cov1} on the mean parameter for all states, \code{cov2}
#' on the mean parameter for state 1, \code{cov2} on the sd parameter for all states, and \code{cov1} on the sd parameter for state 2.
#'
#' \item State- and parameter-specific formulas can be specified for transition probabilities in \code{formula} using special formula functions.
#' These special functions can take the names \code{paste0("state",1:nbStates)} (where the integer indicates the current state from which transitions occur),
#' \code{paste0("toState",1:nbStates)} (where the integer indicates the state to which transitions occur),
#' or \code{paste0("betaCol",nbStates*(nbStates-1))} (where the integer indicates the column of the \code{beta} matrix).  For example with \code{nbStates=3},
#' \code{formula=~cov1+betaCol1(cov2)+state3(cov3)+toState1(cov4)} includes \code{cov1} on all transition probability parameters, \code{cov2} on the \code{beta} column corresponding
#' to the transition from state 1->2, \code{cov3} on transition probabilities from state 3 (i.e., \code{beta} columns corresponding to state transitions 3->1 and 3->2),
#' and \code{cov4} on transition probabilities to state 1 (i.e., \code{beta} columns corresponding to state transitions 2->1 and 3->1).
#' 
#' \item Cyclical relationships (e.g., hourly, monthly) may be simulated using the \code{consinor(x,period)} special formula function for covariate \code{x}
#' and sine curve period of time length \code{period}. For example, if 
#' the data are hourly, a 24-hour cycle can be simulated using \code{~cosinor(cov1,24)}, where the covariate \code{cov1} is a repeating series
#' of integers \code{0,1,...,23,0,1,...,23,0,1,...} (note that \code{simData} will not do this for you, the appropriate covariate must be specified using the \code{covs} argument; see example below). 
#' The \code{cosinor(x,period)} function converts \code{x} to 2 covariates
#' \code{cosinorCos(x)=cos(2*pi*x/period)} and \code{consinorSin(x)=sin(2*pi*x/period} for inclusion in the model (i.e., 2 additional parameters per state). The amplitude of the sine wave
#' is thus \code{sqrt(B_cos^2 + B_sin^2)}, where \code{B_cos} and \code{B_sin} are the working parameters correponding to \code{cosinorCos(x)} and \code{cosinorSin(x)}, respectively (e.g., see Cornelissen 2014).
#'
#' When the circular-circular regression model is used, the special function \code{angleFormula(cov,strength,by)} can be used in \code{DM} for the mean of angular distributions (i.e. 'vm', 'vmConsensus', and 'wrpcauchy'),
#' where \code{cov} is an angle covariate (e.g. wind direction), \code{strength} is a positive real covariate (e.g. wind speed), and \code{by} is an optional factor variable for individual- or group-level effects (e.g. ID, sex). This allows angle covariates to be weighted based on their strength or importance at time step t as in
#' Rivest et al. (2016).
#' }
#' 
#' @seealso \code{\link{prepData}}, \code{\link{simObsData}}
#'
#' @examples
#' # 1. Pass a fitted model to simulate from
#' # (m is a momentuHMM object - as returned by fitHMM - automatically loaded with the package)
#' # We keep the default nbAnimals=1.
#' m <- example$m
#' obsPerAnimal=c(50,100)
#' data <- simData(model=m,obsPerAnimal=obsPerAnimal)
#'
#' \dontrun{
#' # 2. Pass the parameters of the model to simulate from
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean_1, mean_2, sd_1, sd_2, zeromass_1, zeromass_2
#' anglePar <- c(pi,0,0.5,2) # mean_1, mean_2, concentration_1, concentration_2
#' omegaPar <- c(1,10,10,1) # shape1_1, shape1_2, shape2_1, shape2_2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' omegaDist <- "beta"
#' data <- simData(nbAnimals=4,nbStates=2,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),
#'                 Par=list(step=stepPar,angle=anglePar,omega=omegaPar),nbCovs=2,
#'                 zeroInflation=list(step=TRUE),
#'                 obsPerAnimal=obsPerAnimal)
#'
#' # 3. Include covariates
#' # (note that it is useless to specify "nbCovs", which are overruled
#' # by the number of columns of "cov")
#' cov <- data.frame(temp=log(rnorm(500,20,5)))
#' stepPar <- c(log(10),0.1,log(100),-0.1,log(5),log(25)) # working scale parameters for step DM
#' anglePar <- c(pi,0,0.5,2) # mean_1, mean_2, concentration_1, concentration_2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=2,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                 Par=list(step=stepPar,angle=anglePar),
#'                 DM=list(step=list(mean=~temp,sd=~1)),
#'                 covs=cov,
#'                 obsPerAnimal=obsPerAnimal)
#'                 
#' # 4. Include example 'forest' spatial covariate raster layer
#' # nbAnimals and obsPerAnimal kept small to reduce example run time
#' spatialCov<-list(forest=forest)
#' data <- simData(nbAnimals=1,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                 Par=list(step=c(100,1000,50,100),angle=c(0,0,0.1,5)),
#'                 beta=matrix(c(5,-10,-25,50),nrow=2,ncol=2,byrow=TRUE),
#'                 formula=~forest,spatialCovs=spatialCov,
#'                 obsPerAnimal=250,states=TRUE,
#'                 retrySims=100)
#'                 
#' # 5. Specify design matrix for 'omega' data stream
#' # natural scale parameters for step and angle
#' stepPar <- c(1,10,1,5) # shape_1, shape_2, scale_1, scale_2
#' anglePar <- c(pi,0,0.5,0.7) # mean_1, mean_2, concentration_1, concentration_2
#' 
#' # working scale parameters for omega DM
#' omegaPar <- c(log(1),0.1,log(10),-0.1,log(10),-0.1,log(1),0.1)
#' 
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' omegaDist <- "beta"
#' 
#' data <- simData(nbStates=2,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),
#'                 Par=list(step=stepPar,angle=anglePar,omega=omegaPar),nbCovs=2,
#'                 DM=list(omega=list(shape1=~cov1,shape2=~cov2)),
#'                 obsPerAnimal=obsPerAnimal,states=TRUE)
#'                 
#' # 6. Include temporal irregularity and location measurement error
#' lambda <- 2 # expect 2 observations per time step
#' errorEllipse <- list(M=50,m=25,r=180)
#' obsData <- simData(model=m,obsPerAnimal=obsPerAnimal,
#'                    lambda=lambda, errorEllipse=errorEllipse)
#'                    
#' # 7. Cosinor and state-dependent formulas
#' nbStates<-2
#' dist<-list(step="gamma")
#' Par<-list(step=c(100,1000,50,100))
#' 
#' # include 24-hour cycle on all transition probabilities
#' # include 12-hour cycle on transitions from state 2
#' formula=~cosinor(hour24,24)+state2(cosinor(hour12,12))
#' 
#' # specify appropriate covariates
#' covs<-data.frame(hour24=0:23,hour12=0:11)
#' 
#' beta<-matrix(c(-1.5,1,1,NA,NA,-1.5,-1,-1,1,1),5,2)
#' # row names for beta not required but can be helpful
#' rownames(beta)<-c("(Intercept)",
#'                   "cosinorCos(hour24, 24)",
#'                   "cosinorSin(hour24, 24)",
#'                   "cosinorCos(hour12, 12)",
#'                   "cosinorSin(hour12, 12)")
#' data.cos<-simData(nbStates=nbStates,dist=dist,Par=Par,
#'                   beta=beta,formula=formula,covs=covs)     
#'                   
#' # 8. Piecewise constant B-spline on step length mean and angle concentration
#' library(splines2)
#' nObs <- 1000 # length of simulated track
#' cov <- data.frame(time=1:nObs) # time covariate for splines
#' dist <- list(step="gamma",angle="vm")
#' stepDM <- list(mean=~bSpline(time,df=2,degree=0),sd=~1)
#' angleDM <- list(mean=~1,concentration=~bSpline(time,df=2,degree=0))
#' DM <- list(step=stepDM,angle=angleDM)
#' Par <- list(step=c(log(1000),1,-1,log(100)),angle=c(0,log(10),2,-5))
#'
#' data.spline<-simData(obsPerAnimal=nObs,nbStates=1,dist=dist,Par=Par,DM=DM,covs=cov)        
#' 
#' # 9. Initial state (delta) based on covariate
#' nObs <- 100
#' dist <- list(step="gamma",angle="vm")
#' Par <- list(step=c(100,1000,50,100),angle=c(0,0,0.01,0.75))
#' 
#' # create sex covariate
#' cov <- data.frame(sex=factor(rep(c("F","M"),each=nObs))) # sex covariate
#' formulaDelta <- ~ sex + 0
#' 
#' # Female begins in state 1, male begins in state 2
#' delta <- matrix(c(-100,100),2,1,dimnames=list(c("sexF","sexM"),"state 2")) 
#'
#' data.delta<-simData(nbAnimals=2,obsPerAnimal=nObs,nbStates=2,dist=dist,Par=Par,
#'                     delta=delta,formulaDelta=formulaDelta,covs=cov,
#'                     beta=matrix(-1.5,1,2),states=TRUE)        
#' }                
#' @references
#' 
#' Cornelissen, G. 2014. Cosinor-based rhythmometry. Theoretical Biology and Medical Modelling 11:16.
#'
#' McClintock BT, London JM, Cameron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. 
#' Methods in Ecology and Evolution 6(3):266-277.
#' 
#' Rivest, LP, Duchesne, T, Nicosia, A, Fortin, D, 2016. A general angular regression model for the analysis of data on animal movement in ecology. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 65(3):445-463.
#'
#' @export
#' @importFrom stats rnorm runif rmultinom step terms.formula
# @importFrom raster cellFromXY getValues is.factor levels
# @importFrom moveHMM simData
#' @importFrom CircStats rvm
#' @importFrom Brobdingnag as.brob sum
#' @importFrom mvtnorm rmvnorm
#' @importFrom extraDistr rcat

simData <- function(nbAnimals=1,nbStates=2,dist,
                    Par,beta=NULL,delta=NULL,
                    formula=~1,formulaDelta=NULL,mixtures=1,formulaPi=NULL,
                    covs=NULL,nbCovs=0,
                    spatialCovs=NULL,
                    zeroInflation=NULL,
                    oneInflation=NULL,
                    circularAngleMean=NULL,
                    centers=NULL,
                    centroids=NULL,
                    angleCovs=NULL,
                    obsPerAnimal=c(500,1500),
                    initialPosition=c(0,0),
                    DM=NULL,userBounds=NULL,workBounds=NULL,betaRef=NULL,mvnCoords=NULL,stateNames=NULL,
                    model=NULL,states=FALSE,
                    retrySims=0,
                    lambda=NULL,
                    errorEllipse=NULL,
                    ncores=1)
{
  
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    
    if(is.miHMM(model)){
      model <- model$miSum
    } 
    if(inherits(model,"momentuHierHMM") | inherits(model,"hierarchical")) stop("model can not be a 'momentuHierHMM' or 'hierarchical' object; use simHierData instead")
    
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
    } else {
      if(!inherits(model,"momentuHMM")) stop("model must be a 'momentuHMM' object")
    }
    
    userBounds <- model$conditions$bounds
    workBounds <- model$conditions$workBounds
    mvnCoords <- model$conditions$mvnCoords
    stateNames<-model$stateNames
    estAngleMean<-model$conditions$estAngleMean
    circularAngleMean<-model$conditions$circularAngleMean
    DM <- model$conditions$DM
    betaRef <- model$conditions$betaRef
    zeroInflation <- model$conditions$zeroInflation
    oneInflation <- model$conditions$oneInflation
    formula <- model$conditions$formula
    if(is.null(model$condition$formulaDelta)){
      formulaDelta <- formDelta <- ~1
    } else formulaDelta <- formDelta <- model$condition$formulaDelta
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
      }
    }
    for(i in distnames[which(dist %in% angledists)]){
      if(!estAngleMean[[i]]){
        estAngleMean[[i]]<-TRUE
        userBounds[[i]]<-rbind(matrix(rep(c(-pi,pi),nbStates),nbStates,2,byrow=TRUE),userBounds[[i]])
        workBounds[[i]]<-rbind(matrix(rep(c(-Inf,Inf),nbStates),nbStates,2,byrow=TRUE),workBounds[[i]])
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
    nbCovsDelta <- ncol(model$covsDelta)-1
    if(nbStates>1){
      beta <- nw2w(model$mle$beta,model$conditions$workBounds$beta)
      foo <- length(model$mod$estimate)-length(g0)-length(theta)-(nbCovsDelta+1)*(nbStates-1)*mixtures+1
      delta <- matrix(model$mod$estimate[foo:(length(model$mod$estimate)-length(g0)-length(theta))],nrow=(nbCovsDelta+1)*mixtures) 
    } else {
      formulaDelta <- NULL
    }
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
          form<-stats::formula(paste("~",cov))
          varform<-all.vars(form)
          if(any(varform %in% factorcovs) && !all(varform %in% factorterms)){
            factorvar<-factorcovs %in% (varform[!(varform %in% factorterms)])
            covsCol[jj]<-rep(factorterms,times=unlist(lapply(model$data[factorterms],nlevels)))[which(factorvar)]
          } 
        }
      }
      covsCol<-unique(covsCol)
      covsCol <- covsCol[!(covsCol %in% "ID")]
      
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
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
    distnames<-names(dist)
    if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
    Par <- Par[distnames]
    
    if(is.null(formulaPi)){
      formPi <- ~1
    } else formPi <- formulaPi
    if(is.null(formulaDelta)){
      formDelta <- ~1
    } else formDelta <- formulaDelta
    mHind <- (requireNamespace("moveHMM", quietly = TRUE) && is.null(DM) & is.null(userBounds) & is.null(workBounds) & is.null(spatialCovs) & is.null(centers) & is.null(centroids) & ("step" %in% names(dist)) & all(unlist(initialPosition)==c(0,0)) & is.null(lambda) & is.null(errorEllipse) & !is.list(obsPerAnimal) & is.null(covs) & !nbCovs & !length(attr(stats::terms.formula(formula),"term.labels")) & !length(attr(stats::terms.formula(formDelta),"term.labels")) & is.null(delta) & is.null(betaRef) & is.null(mvnCoords) & mixtures==1) # indicator for moveHMM::simData

    if(all(names(dist) %in% c("step","angle")) & all(unlist(dist) %in% moveHMMdists) & mHind){
      zi <- FALSE
      if(!is.null(zeroInflation$step)) zi <- zeroInflation$step
      if(is.null(dist$angle)) dist$angle<-"none"
      data <- moveHMM::simData(nbAnimals, nbStates, dist$step, dist$angle, Par$step, Par$angle, beta, covs, nbCovs, zi, obsPerAnimal, model, states)
      attr(data,"class") <- "data.frame"
      data$ID <- as.factor(data$ID)
      return(momentuHMMData(data))
    }
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
    if(sum(unlist(dist) %in% rwdists)>1) stop("sorry, simData currently only supports a single multivariate normal random walk distribution (and it must correspond to location data)")
  } else if(any(unlist(dist) %in% rwdists)) stop("sorry, simData currently only supports random walk distributions for multivariate location data identified through the mvnCoords argument")
  if(any(unlist(dist)=="rw_norm")) stop("sorry, 'rw_norm' is not currently supported by simData")
  
  estAngleMean <- vector('list',length(distnames))
  names(estAngleMean) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% angledists) estAngleMean[[i]]<-TRUE
    else estAngleMean[[i]]<-FALSE
  }
  
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,stateNames,checkInflation = TRUE)
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
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package \"raster\" needed for spatial covariates. Please install it.",
           call. = FALSE)
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

  if(is.list(obsPerAnimal)){
    if(length(obsPerAnimal)!=nbAnimals) stop("obsPerAnimal must be a list of length ",nbAnimals)
    for(i in 1:length(obsPerAnimal)){
      if(length(which(obsPerAnimal[[i]]<1))>0)
        stop("obsPerAnimal elements should have positive values.")
      if(length(obsPerAnimal[[i]])==1)
        obsPerAnimal[[i]] <- rep(obsPerAnimal[[i]],2)
      else if(length(obsPerAnimal[[i]])!=2)
        stop("obsPerAnimal elements should be of length 1 or 2.")
    }
  } else {
    if(length(which(obsPerAnimal<1))>0) 
      stop("obsPerAnimal should have positive values.")
    if(length(obsPerAnimal)==1)
      obsPerAnimal <- rep(obsPerAnimal,2)
    else if(length(obsPerAnimal)!=2)
      stop("obsPerAnimal should be of length 1 or 2.")
    tmpObs<-obsPerAnimal
    obsPerAnimal<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      obsPerAnimal[[i]]<-tmpObs
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
  
  if(!is.null(betaRef)){
    if(length(betaRef)!=nbStates) stop("betaRef must be a vector of length ",nbStates)
    if(!is.numeric(betaRef)) stop("betaRef must be a numeric vector")
    if(min(betaRef)<1 | max(betaRef)>nbStates) stop("betaRef elements must be between 1 and ",nbStates)
  } else {
    betaRef <- 1:nbStates
  }
  betaRef <- as.integer(betaRef)

  #######################################
  ## Prepare parameters for simulation ##
  #######################################
  # define number of observations for each animal
  allNbObs <- rep(NA,nbAnimals)
  for(zoo in 1:nbAnimals) {
    if(obsPerAnimal[[zoo]][1]!=obsPerAnimal[[zoo]][2])
      allNbObs[zoo] <- sample(obsPerAnimal[[zoo]][1]:obsPerAnimal[[zoo]][2],size=1)
    else
      allNbObs[zoo] <- obsPerAnimal[[zoo]][1]
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
  } else centerNames <- centerCovs <- NULL
  
  centroidInd<-NULL
  if(!is.null(centroids)){
    if(!is.list(centroids)) stop("centroids must be a named list")
    centroidNames <- character()
    for(j in 1:length(centroids)){
      if(!is.data.frame(centroids[[j]])) stop("each element of centroids must be a data frame")
      if(dim(centroids[[j]])[1]<max(unlist(obsPerAnimal)) | dim(centroids[[j]])[2]!=2) stop("each element of centroids must be a data frame consisting of at least ",max(unlist(obsPerAnimal))," rows (i.e., the maximum number of observations per animal) and 2 columns (i.e., x- and y-coordinates)")
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
  } else centroidNames <- centroidCovs <- NULL
  
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

  # initial state distribution
  if(is.null(delta)) delta0 <- matrix(rep(1,nbStates)/nbStates,mixtures,nbStates)
  else delta0 <- delta
  
  zeroMass<-oneMass<-vector('list',length(inputs$dist))
  names(zeroMass)<-names(oneMass)<-distnames

  allStates <- NULL
  allSpatialcovs<-NULL
  
  #make sure 'step' preceeds 'angle'
  if(all(c("step","angle") %in% distnames)){
    distnames<-c("step","angle",distnames[!(distnames %in% c("step","angle"))])
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
        data$x<-numeric()
        data$y<-numeric()
      }
  } else if("step" %in% distnames){
      if(inputs$dist[["step"]] %in% stepdists){
        data$x<-numeric()
        data$y<-numeric()
      }    
  } else if(!is.null(mvnCoords)){
    data[[paste0(mvnCoords,".x")]]<-numeric()
    data[[paste0(mvnCoords,".y")]]<-numeric()
    if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) data[[paste0(mvnCoords,".z")]]<-numeric()
  } else if(nbSpatialCovs | length(centerInd) | length(centroidInd) | length(angleCovs)) stop("spatialCovs, angleCovs, centers, and/or centroids cannot be specified without valid step length and turning angle distributions")
  
  rwInd <- any(unlist(dist) %in% rwdists)
  
  #if(is.null(formula)) {
  #  if(allNbCovs) formula <- stats::formula(paste0("~",paste0(c(colnames(allCovs),spatialcovnames),collapse="+")))
  #  else formula <- stats::formula(~1)
  #}
  
  printMessage(nbStates,dist,p,DM,formula,formDelta,formPi,mixtures,"Simulating",FALSE)
  
  if(length(all.vars(formula)))
    if(!all(all.vars(formula) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formula' covariate(s) not found")
  if(length(all.vars(formPi)))
    if(!all(all.vars(formPi) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formulaPi' covariate(s) not found")
  if(length(all.vars(formDelta)))
    if(!all(all.vars(formDelta) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formulaDelta' covariate(s) not found")
  if(("ID" %in% all.vars(formula) | "ID" %in% all.vars(formPi) | "ID" %in% all.vars(formDelta)) & nbAnimals<2) stop("ID cannot be a covariate when nbAnimals=1")
  
  newForm <- newFormulas(formula,nbStates,betaRef)
  formulaStates <- newForm$formulaStates
  formterms <- newForm$formterms
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
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
  
  if(mixtures>1){
    if(!is.null(beta)){
      if(!is.list(beta)) stop("beta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
  }
  
  # build design matrix for recharge model
  if(!is.null(recharge)){
    g0covs <- stats::model.matrix(recharge$g0,tmpCovs[1,,drop=FALSE])
    nbG0covs <- ncol(g0covs)-1
    recovs <- stats::model.matrix(recharge$theta,tmpCovs)
    nbRecovs <- ncol(recovs)-1
    if(!nbRecovs) stop("invalid recharge model -- theta must include an intercept and at least 1 covariate")
    tmpcovs <- cbind(tmpCovs,rep(0,nrow(tmpCovs)))
    colnames(tmpcovs) <- c(colnames(tmpCovs),"recharge")
    tmpCovs <- tmpcovs
    if(!is.null(beta)){
      if(!is.list(beta)) stop("beta must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
    }
    for(j in 1:nbStates){
      formulaStates[[j]] <- stats::as.formula(paste0(Reduce( paste, deparse(formulaStates[[j]]) ),"+recharge"))
    }
    formterms <- c(formterms,"recharge")
    newformula <- stats::as.formula(paste0(Reduce( paste, deparse(newformula) ),"+recharge"))
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
  
  nbBetaCovs <- ncol(stats::model.matrix(newformula,tmpCovs))

  if(is.null(beta0$beta))
    beta0$beta <- matrix(rnorm(nbStates*(nbStates-1)*nbBetaCovs*mixtures),nrow=nbBetaCovs*mixtures)
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
      noBeta<-which(match(colnames(stats::model.matrix(newformula,tmpCovs)),colnames(stats::model.matrix(formulaStates[[state]],tmpCovs)),nomatch=0)==0)
      if(length(noBeta)) beta0$beta[noBeta,state] <- 0
    }
  }
  
  covsPi <- stats::model.matrix(formPi,tmpCovs)
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
  
  covsDelta <- stats::model.matrix(formDelta,tmpCovs)
  nbCovsDelta <- ncol(covsDelta)-1
  if(!nbCovsDelta & is.null(formulaDelta)){
    if(mixtures==1){
      if(length(delta0) != (nbCovsDelta+1)*nbStates)
        stop(paste("delta has the wrong length: it should have",nbStates*mixtures,"elements."))
      deltaB <- matrix(log(delta0[-1]/delta0[1]),1)
    } else {
      if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates | nrow(delta0)!=mixtures))
        stop(paste("delta has wrong dimensions: it should have",mixtures,"rows and",
                   nbStates,"columns."))
      deltaB <- apply(delta0,1,function(x) log(x[-1]/x[1]))
    }

  } else {
    if(is.null(dim(delta0)) || (ncol(delta0)!=nbStates-1 | nrow(delta0)!=(nbCovsDelta+1)*mixtures))
      stop(paste("delta has wrong dimensions: it should have",(nbCovsDelta+1)*mixtures,"rows and",
                 nbStates-1,"columns."))
    deltaB <- delta0
  }
  
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
    registerDoParallel(cores=ncores)
    withCallingHandlers(simDat <- foreach(zoo=1:nbAnimals,.combine='comb') %dorng% {
      
      message("\r        Simulating individual ",zoo,"... ",sep="")
      
      # number of observations for animal zoo
      nbObs <- allNbObs[zoo]
      d <- data.frame(ID=factor(rep(zoo,nbObs)))
      
      ###############################
      ## Simulate covariate values ##
      ###############################
      subCovs<-data.frame(ID=rep(factor(zoo,levels=1:nbAnimals),nbObs))
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
  
      X <- matrix(0,nrow=nbObs,ncol=2)
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
          g0 <- stats::model.matrix(recharge$g0,subCovs[1,,drop=FALSE]) %*% wng0
          subCovs[,"recharge"] <- cumsum(c(g0,stats::model.matrix(recharge$theta,subCovs[-nrow(subCovs),])%*%wntheta))
        }
        DMcov <- stats::model.matrix(newformula,subCovs)
        
        # format parameters
        DMinputs<-getDM(subCovs,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,zeroInflation,oneInflation,inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,p$Bndind,inputs$dist)
        if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")
        
        ncmean <- get_ncmean(distnames,fullDM,inputs$circularAngleMean,nbStates)
        nc <- ncmean$nc
        meanind <- ncmean$meanind

        covsDelta <- stats::model.matrix(formDelta,subCovs[1,,drop=FALSE])
        covsPi <- stats::model.matrix(formPi,subCovs[1,,drop=FALSE])
        fullsubPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,fullDM,DMind,nbObs,inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
        
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
          g0 <- stats::model.matrix(recharge$g0,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wng0
          subCovs[1,"recharge"] <- g0 #+ stats::model.matrix(recharge$theta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wntheta
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
        maxlag <- getDM(cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]),inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,zeroInflation,oneInflation,inputs$circularAngleMean,wlag=TRUE)$lag
        
        covsPi <- stats::model.matrix(formPi,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        pie <- mlogit(wnpi,covsPi,nbCovsPi,1,mixtures)
        
        # assign individual to mixture
        if(mixtures>1) mix[zoo] <- sample.int(mixtures,1,prob=pie)
        
        g <- stats::model.matrix(newformula,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
        covsDelta <- stats::model.matrix(formDelta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        
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
      
      for (k in 1:(nbObs-1)){
        
        if(nbSpatialCovs |  length(centerInd) | length(centroidInd) | length(angleCovs) | rwInd){
          # format parameters
          DMinputs<-getDM(cbind(subCovs[k-maxlag:0,,drop=FALSE],subSpatialcovs[k-maxlag:0,,drop=FALSE]),inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par,zeroInflation,oneInflation,inputs$circularAngleMean,wlag=TRUE)
          fullDM <- DMinputs$fullDM
          DMind <- DMinputs$DMind
          wpar <- n2w(Par,bounds,beta0,deltaB,nbStates,inputs$estAngleMean,inputs$DM,p$Bndind,inputs$dist)
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
          subPar <- w2n(wpar,bounds,parSize,nbStates,nbBetaCovs-1,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary=FALSE,fullDM,DMind,1,inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
        } else {
          subPar <- lapply(fullsubPar[distnames],function(x) x[,k,drop=FALSE])#fullsubPar[,k,drop=FALSE]
        }
        
        for(i in distnames){
          
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
                if(genData$step[k]>0){
                  phi <- phi + genData[[i]][k]
                } #else if(genData$step[k]==0) {
                  #genData[[i]][k] <- NA # angle = NA if step = 0
                #}
                m <- genData$step[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
                X[k+1,] <- X[k,] + m
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
            X[k+1,] <- genData[[i]][k,]
            d[[i]] <- X
          } else d[[i]] <- genData[[i]]
          
        }
        # get next state
        gamma <- matrix(0,nbStates,nbStates)
        gamma[cbind(1:nbStates,betaRef)] <- 1
        gamma <- t(gamma)
        
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
                subSpatialcovs[k+1,j] <- circAngles(subAnglecovs[k:(k+1),spatialcovnames[j]],data.frame(x=X[k:(k+1),1],y=X[k:(k+1),2]))[2] 
              }
            }
          }
          
          for(j in angleCovs[which(angleCovs %in% names(subCovs))]){
            subAnglecovs[k+1,j] <- subCovs[k+1,j]
            subCovs[k+1,j] <- circAngles(subAnglecovs[k:(k+1),j],data.frame(x=X[k:(k+1),1],y=X[k:(k+1),2]))[2] 
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
            subCovs[k+1,"recharge"] <- subCovs[k,"recharge"] + stats::model.matrix(recharge$theta,cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE])) %*% wntheta
          }
          for(i in distnames){
            if(dist[[i]] %in% rwdists){
              if(dist[[i]] %in% c("rw_mvnorm2")) subCovs[k+1,paste0(i,c(".x_tm1",".y_tm1"))] <- X[k+1,]
              else if(dist[[i]] %in% c("rw_mvnorm3")) subCovs[k+1,paste0(i,c(".x_tm1",".y_tm1",".z_tm1"))] <- X[k+1,]
            }
          }
          g <- stats::model.matrix(newformula,cbind(subCovs[k+1,,drop=FALSE],subSpatialcovs[k+1,,drop=FALSE])) %*% wnbeta[(mix[zoo]-1)*nbBetaCovs+1:nbBetaCovs,]
        } else {
          g <- gFull[k+1,,drop=FALSE]
        }
        gamma[!gamma] <- exp(g)
        gamma <- t(gamma)
        gamma <- gamma/apply(gamma,1,sum)
        Z[k+1] <- sample(1:nbStates,size=1,prob=gamma[Z[k],])  
      }
      #allStates <- c(allStates,Z)
      #if(nbSpatialCovs>0) {
        #allSpatialcovs <- rbind(allSpatialcovs,subSpatialcovs)
      #}
      
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
          d$x=c(initialPosition[[zoo]][1],initialPosition[[zoo]][1]+cumsum(d$step)[-nrow(d)])
          d$y=rep(initialPosition[[zoo]][2],nrow(d))
        }    
      } else if(!is.null(mvnCoords)){
        d[[paste0(mvnCoords,".x")]] <- d[[mvnCoords]][,1]
        d[[paste0(mvnCoords,".y")]] <- d[[mvnCoords]][,2]
        if(dist[[mvnCoords]] %in% c("mvnorm3","rw_mvnorm3")) d[[paste0(mvnCoords,".z")]] <- d[[mvnCoords]][,3]
        d[[mvnCoords]] <- NULL
      }
      for(j in angleCovs[which(angleCovs %in% names(subCovs))])
        allCovs[cumNbObs[zoo]+1:nbObs,j] <- subCovs[,j]
      if(length(centerInd)) centerCovs <- subCovs[,centerNames]
      if(length(centroidInd)) centroidCovs <- subCovs[,centroidNames]
      #data <- rbind(data,d)
      return(list(data=d,allCovs=allCovs[cumNbObs[zoo]+1:nbObs,,drop=FALSE],allSpatialcovs=subSpatialcovs,centerCovs=centerCovs,centroidCovs=centroidCovs,allStates=matrix(Z,ncol=1)))
    }
    ,warning=muffleRNGwarning)
    stopImplicitCluster()
    
    if(nbCovs>0)
      simDat$data <- cbind(simDat$data,simDat$allCovs)
    
    if(nbSpatialCovs>0){
      colnames(simDat$allSpatialcovs)<-spatialcovnames
      for(j in spatialcovnames){
        if(any(raster::is.factor(spatialCovs[[j]]))){
          simDat$allSpatialcovs[[j]] <- factor(simDat$allSpatialcovs[[j]],levels=unique(unlist(raster::levels(spatialCovs[[j]]))))
        }
      }
      simDat$data <- cbind(simDat$data,simDat$allSpatialcovs)
    }
    
    if(length(centerInd)){
      simDat$data <- cbind(simDat$data,simDat$centerCovs)
      for(j in which(grepl(".angle",names(simDat$data)))){
        if(names(simDat$data[j]) %in% centerNames)
          class(simDat$data[[j]]) <- c(class(simDat$data[[j]]), "angle")
      }
    }
    
    if(length(centroidInd)){
      simDat$data <- cbind(simDat$data,simDat$centroidCovs)
      for(j in which(grepl(".angle",names(simDat$data)))){
        if(names(simDat$data[j]) %in% centroidNames)
          class(simDat$data[[j]]) <- c(class(simDat$data[[j]]), "angle")
      }
    }
    
    # include states sequence in the data
    if(states)
      simDat$data <- cbind(simDat$data,states=simDat$allStates)
    
    for(i in distnames){
      if(inputs$dist[[i]] %in% angledists)
        class(simDat$data[[i]]) <- c(class(simDat$data[[i]]), "angle")
    }
    
    for(i in angleCovs){
      class(simDat$data[[i]]) <- c(class(simDat$data[[i]]), "angle")
    }
    
    if(!is.null(mvnCoords)){
      attr(simDat$data,'coords') <- paste0(mvnCoords,c(".x",".y"))
      #tmpNames <- colnames(data)[-which(colnames(data)=="ID")]
      #prepDat <- prepData(data,coordNames = paste0(mvnCoords,c(".x",".y")))
      #data$step <- prepDat$step
      #data$angle <- prepDat$angle
      #data$x <- prepDat$x
      #data$y <- prepDat$y
      #data <- data[,c("ID","step","angle",tmpNames,"x","y")]
    }
    
    # account for observation error (if any)
    out<-simObsData(momentuHMMData(simDat$data),lambda,errorEllipse)
    
    message("DONE")
    return(out)
  } else {
    simCount <- 0
    cat("Attempting to simulate tracks within spatial extent(s) of raster layers(s). Press 'esc' to force exit from 'simData'\n",sep="")
    while(simCount < retrySims){
      if(ncores==1) cat("    Attempt ",simCount+1," of ",retrySims,"...\n",sep="")
      else cat("\r    Attempt ",simCount+1," of ",retrySims,"...",sep="")
      tmp<-suppressMessages(tryCatch(simData(nbAnimals,nbStates,dist,
                          Par,beta,delta,
                          formula,formulaDelta,mixtures,formulaPi,
                          covs,nbCovs,
                          spatialCovs,
                          zeroInflation,
                          oneInflation,
                          circularAngleMean,
                          centers,
                          centroids,
                          angleCovs,
                          obsPerAnimal,
                          initialPosition,
                          DM,userBounds,workBounds,betaRef,mvnCoords,stateNames,
                          model,states,
                          retrySims=0,
                          lambda,
                          errorEllipse,
                          ncores),error=function(e) e))
      if(inherits(tmp,"error")){
        if(ncores==1) cat("FAILED\n")
        if(grepl("Try expanding the extent of the raster",tmp)) simCount <- simCount+1
        else stop(tmp)
      } else {
        simCount <- retrySims
        if(ncores>1) message("\nDONE")
        else message("DONE")
        return(tmp)
      }
    }
    if(ncores>1) message("\nFAILED")
    else message("FAILED")
    stop(tmp)
  }
}
