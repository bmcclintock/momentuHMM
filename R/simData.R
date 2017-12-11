
#' Simulation tool
#'
#' Simulates data from a (multivariate) hidden Markov model. Movement data can be generated with or without observation error attributable to temporal irregularity or location measurement error.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'pois', 'vm', 'weibull', and 'wrpcauchy'. For example,
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
#' @param delta Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}. If \code{formulaDelta} includes covariates, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix, where k is the number of covariates and the columns correspond to states 2:\code{nbStates}. See details below.
#' @param formula Regression formula for the transition probability covariates. Default: \code{~1} (no covariate effect). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}), 
#' and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities.
#' @param formulaDelta Regression formula for the initial distribution. Default: \code{~1} (no covariate effect). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}).
#' @param covs Covariate values to include in the simulated data, as a dataframe. The names of any covariates specified by \code{covs} can
#' be included in \code{formula} and/or \code{DM}. Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL} (the default), and specifying \code{nbCovs>0}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1' and 'cov2' for \code{nbCovs=2}) and can be included in \code{formula} and/or \code{DM}.
#' @param spatialCovs List of \code{\link[raster]{RasterLayer-class}} objects for spatially-referenced covariates. Covariates specified by \code{spatialCovs} are
#' extracted from the raster layer(s) based on the simulated location data for each time step (if applicable). The names of the raster layer(s) can be included in 
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
#' @param obsPerAnimal Either the number of the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated for each animal are uniformously picked from this interval. Alternatively, \code{obsPerAnimal} can be specified as
#' a list of length \code{nbAnimals} with each element providing the number of observations (if single value) or the bounds (if vector of two values) for each individual.
#' Default: \code{c(500,1500)}.
#' @param initialPosition 2-vector providing the x- and y-coordinates of the initial position for all animals. Alternatively, \code{initialPosition} can be specified as
#' a list of length \code{nbAnimals} with each element a 2-vector providing the x- and y-coordinates of the initial position for each individual.
#' Default: \code{c(0,0)}.
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
#' and state-specific formulas (see details). Any formula terms that are not state-specific are included on the parameters for all \code{nbStates} states.
#' @param cons An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. While there could be other uses, primarily intended to constrain specific parameters to be positive. For example, 
#' \code{cons=list(step=c(1,2,1,1))} raises the second parameter to the second power. Default=NULL, which simply raises all parameters to 
#' the power of 1. \code{cons} is ignored for any given data stream unless \code{DM} is specified.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2,dimnames=list(c("mean_1",
#' "mean_2","concentration_1","concentration_2"))))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workcons An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. Warning: use of \code{workcons} is recommended only for advanced users implementing unusual parameter constraints 
#' through a combination of \code{DM}, \code{cons}, and \code{workcons}. \code{workcons} is ignored for any given data stream unless \code{DM} is specified.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param model A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for \code{nbAnimals},
#' \code{obsPerAnimal}, \code{states}, \code{initialPosition}, \code{lambda}, \code{errorEllipse}, and, if covariate values different from those in the data should be specified, 
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
#' as a \code{\link{momentuHMMData}} object suitable for analysis using \code{\link{fitHMM}}.
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
#' \item When covariates are not included in \code{formulaDelta} (i.e. \code{formulaDelta=~1}), then \code{delta} is specified as a vector of length \code{nbStates} that 
#' sums to 1.  When covariates are included in \code{formulaDelta}, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#' 
#' \item State-specific formulas can be specified in \code{DM} using special formula functions. These special functions can take
#' the names \code{paste0("state",1:nbStates)} (where the integer indicates the state-specific formula).  For example, 
#' \code{DM=list(step=list(mean=~cov1+state1(cov2),sd=~cov2+state2(cov1)))} includes \code{cov1} on the mean parameter for all states, \code{cov2}
#' on the mean parameter for state 1, \code{cov2} on the sd parameter for all states, and \code{cov1} on the sd parameter for state 2.
#'
#' \item State- and parameter-specific formulas can be specified for transition probabilities in \code{formula} using special formula functions.
#' These special functions can take the names \code{paste0("state",1:nbStates)} (where the integer indicates the current state from which transitions occur),
#' or \code{paste0("betaCol",nbStates*(nbStates-1))} (where the integer indicates the column of the \code{beta} matrix).  For example with \code{nbStates=3},
#' \code{formula=~cov1+betaCol1(cov2)+state3(cov3)} includes \code{cov1} on all transition probability parameters, \code{cov2} on the \code{beta} column corresponding
#' to the transition from state 1->2, and \code{cov3} on transition probabilities from state 3 (i.e., \code{beta} columns corresponding to state transitions 3->1 and 3->2).
#' 
#' \item Cyclical relationships (e.g., hourly, monthly) may be simulated using the \code{consinor(x,period)} special formula function for covariate \code{x}
#' and sine curve period of time length \code{period}. For example, if 
#' the data are hourly, a 24-hour cycle can be simulated using \code{~cosinor(cov1,24)}, where the covariate \code{cov1} is a repeating series
#' of integers \code{0,1,...,23,0,1,...,23,0,1,...} (note that \code{simData} will not do this for you, the appropriate covariate must be specified using the \code{covs} argument; see example below). 
#' The \code{cosinor(x,period)} function converts \code{x} to 2 covariates
#' \code{cosinorCos(x)=cos(2*pi*x/period)} and \code{consinorSin(x)=sin(2*pi*x/period} for inclusion in the model (i.e., 2 additional parameters per state). The amplitude of the sine wave
#' is thus \code{sqrt(B_cos^2 + B_sin^2)}, where \code{B_cos} and \code{B_sin} are the working parameters correponding to \code{cosinorCos(x)} and \code{cosinorSin(x)}, respectively (e.g., see Cornelissen 2014).
#'
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal', then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
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
#'                 retrySims=10)
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
#' stepDM <- list(mean=~bSpline(time,df=3,degree=0),sd=~1)
#' angleDM <- list(mean=~1,concentration=~bSpline(time,df=3,degree=0))
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
#'                     delta=delta,formulaDelta=formulaDelta,covs=cov)        
#'                 
#' @references
#' 
#' Cornelissen, G. 2014. Cosinor-based rhythmometry. Theoretical Biology and Medical Modelling 11:16.
#'
#' McClintock BT, London JM, Cameron MF, Boveng PL. 2015. Modelling animal movement using the Argos satellite telemetry location error ellipse. 
#' Methods in Ecology and Evolution 6(3):266-277.
#'
#' @export
#' @importFrom stats rnorm runif rmultinom step terms.formula
#' @importFrom raster cellFromXY getValues
#' @importFrom moveHMM simData
#' @importFrom CircStats rvm
#' @importFrom LaplacesDemon rbern
#' @importFrom Brobdingnag as.brob sum

simData <- function(nbAnimals=1,nbStates=2,dist,
                    Par,beta=NULL,delta=NULL,
                    formula=~1,formulaDelta=~1,
                    covs=NULL,nbCovs=0,
                    spatialCovs=NULL,
                    zeroInflation=NULL,
                    oneInflation=NULL,
                    circularAngleMean=NULL,
                    centers=NULL,
                    centroids=NULL,
                    obsPerAnimal=c(500,1500),
                    initialPosition=c(0,0),
                    DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL,stateNames=NULL,
                    model=NULL,states=FALSE,
                    retrySims=0,
                    lambda=NULL,
                    errorEllipse=NULL)
{
  ##############################
  ## Check if !is.null(model) ##
  ##############################
  if(!is.null(model)) {
    
    if(is.miHMM(model)){
      model <- model$miSum
    }
    
    model <- delta_bc(model)
    
    if(is.miSum(model)){
      model$mle <- lapply(model$Par$real,function(x) x$est)
      model$mle$beta <- model$Par$beta$beta$est
      model$mle$delta <- model$Par$real$delta$est
      model$mod <- list()
      model$mod$estimate <- model$MIcombine$coefficients
    }
    
    # extract simulation parameters from model
    nbStates <- length(model$stateNames)
    dist<-model$conditions$dist
    distnames<-names(dist)
    userBounds <- model$conditions$bounds
    stateNames<-model$stateNames
    estAngleMean<-model$conditions$estAngleMean
    circularAngleMean<-model$conditions$circularAngleMean
    DM <- model$conditions$DM
    cons <- model$conditions$cons
    workcons <- model$conditions$workcons
    zeroInflation <- model$conditions$zeroInflation
    oneInflation <- model$conditions$oneInflation
    formula <- model$conditions$formula
    formulaDelta <- model$condition$formulaDelta
  
    Par <- model$mle[distnames]
    parindex <- c(0,cumsum(unlist(lapply(model$conditions$fullDM,ncol)))[-length(model$conditions$fullDM)])
    names(parindex) <- distnames
    for(i in distnames){
      if(!is.null(DM[[i]])){
        Par[[i]] <- model$mod$estimate[parindex[[i]]+1:ncol(model$conditions$fullDM[[i]])]
        names(Par[[i]])<-colnames(model$conditions$fullDM[[i]])
        cons[[i]]<-rep(1,length(cons[[i]]))
        workcons[[i]]<-rep(0,length(workcons[[i]]))
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
        model$conditions$estAngleMean[[i]]<-estAngleMean[[i]]
        model$conditions$userBounds[[i]]<-userBounds[[i]]
        model$conditions$cons[[i]]<-cons[[i]]
        model$conditions$workcons[[i]]<-workcons[[i]]
        model$conditions$DM[[i]]<-DM[[i]]
      }
    }
    beta <- model$mle$beta
    delta <- model$mle$delta
    if(!length(attr(terms.formula(formulaDelta),"term.labels"))){
      delta <- delta[1,]
      rownames(delta)<-NULL
    } else {
      nbCovsDelta <- ncol(model$covsDelta)-1
      foo <- length(model$mod$estimate)-(nbCovsDelta+1)*(nbStates-1)+1
      delta <- matrix(model$mod$estimate[foo:length(model$mod$estimate)],nrow=nbCovsDelta+1)
    }  
    Par<-lapply(Par,function(x) c(t(x)))
    
    if(states) model$data$states <- NULL

    if(is.null(covs)) {
      p<-parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds)
      covNames<-character()
      for(i in distnames){
        covNames<-c(covNames,getCovNames(model,p,i)$DMterms)
      }
      if(!is.null(model$rawCovs)){
        covNames <- c(colnames(model$rawCovs),covNames)
      }
      covNames <- c(covNames,colnames(model$covsDelta)[-1])
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
      
      

      if(length(covsCol)) covs <- model$data[covsCol]
    }
    # else, allow user to enter new values for covariates

  } else {
    if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
    if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
    distnames<-names(dist)
    if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
    Par <- Par[distnames]
    
    mHind <- (is.null(DM) & is.null(userBounds) & is.null(spatialCovs) & is.null(centers) & is.null(centroids) & ("step" %in% names(dist)) & all(unlist(initialPosition)==c(0,0)) & is.null(lambda) & is.null(errorEllipse) & !is.list(obsPerAnimal) & is.null(covs) & !nbCovs & !length(attr(terms.formula(formula),"term.labels")) & !length(attr(terms.formula(formulaDelta),"term.labels")) & is.null(delta)) # indicator for moveHMM::simData
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

  spatialcovnames<-NULL
  if(!is.null(spatialCovs)){
    if(!is.list(spatialCovs)) stop('spatialCovs must be a list')
    spatialcovnames<-names(spatialCovs)
    if(is.null(spatialcovnames)) stop('spatialCovs must be a named list')
    nbSpatialCovs<-length(spatialcovnames)
    if(!("step" %in% distnames)) stop("spatialCovs can only be included when 'step' distribution is specified") 
    else if(!(dist[["step"]] %in% stepdists)) stop("spatialCovs can only be included when valid 'step' distributions are specified") 
    for(j in 1:nbSpatialCovs){
      if(class(spatialCovs[[j]])!="RasterLayer") stop("spatialCovs$",spatialcovnames[j]," must be of class 'RasterLayer'")
      if(any(is.na(raster::getValues(spatialCovs[[j]])))) stop("missing values are not permitted in spatialCovs")
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
      if(length(initialPosition[[i]])!=2 | !is.numeric(initialPosition[[i]])) stop("each element of initialPosition must be a numeric vector of length 2")
    }
  } else {
    if(length(initialPosition)!=2 | !is.numeric(initialPosition)) stop("initialPosition must be a numeric vector of length 2")
    tmpPos<-initialPosition
    initialPosition<-vector('list',nbAnimals)
    for(i in 1:nbAnimals){
      initialPosition[[i]]<-tmpPos
    }
  }

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

  # initial state distribution
  if(is.null(delta)) delta <- rep(1,nbStates)/nbStates

  zeroMass<-oneMass<-vector('list',length(dist))
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
    data[[i]]<-numeric()
  }
  if("angle" %in% distnames){ 
    if(dist[["angle"]] %in% angledists & ("step" %in% distnames))
      if(dist[["step"]] %in% stepdists){
        data$x<-numeric()
        data$y<-numeric()
      }
  } else if("step" %in% distnames){
      if(dist[["step"]] %in% stepdists){
        data$x<-numeric()
        data$y<-numeric()
      }    
  } else if(nbSpatialCovs | length(centerInd) | length(centroidInd)) stop("spatialCovs, centers, and/or centroids cannot be specified without valid step length and turning angle distributions")
  
  #if(is.null(formula)) {
  #  if(allNbCovs) formula <- formula(paste0("~",paste0(c(colnames(allCovs),spatialcovnames),collapse="+")))
  #  else formula <- formula(~1)
  #}
  
  message("=======================================================================")
  message("Simulating HMM with ",nbStates," states and ",length(distnames)," data streams")
  message("-----------------------------------------------------------------------\n")
  for(i in distnames){
    pNames<-p$parNames[[i]]
    if(inputs$circularAngleMean[[i]]) 
      pNames[1]<-paste0("circular ",pNames[1])
    if(is.null(DM[[i]])){
      message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=~1",collapse=", "),")")
    } else if(is.list(DM[[i]])){
      message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=",DM[[i]],collapse=", "),")")
    } else message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,": custom",collapse=", "),")")
  }
  message("\n Transition probability matrix formula: ",paste0(formula,collapse=""))
  message("\n Initial distribution formula: ",paste0(formulaDelta,collapse=""))
  message("=======================================================================")
  
  if(length(all.vars(formula)))
    if(!all(all.vars(formula) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formula' covariate(s) not found")
  if(length(all.vars(formulaDelta)))
    if(!all(all.vars(formulaDelta) %in% c("ID",names(allCovs),centerNames,centroidNames,spatialcovnames)))
      stop("'formulaDelta' covariate(s) not found")
  if(("ID" %in% all.vars(formula) | "ID" %in% all.vars(formulaDelta)) & nbAnimals<2) stop("ID cannot be a covariate when nbAnimals=1")
  
  stateForms<- terms(formula, specials = paste0("state",1:nbStates))
  newformula<-formula
  if(nbStates>1){
    if(length(unlist(attr(stateForms,"specials")))){
      newForm<-attr(stateForms,"term.labels")[-unlist(attr(stateForms,"specials"))]
      for(i in 1:nbStates){
        if(!is.null(attr(stateForms,"specials")[[paste0("state",i)]])){
          for(j in 1:(nbStates-1)){
            newForm<-c(newForm,gsub(paste0("state",i),paste0("betaCol",(i-1)*(nbStates-1)+j),attr(stateForms,"term.labels")[attr(stateForms,"specials")[[paste0("state",i)]]]))
          }
        }
      }
      newformula<-as.formula(paste("~",paste(newForm,collapse="+")))
    }
    formulaStates<-stateFormulas(newformula,nbStates*(nbStates-1),spec="betaCol")
    if(length(unlist(attr(terms(newformula, specials = c(paste0("betaCol",1:(nbStates*(nbStates-1))),"cosinor")),"specials")))){
      allTerms<-unlist(lapply(formulaStates,function(x) attr(terms(x),"term.labels")))
      newformula<-as.formula(paste("~",paste(allTerms,collapse="+")))
      formterms<-attr(terms.formula(newformula),"term.labels")
    } else {
      formterms<-attr(terms.formula(newformula),"term.labels")
      newformula<-formula
    }
  }
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
  nbBetaCovs <- ncol(model.matrix(newformula,tmpCovs))

  if(is.null(beta))
    beta <- matrix(rnorm(nbStates*(nbStates-1)*nbBetaCovs),nrow=nbBetaCovs)
  else {
    if(ncol(beta)!=nbStates*(nbStates-1) | nrow(beta)!=nbBetaCovs) {
      error <- paste("beta has wrong dimensions: it should have",nbBetaCovs,"rows and",
                     nbStates*(nbStates-1),"columns.")
      stop(error)
    }
  }
  if(nbStates>1){
    for(state in 1:(nbStates*(nbStates-1))){
      noBeta<-which(match(colnames(model.matrix(newformula,tmpCovs)),colnames(model.matrix(formulaStates[[state]],tmpCovs)),nomatch=0)==0)
      if(length(noBeta)) beta[noBeta,state] <- 0
    }
  }
  
  covsDelta <- model.matrix(formulaDelta,tmpCovs)
  nbCovsDelta <- ncol(covsDelta)-1
  if(!nbCovsDelta){
    if(length(delta) != (nbCovsDelta+1)*nbStates)
      stop(paste("delta has the wrong length: it should have",nbStates,"elements."))
    deltaB <- log(delta[-1]/delta[1])
  } else {
    if(is.null(dim(delta)))
      stop(paste("delta has wrong dimensions: it should have",nbCovsDelta+1,"rows and",
                 nbStates-1,"columns."))
    if(ncol(delta)!=nbStates-1 | nrow(delta)!=nbCovsDelta+1)
      stop(paste("delta has wrong dimensions: it should have",nbCovsDelta+1,"rows and",
                 nbStates-1,"columns."))
    deltaB <- delta
  }
  
  if(!nbSpatialCovs | !retrySims){
  
    ###########################
    ## Loop over the animals ##
    ###########################
    for (zoo in 1:nbAnimals) {
  
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
  
      X <- matrix(0,nrow=nbObs,ncol=2)
      X[1,] <- initialPosition[[zoo]] # initial position of animal
  
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
      
      if(!nbSpatialCovs & !length(centerInd) & !length(centroidInd)) {
        DMcov <- model.matrix(newformula,subCovs)
        gFull <-  DMcov %*% beta
        
        # format parameters
        DMinputs<-getDM(subCovs,inputs$DM,dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
        fullDM <- DMinputs$fullDM
        DMind <- DMinputs$DMind
        wpar <- n2w(Par,bounds,beta,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
        nc <- meanind <- vector('list',length(distnames))
        names(nc) <- names(meanind) <- distnames
        for(i in distnames){
          nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
          if(inputs$circularAngleMean[[i]]) meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
        }
        covsDelta <- model.matrix(formulaDelta,subCovs[1,,drop=FALSE])
        fullsubPar <- w2n(wpar,bounds,parSize,nbStates,length(attr(terms.formula(newformula),"term.labels")),inputs$estAngleMean,inputs$circularAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nbObs,dist,p$Bndind,nc,meanind,covsDelta)
        g <- gFull[1,,drop=FALSE]
        delta <- fullsubPar$delta
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
            subCovs[1,centerNames[(j-1)*2+1:2]]<-distAngle(X[1,],X[1,],centers[centerInd[j],])
          }
        }
        
        if(length(centroidInd)){
          for(j in 1:centroidInd){
            subCovs[1,centroidNames[(j-1)*2+1:2]]<-distAngle(X[1,],X[1,],as.numeric(centroids[[j]][1,]))
          }
        }
        g <- model.matrix(newformula,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE])) %*% beta
        covsDelta <- model.matrix(formulaDelta,cbind(subCovs[1,,drop=FALSE],subSpatialcovs[1,,drop=FALSE]))
        delta <- c(rep(0,nbCovsDelta+1),deltaB)
        deltaXB <- covsDelta %*% matrix(delta,nrow=nbCovsDelta+1)
        expdelta <- exp(deltaXB)
        delta <- expdelta/rowSums(expdelta)
        for(i in which(!is.finite(rowSums(delta)))){
          tmp <- exp(Brobdingnag::as.brob(deltaXB[i,]))
          delta[i,] <- as.numeric(tmp/Brobdingnag::sum(tmp))
        }
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
        
        if(nbSpatialCovs |  length(centerInd) | length(centroidInd)){
          # format parameters
          DMinputs<-getDM(cbind(subCovs[k,,drop=FALSE],subSpatialcovs[k,,drop=FALSE]),inputs$DM,dist,nbStates,p$parNames,p$bounds,Par,cons,workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
          fullDM <- DMinputs$fullDM
          DMind <- DMinputs$DMind
          wpar <- n2w(Par,bounds,beta,deltaB,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
          nc <- meanind <- vector('list',length(distnames))
          names(nc) <- names(meanind) <- distnames
          for(i in distnames){
            nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
            if(inputs$circularAngleMean[[i]]) meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
          }
          subPar <- w2n(wpar,bounds,parSize,nbStates,length(attr(terms.formula(newformula),"term.labels")),inputs$estAngleMean,inputs$circularAngleMean,stationary=FALSE,DMinputs$cons,fullDM,DMind,DMinputs$workcons,1,dist,p$Bndind,nc,meanind,covsDelta)
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
  
          for(j in 1:(parSize[[i]]-zeroInflation[[i]]-oneInflation[[i]]))
            genArgs[[i]][[j+1]] <- subPar[[i]][(j-1)*nbStates+Z[k]]
          
          if(dist[[i]] %in% angledists){
            
            genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
            if(genData[[i]][k] >  pi) genData[[i]][k] <- genData[[i]][k]-2*pi
            if(genData[[i]][k] < -pi) genData[[i]][k] <- genData[[i]][k]+2*pi
  
            if(i=="angle" & ("step" %in% distnames)){
              if(dist[["step"]] %in% stepdists) {
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
            
            if(dist[[i]]=="gamma") {
              shape <- genArgs[[i]][[2]]^2/genArgs[[i]][[3]]^2
              scale <- genArgs[[i]][[3]]^2/genArgs[[i]][[2]]
              genArgs[[i]][[2]] <- shape
              genArgs[[i]][[3]] <- 1/scale # rgamma expects rate=1/scale
            }
      
            probs <- c(1.-zeroMass[[i]][Z[k]]-oneMass[[i]][Z[k]],zeroMass[[i]][Z[k]],oneMass[[i]][Z[k]])
            rU <- which(rmultinom(1,1,prob=probs)==1)
            if(rU==1)
              genData[[i]][k] <- do.call(Fun[[i]],genArgs[[i]])
            else if(rU==2)
              genData[[i]][k] <- 0
            else
              genData[[i]][k] <- 1
          }
          
          d[[i]] <- genData[[i]]
          
        }
        # get next state
        gamma <- diag(nbStates)
        if(nbSpatialCovs | length(centerInd) | length(centroidInd)){
          if(nbSpatialCovs){
            for(j in 1:nbSpatialCovs){
              getCell<-raster::cellFromXY(spatialCovs[[j]],c(X[k+1,1],X[k+1,2]))
              if(is.na(getCell)) stop("Movement is beyond the spatial extent of the ",spatialcovnames[j]," raster. Try expanding the extent of the raster.")
              subSpatialcovs[k+1,j]<-spatialCovs[[j]][getCell]
            }
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
          g <- model.matrix(newformula,cbind(subCovs[k+1,,drop=FALSE],subSpatialcovs[k+1,,drop=FALSE])) %*% beta
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
      
      if("angle" %in% distnames){ 
        if(dist[["angle"]] %in% angledists & ("step" %in% distnames))
          if(dist[["step"]] %in% stepdists){
            d$angle[1] <- NA # the first angle value is arbitrary
            step0 <- which(d$step==0)
            d$angle[c(step0,step0+1)] <- NA
            #if(length(centerInd)) subCovs[1,centerNames[seq(2,2*length(centerInd),2)]] <- NA
            d$x=X[,1]
            d$y=X[,2]
          }
      } else if("step" %in% distnames){
        if(dist[["step"]] %in% stepdists){
          d$x=c(0,cumsum(d$step)[-nrow(d)])
          d$y=X[,2]
        }    
      }
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
      if(dist[[i]] %in% angledists)
        class(data[[i]]) <- c(class(data[[i]]), "angle")
    }
    
    # account for observation error (if any)
    out<-simObsData(momentuHMMData(data),lambda,errorEllipse)
    
    message("DONE")
    return(out)
  } else {
    simCount <- 0
    cat("Attempting to simulate tracks within spatial extent(s) of raster layers(s). Press 'esc' to force exit from 'simData'\n",sep="")
    while(simCount < retrySims){
      cat("\r    Attempt ",simCount+1," of ",retrySims,"...",sep="")
      tmp<-suppressMessages(tryCatch(simData(nbAnimals,nbStates,dist,
                          Par,beta,delta,
                          formula,formulaDelta,
                          covs,nbCovs,
                          spatialCovs,
                          zeroInflation,
                          oneInflation,
                          circularAngleMean,
                          centers,
                          centroids,
                          obsPerAnimal,
                          initialPosition,
                          DM,cons,userBounds,workcons,stateNames,
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
