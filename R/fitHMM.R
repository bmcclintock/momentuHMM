
#' Fit a multivariate HMM to the data
#'
#' Fit a (multivariate) hidden Markov model to the data provided, using numerical optimization of the log-likelihood
#' function.
#'
#' @param data An object \code{momentuHMMData}.
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'gamma','weibull','exp','lnorm','beta','pois','wrpcauchy', and 'vm'. For example,
#' \code{dist=list(step='gamma', angle='vm', dives='pois')} indicates 3 data streams ('step', 'angle', and 'dives')
#' and their respective probability distributions ('gamma', 'vm', and 'pois').  The names of the data streams 
#' (e.g., 'step', 'angle', 'dives') must match component names in \code{data}.
#' @param Par0 A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be in the order expected by the pdfs of \code{dist}, 
#' and any zero-mass and/or one-mass parameters should be the last (if both are present, then zero-mass parameters must preceed one-mass parameters). 
#' Note that zero-mass parameters are mandatory if there are zeros in 
#' data streams with a 'gamma','weibull','exp','lnorm', or 'beta' distribution, and one-mass parameters are mandatory if there are ones in 
#' data streams with a 'beta' distribution.
#' For example, for a 2-state model using the Von Mises (vm) distribution for a data stream named 'angle' and 
#' the zero-inflated gamma distribution for a data stream named 'step', the vector of initial parameters would be something like: 
#' \code{Par0=list(step=c(mean_1,mean_2,sd_1,sd_2,zeromass_1,zeromass_2), angle=c(mean_1,mean_2,concentration_1,concentration_2))}.
#' 
#' If \code{DM} is not specified for a given data stream, then \code{Par0} is on the natural (i.e., real) scale of the parameters.  
#' However, if \code{DM} is specified for a given data stream, then \code{Par0} must be on the working (i.e., beta) scale of the 
#' parameters, and the length of \code{Par0} must match the number of columns in the design matrix.  See details below.
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities (more
#' information in 'Details').
#' Default: \code{NULL}. If not specified, \code{beta0} is initialized such that the diagonal elements
#' of the transition probability matrix are dominant.
#' @param delta0 Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). For example, \code{estAngleMean=list(angle=TRUE)} indicates the angle mean is to be 
#' estimated for 'angle'.  Default is \code{NULL}, which assumes any angle means are fixed to zero and are not to be estimated. 
#' Any \code{estAngleMean} elements corresponding to data streams that do not have angular distributions are ignored.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  For example, 
#' \code{circularAngleMean=list(angle=TRUE)} indicates the angle mean is be estimated for 'angle' using circular-circular 
#' regression.  Whenever circular-circular regression is used for an angular data stream, a corresponding design matrix (\code{DM}) 
#' must be specified for the data stream, and the previous movement direction (i.e., a turning angle of zero) is automatically used 
#' as the reference angle (i.e., the intercept). Any circular-circular regression covariates in \code{data} should therefore be relative to the previous 
#' direction of movement (instead of true north; see \code{\link{prepData}} and \code{\link{circAngles}}).  See Duchesne et al. (2015) for specifics on the circular-circular regression model 
#' using previous movement direction as the reference angle. Default is \code{NULL}, which assumes circular-linear regression is 
#' used for any angular distributions for which the mean angle is to be estimated. \code{circularAngleMean} elements corresponding to angular data 
#' streams are ignored unless the corresponding element of \code{estAngleMean} is \code{TRUE}. Any \code{circularAngleMean} elements 
#' corresponding to data streams that do not have angular distributions are ignored.
#' @param formula Regression formula for the transition probability covariates. Default: \code{~1} (no covariate effects).
#' @param stationary \code{FALSE} if there are covariates. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param verbose Determines the print level of the optimizer. The default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed.
#' @param nlmPar List of parameters to pass to the optimization function \code{nlm} (which should be either
#' '\code{gradtol}', '\code{stepmax}', '\code{steptol}', or '\code{iterlim}' -- see \code{nlm}'s documentation
#' for more detail)
#' @param fit \code{TRUE} if an HMM should be fitted to the data, \code{FALSE} otherwise.
#' If fit=\code{FALSE}, a model is returned with the MLE replaced by the initial parameters given in
#' input. This option can be used to assess the initial parameters, parameter bounds, etc. Default: \code{TRUE}.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a matrix.  For example, for a 2-state 
#' model using the gamma distribution for a data stream named 'step', \code{DM=list(step=list(mean=~cov1, sd=~1))} specifies the mean 
#' parameters as a function of the covariate 'cov1' for each state.  This model could equivalently be specified as a 4x6 matrix using 
#' character strings for the covariate: 
#' \code{DM=list(step=matrix(c(1,0,0,0,'cov1',0,0,0,0,1,0,0,0,'cov1',0,0,0,0,1,0,0,0,0,1),4,6))}
#' where the 4 rows correspond to the state-dependent paramaters (mean_1,mean_2,sd_1,sd_2) and the 6 columns correspond to the regression 
#' coefficients. 
#' 
#' Design matrices specified using formulas automatically include effects on the parameters for all \code{nbStates} states. Thus 
#' for models including covariates on a subset of the state-dependent parameters, the design matrix must be specified manually using matrices or through
#' a combination of formulas and \code{fixPar} (i.e., fixing the working parameters to zero ).
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
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any). Default: NULL (states are not known). This should be a vector with length the number
#' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
#' the state is not known.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. Default: NULL 
#' (no parameters are fixed). Each element of \code{fixPar} should be a vector of the same length as the corresponding vector of 
#' \code{Par0}; each element should either be numeric (the fixed value of the parameter) or NA if the parameter is to be estimated. 
#'
#' @return A \code{\link{momentuHMM}} object, i.e. a list of:
#' \item{mle}{A named list of the maximum likelihood estimates of the parameters of the model (if the numerical algorithm
#' has indeed identified the global maximum of the likelihood function). Elements are included for the parameters of each
#' data strea, as well as \code{beta} (transition probabilities regression coefficients - more information
#' in 'Details'), \code{gamma} (transition probabilities on real scale, based on mean covariate values if \code{formula}
#' includes covariates), and \code{delta} (initial distribution).}
#' \item{CIreal}{Standard errors and 95\% confidence intervals on the real (i.e., natural) scale of parameters}
#' \item{CIbeta}{Standard errors and 95\% confidence intervals on the beta (i.e., working) scale of parameters}
#' \item{data}{The momentuHMMData object}
#' \item{mod}{The object returned by the numerical optimizer \code{nlm}}
#' \item{conditions}{Conditions used to fit the model, e.g., \code{bounds} (parameter bounds), distributions, \code{zeroInflation},
#' \code{estAngleMean}, \code{stationary}, \code{formula}, \code{DM}, \code{fullDM} (full design matrix), etc.)}
#' \item{rawCovs}{Raw covariate values for transition probabilities, as found in the data (if any). Used in \code{\link{plot.momentuHMM}}.}
#' \item{stateNames}{The names of the states.}
#' \item{knownStates}{Vector of values of the state process which are known.}
#'
#' @details
#' \itemize{
#' \item The matrix \code{beta} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix -
#' filled in row-wise).
#' In a covariate-free model (default), \code{beta} has one row, for the intercept.
#'
#' \item The choice of initial parameters (particularly \code{Par0} and \code{beta0}) is crucial to fit a model. The algorithm might not find
#' the global optimum of the likelihood function if the initial parameters are poorly chosen.
#' 
#' \item If \code{DM} is specified for a particular data stream, then the initial values are specified on 
#' the working (i.e., beta) scale of the parameters. The working scale of each parameter is determined by the link function used.
#' If a parameter P is bound by (0,Inf) then the working scale is the log(P) scale.  If the parameter bounds are (-pi,pi) then the working 
#' scale is tan(P/2) unless circular-circular regression is used. Otherwise if the parameter bounds are finite then logit(P) is the working scale.  
#' The function \code{\link{getParDM}} is intended to help with obtaining initial values on the working scale when specifying a design matrix and other 
#' parameter constraints (see example below). When circular-circular regression is specified using \code{circularAngleMean}, the working scale 
#' for the mean turning angle is not as easily interpretable, but the 
#' link function is atan2(sin(X)*B,1+cos(X)*B), where X are the angle covariates and B the angle coefficients (see Duchesne et al. 2015). 
#' Under this formulation, the reference turning angle is 0 (i.e., movement in the same direction as the previous time step). 
#' In other words, the mean turning angle is zero when the coefficient(s) B=0.
#' 
#' \item Circular-circular regression in \code{momentuHMM} is designed for turning angles (not bearings) as computed by \code{\link{simData}} and \code{\link{prepData}}. 
#' Any circular-circular regression covariates for time step t should therefore be relative to the previous 
#' direction of movement for time step t-1 (instead of true north).  In other words, circular-circular regression covariates for time step t should be the turning angle
#' between the direction of movement for time step t-1 and the bearing of the covariate relative to true north for time step t.  If provided bearings relative to true north, 
#' \code{\link{circAngles}} or \code{\link{prepData}} can perform this calculation for you.  
#' }
#' 
#' @seealso \code{\link{getParDM}}, \code{\link{prepData}}, \code{\link{simData}}
#'
#' @examples
#' nbStates <- 2
#' stepDist <- "gamma" # step distribution
#' angleDist <- "vm" # turning angle distribution
#' 
#' # extract data from momentuHMM example
#' data <- example$data
#'
#' ### 1. fit the model to the simulated data
#' # define initial values for the parameters
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
#' anglePar <- kappa0 # not estimating angle mean, so not included
#' formula <- ~cov1+cos(cov2)
#'
#' m <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'             Par0=list(step=stepPar,angle=anglePar),formula=formula)
#'
#' print(m)
#' 
#' \dontrun{
#' ### 2. fit the exact same model to the simulated data using DM formulas
#' # Get initial values for the parameters on working scale
#' Par0 <- getParDM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'         Par=list(step=stepPar,angle=anglePar),
#'         DM=list(step=list(mean=~1,sd=~1),angle=list(concentration=~1)))
#'
#' mDMf <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'               Par0=Par0,formula=formula,
#'               DM=list(step=list(mean=~1,sd=~1),angle=list(concentration=~1)))
#'
#' print(mDMf)
#' 
#' ### 3. fit the exact same model to the simulated data using DM matrices
#' # define DM
#' DMm <- list(step=diag(4),angle=diag(2))
#' 
#' # user-specified dimnames not required but are recommended
#' dimnames(DMm$step) <- list(c("mean_1","mean_2","sd_1","sd_2"),
#'                    c("mean_1:(Intercept)","mean_2:(Intercept)",
#'                    "sd_1:(Intercept)","sd_2:(Intercept)"))
#' dimnames(DMm$angle) <- list(c("concentration_1","concentration_2"),
#'                     c("concentration_1:(Intercept)","concentration_2:(Intercept)"))
#'                   
#' mDMm <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'               Par0=Par0,formula=formula,
#'               DM=DMm)
#'
#' print(mDMm)
#' 
#' ### 4. fit step mean parameter covariate model to the simulated data using DM
#' stepDMf <- list(mean=~cov1,sd=~1)
#' Par0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
#'                  Par=list(step=stepPar,angle=anglePar),
#'                  DM=list(step=stepDMf,angle=DMm$angle))
#' mDMfcov <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'               Par0=Par0,
#'               formula=formula,
#'               DM=list(step=stepDMf,angle=DMm$angle))
#'
#' print(mDMfcov)
#' 
#' ### 5. fit the exact same step mean parameter covariate model using DM matrix
#' stepDMm <- matrix(c(1,0,0,0,"cov1",0,0,0,0,1,0,0,0,"cov1",0,0,
#'                  0,0,1,0,0,0,0,1),4,6,dimnames=list(c("mean_1","mean_2","sd_1","sd_2"),
#'                  c("mean_1:(Intercept)","mean_1:cov1","mean_2:(Intercept)","mean_2:cov1",
#'                  "sd_1:(Intercept)","sd_2:(Intercept)")))
#' Par0 <- getParDM(data,nbStates,list(step=stepDist,angle=angleDist),
#'                  Par=list(step=stepPar,angle=anglePar),
#'                  DM=list(step=stepDMm,angle=DMm$angle))
#' mDMmcov <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'               Par0=Par0,
#'               formula=formula,
#'               DM=list(step=stepDMm,angle=DMm$angle))
#'
#' print(mDMmcov)
#' 
#' ### 6. fit circular-circular angle mean covariate model to the simulated data using DM
#'
#' # Generate fake circular covariate using circAngles
#' data$cov3 <- circAngles(refAngle=2*atan(rnorm(nrow(data))),data)
#' 
#' # Fit circular-circular regression model for angle mean
#' # Note no intercepts are estimated for angle means because these are by default
#' # the previous movement direction (i.e., a turning angle of zero)
#' mDMcircf <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                  Par0=list(step=stepPar,angle=c(0,0,Par0$angle)),
#'                   formula=formula,
#'                   estAngleMean=list(angle=TRUE),
#'                   circularAngleMean=list(angle=TRUE),
#'                   DM=list(angle=list(mean=~cov3,concentration=~1)))
#'                   
#' print(mDMcircf)
#'                   
#' ### 7. fit the exact same circular-circular angle mean model using DM matrices
#' 
#' # Note no intercept terms are included in DM for angle means because the intercept is
#' # by default the previous movement direction (i.e., a turning angle of zero)
#' mDMcircm <- fitHMM(data=data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                  Par0=list(step=stepPar,angle=c(0,0,Par0$angle)),
#'                   formula=formula,
#'                   estAngleMean=list(angle=TRUE),
#'                   circularAngleMean=list(angle=TRUE),
#'                   DM=list(angle=matrix(c("cov3",0,0,0,0,"cov3",0,0,0,0,1,0,0,0,0,1),4,4)))
#'                   
#' print(mDMcircm)
#' }
#' 
#' @references
#' 
#' Duchesne, T., Fortin, D., Rivest L-P. 2015. Equivalence between step selection functions and 
#' biased correlated random walks for statistical inference on animal movement. PLoS ONE 10 (4):
#' e0122947.
#' 
#' Langrock R., King R., Matthiopoulos J., Thomas L., Fortin D., Morales J.M. 2012.
#' Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions.
#' Ecology, 93 (11), 2336-2342.
#' 
#' McClintock B.T., King R., Thomas L., Matthiopoulos J., McConnell B.J., Morales J.M. 2012. A general 
#' discrete-time modeling framework for animal movement using multistate random walks. Ecological 
#' Monographs, 82 (3), 335-349.
#' 
#' McClintock B.T., Russell D.J., Matthiopoulos J., King R. 2013. Combining individual animal movement 
#' and ancillary biotelemetry data to investigate population-level activity budgets. Ecology, 94 (4), 838-849.
#' 
#' Patterson T.A., Basson M., Bravington M.V., Gunn J.S. 2009.
#' Classifying movement behaviour in relation to environmental conditions using hidden Markov models.
#' Journal of Animal Ecology, 78 (6), 1113-1123.
#' 
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats model.matrix get_all_vars nlm terms terms.formula
#' @importFrom CircStats dwrpcauchy dvm pvm
#'
#' @useDynLib momentuHMM

fitHMM <- function(data,nbStates,dist,
                   Par0,beta0=NULL,delta0=NULL,
                   estAngleMean=NULL,circularAngleMean=NULL,
                   formula=~1,stationary=FALSE,
                   verbose=0,nlmPar=NULL,fit=TRUE,
                   DM=NULL,cons=NULL,userBounds=NULL,workcons=NULL,
                   stateNames=NULL,knownStates=NULL,fixPar=NULL)
{
  
  #####################
  ## Check arguments ##
  #####################
  
  # check that the data is a momentuHMMData object
  if(!is.momentuHMMData(data))
    stop("'data' must be a momentuHMMData object (as output by prepData or simData)")
  
  if(length(data)<1 | any(dim(data)<1))
    stop("The data input is empty.")

  if(!is.formula(formula))
    stop("Check the argument 'formula'.")
  
  if(nbStates<1) stop('nbStates must be >0')

  # check that there is no response varibale in the formula
  if(attr(terms(formula),"response")!=0)
    stop("The response variable should not be specified in the formula.")
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par0) | is.null(names(Par0))) stop("'Par0' must be a named list")
  distnames<-names(dist)
  if(any(is.na(match(distnames,names(data))))) stop(paste0(distnames[is.na(match(distnames,names(data)))],collapse=", ")," not found in data")
  if(!all(distnames %in% names(Par0))) stop(paste0(distnames[which(!(distnames %in% names(Par0)))],collapse=", ")," missing in 'Par0'")
  Par0 <- Par0[distnames]
  
  # build design matrix for t.p.m.
  covsCol <- get_all_vars(formula,data)#rownames(attr(terms(formula),"factors"))#attr(terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
  covs <- model.matrix(formula,data)
  if(nrow(covs)!=nrow(data)) stop("covariates cannot contain missing values")
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  # check that stationary==FALSE if there are covariates
  if(nbCovs>0 & stationary==TRUE)
    stop("stationary can't be set to TRUE if there are covariates.")
  
  if(length(knownStates) > 0){
    if(length(knownStates) != nrow(data)) 
      stop("'knownStates' should be of same length as the data, i.e. ",nrow(data))
    if(!all(is.na(knownStates))) {
      if(max(knownStates, na.rm = TRUE) > nbStates | min(knownStates, na.rm = TRUE) < 1 | !isTRUE(all.equal(knownStates,as.integer(knownStates)))) 
        stop("'knownStates' should only contain integers between 1 and ", nbStates, " (or NAs)")
    }
  }

  if(length(covsCol)>0) {
    rawCovs <- covsCol
  }
  else {
    rawCovs <- NULL
  }
  
  # check that observations are within expected bounds
  for(i in which(unlist(lapply(dist,function(x) x %in% nonnegativedists))))
    if(length(which(data[[distnames[[i]]]]<0))>0)
      stop(distnames[[i]]," data should be non-negative")
  
  for(i in which(unlist(lapply(dist,function(x) x %in% angledists))))
    if(length(which(data[[distnames[[i]]]] < -pi | data[[distnames[[i]]]] > pi))>0)
      stop(distnames[[i]]," data should be between -pi and pi")
  
  for(i in which(unlist(lapply(dist,function(x) x %in% "beta"))))
    if(length(which(data[[distnames[[i]]]]<0 | data[[distnames[[i]]]]>1))>0)
      stop(distnames[[i]]," data should be between 0 and 1")

  for(i in which(unlist(lapply(dist,function(x) x %in% "pois"))))
    if(!isTRUE(all.equal(data[[distnames[[i]]]],as.integer(data[[distnames[[i]]]]))))
      stop(distnames[[i]]," data should be non-negative integers")
  
  # determine whether zero-inflation or one-inflation should be included
  zeroInflation <- oneInflation <- vector('list',length(distnames))
  names(zeroInflation) <- names(oneInflation) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% zeroInflationdists){
      if(length(which(data[[i]]==0))>0) {
        zeroInflation[[i]]<-TRUE
      }
      else 
        zeroInflation[[i]]<-FALSE
    }
    else zeroInflation[[i]]<-FALSE
    if(dist[[i]] %in% oneInflationdists){
      if(length(which(data[[i]]==1))>0) {
        oneInflation[[i]]<-TRUE
      }
      else 
        oneInflation[[i]]<-FALSE
    }
    else oneInflation[[i]]<-FALSE
    if(is.null(DM[[i]])){
      if(zeroInflation[[i]]){
        # check that zero-mass is in the open interval (0,1)
        zm0 <- Par0[[i]][(length(Par0[[i]])-nbStates+1):length(Par0[[i]])]
        zm0[which(zm0==0)] <- 1e-8
        zm0[which(zm0==1)] <- 1-1e-8
        Par0[[i]][(length(Par0[[i]])-nbStates*oneInflation[[i]]-nbStates+1):(length(Par0[[i]])-nbStates*oneInflation[[i]])] <- zm0
      }
      if(oneInflation[[i]]){
        # check that one-mass is in the open interval (0,1)
        om0 <- Par0[[i]][(length(Par0[[i]])-nbStates+1):length(Par0[[i]])]
        om0[which(om0==0)] <- 1e-8
        om0[which(om0==1)] <- 1-1e-8
        Par0[[i]][(length(Par0[[i]])-nbStates+1):length(Par0[[i]])] <- om0
      }
    }
  }

  mHind <- (is.null(DM) & is.null(userBounds) & ("step" %in% distnames) & is.null(fixPar) & !length(attr(terms.formula(formula),"term.labels")) & stationary) # indicator for moveHMMwrap below
  
  inputs <- checkInputs(nbStates,dist,Par0,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,cons,workcons,stateNames)
  p <- inputs$p
  
  DMinputs<-getDM(data,inputs$DM,dist,nbStates,p$parNames,p$bounds,Par0,inputs$cons,inputs$workcons,zeroInflation,oneInflation,inputs$circularAngleMean)
  fullDM <- DMinputs$fullDM
  DMind <- DMinputs$DMind
  
  if(!is.null(beta0)) {
    if(is.null(dim(beta0)))
      stop(paste("beta0 has wrong dimensions: it should have",nbCovs+1,"rows and",
                     nbStates*(nbStates-1),"columns."))
    if(ncol(beta0)!=nbStates*(nbStates-1) | nrow(beta0)!=nbCovs+1)
      stop(paste("beta0 has wrong dimensions: it should have",nbCovs+1,"rows and",
                     nbStates*(nbStates-1),"columns."))
  }

  if(!is.null(delta0))
    if(length(delta0)!=nbStates)
      stop(paste("delta0 has the wrong length: it should have",nbStates,"elements."))
  
  # check that verbose is in {0,1,2}
  if(!(verbose %in% c(0,1,2)))
    stop("verbose must be in {0,1,2}")

  # check elements of nlmPar
  lsPars <- c("gradtol","stepmax","steptol","iterlim")
  if(length(which(!(names(nlmPar) %in% lsPars)))>0)
    stop("Check the names of the element of 'nlmPar'; they should be in
         ('gradtol','stepmax','steptol','iterlim')")

  ####################################
  ## Prepare initial values for nlm ##
  ####################################
  if(is.null(beta0) & nbStates>1) {
    beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                    nrow=nbCovs+1,byrow=TRUE)
  }

  if(is.null(delta0))
    delta0 <- rep(1,nbStates)/nbStates
  if(stationary)
    delta0 <- NULL

  # build the vector of initial working parameters
  wparIndex <- numeric()
  if(!is.null(fixPar)){
    parindex <- c(0,cumsum(unlist(lapply(Par0,length))))
    names(parindex) <- c(distnames,"beta")
    for(i in distnames){
      if(!is.null(fixPar[[i]])){
        if(length(fixPar[[i]])!=length(Par0[[i]])) stop("fixPar$",i," must be of length ",length(Par0[[i]]))
        tmp <- which(!is.na(fixPar[[i]]))
        Par0[[i]][tmp]<-fixPar[[i]][tmp]
        wparIndex <- c(wparIndex,parindex[[i]]+tmp)
      } else {
        fixPar[[i]] <- rep(NA,length(Par0[[i]]))
        #wparIndex <- c(wparIndex,parindex[[i]]+1:length(Par0[[i]]))
      }
    }
    if(nbStates>1){
      if(!is.null(fixPar$beta)){
        if(length(fixPar$beta)!=length(beta0)) stop("fixPar$beta must be of length ",length(beta0))
        tmp <- which(!is.na(fixPar$beta))
        beta0[tmp]<-fixPar$beta[tmp]
        wparIndex <- c(wparIndex,parindex[["beta"]]+tmp)
      } else {
        fixPar$beta <- rep(NA,length(beta0))
      }
    }
    if(!is.null(fixPar$delta)){
      if(stationary & any(!is.na(fixPar$delta))) stop("delta cannot be fixed when stationary=TRUE")
      else if(!stationary) {
        tmp <- which(!is.na(fixPar$delta))
        if(length(tmp)){
          delta0[tmp] <- fixPar$delta[tmp]
          if(length(tmp)!=length(delta0) | sum(delta0)!=1) stop("fixPar$delta must sum to 1")
          wparIndex <- c(wparIndex,parindex[["beta"]]+length(beta0)+tmp)
        }
      }
    } else {
      fixPar$delta <- rep(NA,length(delta0))
    }
  } else {
    fixPar <- vector('list',length(distnames))
    names(fixPar) <- distnames
    for(i in distnames){
      fixPar[[i]] <- rep(NA,length(Par0[[i]]))
    }
    fixPar$beta <- rep(NA,length(beta0))
    fixPar$delta <- rep(NA,length(delta0))
  }
  fixPar <- fixPar[c(distnames,"beta","delta")]
  
  wpar <- n2w(Par0,p$bounds,beta0,delta0,nbStates,inputs$estAngleMean,inputs$DM,DMinputs$cons,DMinputs$workcons,p$Bndind)
  if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")

  
  ##################
  ## Optimization ##
  ##################
  # this function is used to muffle the warning "NA/Inf replaced by maximum positive value" in nlm and "value out of range in 'lgamma'" in nLogLike_rcpp
  h <- function(w) {
    if(any(grepl("NA/Inf replaced by maximum positive value",w)) | any(grepl("value out of range in 'lgamma'",w)))
      invokeRestart("muffleWarning")
  }
  
  message("=======================================================================")
  message("Fitting HMM with ",nbStates," states and ",length(distnames)," data streams")
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
  message("=======================================================================")

  # just use moveHMM if simpler models are specified
  if(all(distnames %in% c("step","angle")) & mHind){
    out<-moveHMMwrap(data,nbStates,dist,Par0,beta0,delta0,inputs$estAngleMean,formula,stationary,verbose,nlmPar,fit)
    mod<-out$mod
    mle<-out$mle
  } 
  else if(fit) {
    # check additional parameters for nlm
    gradtol <- ifelse(is.null(nlmPar$gradtol),1e-6,nlmPar$gradtol)
    typsize = rep(1, length(wpar))
    defStepmax <- max(1000 * sqrt(sum((wpar/typsize)^2)),1000)
    stepmax <- ifelse(is.null(nlmPar$stepmax),defStepmax,nlmPar$stepmax)
    steptol <- ifelse(is.null(nlmPar$steptol),1e-6,nlmPar$steptol)
    iterlim <- ifelse(is.null(nlmPar$iterlim),1000,nlmPar$iterlim)

    # call to optimizer nlm
    withCallingHandlers(mod <- nlm(nLogLike,wpar,nbStates,formula,p$bounds,p$parSize,data,dist,covs,
                                   inputs$estAngleMean,inputs$circularAngleMean,zeroInflation,oneInflation,
                                   stationary,DMinputs$cons,fullDM,DMind,DMinputs$workcons,p$Bndind,knownStates,unlist(fixPar),wparIndex,
                                   print.level=verbose,gradtol=gradtol,
                                   stepmax=stepmax,steptol=steptol,
                                   iterlim=iterlim,hessian=TRUE),
                        warning=h) # filter warnings using function h

    # convert the parameters back to their natural scale
    wpar <- mod$estimate
    mle <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,inputs$circularAngleMean,stationary,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nrow(data),dist,p$Bndind)
  }
  else {
    mod <- NA
    mle <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,inputs$circularAngleMean,stationary,DMinputs$cons,fullDM,DMind,DMinputs$workcons,nrow(data),dist,p$Bndind)
  }

  ####################
  ## Prepare output ##
  ####################
  if(is.null(stateNames)){
    for(i in 1:nbStates)
      stateNames[i] <- paste("state",i)
  }
  
  # name columns and rows of MLEs
  parindex <- c(0,cumsum(unlist(lapply(fullDM,ncol)))[-length(fullDM)])
  names(parindex) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% angledists)
      if(!inputs$estAngleMean[[i]]){
        p$parNames[[i]] <- c("mean",p$parNames[[i]])
      }
    if(DMind[[i]]){
      mle[[i]]<-matrix(mle[[i]][,1],nrow=length(p$parNames[[i]]),ncol=nbStates,byrow=TRUE)
      rownames(mle[[i]]) <- p$parNames[[i]]
      colnames(mle[[i]]) <- stateNames
    } else {
      mle[[i]]<-matrix(wpar[parindex[[i]]+1:ncol(fullDM[[i]])],1)
      rownames(mle[[i]])<-"[1,]"
      colnames(mle[[i]])<-colnames(fullDM[[i]])
      #if(is.null(names(mle[[i]]))) warning("No names for the regression coeffs were provided in DM$",i)
    }
  }

  if(!is.null(mle$beta)) {
    rownames(mle$beta) <- colnames(covs)#c("(Intercept)",attr(terms(formula),"term.labels"))
    columns <- NULL
    for(i in 1:nbStates)
      for(j in 1:nbStates) {
        if(i<j)
          columns[(i-1)*nbStates+j-i] <- paste(i,"->",j)
        if(j<i)
            columns[(i-1)*(nbStates-1)+j] <- paste(i,"->",j)
      }
    colnames(mle$beta) <- columns
  }

  # compute stationary distribution
  if(stationary) {
    gamma <- trMatrix_rcpp(nbStates,mle$beta,covs)[,,1]

    # error if singular system
    tryCatch(
      mle$delta <- solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)),
      error = function(e) {
        stop(paste("A problem occurred in the calculation of the stationary",
                   "distribution. You may want to try different initial values",
                   "and/or the option stationary=FALSE."))
      }
    )
  }

  if(nbStates==1)
    mle$delta <- 1
  names(mle$delta) <- stateNames
  
  # compute t.p.m. if no covariates
  if(nbCovs==0 & nbStates>1) {
    trMat <- trMatrix_rcpp(nbStates,mle$beta,covs)
    mle$gamma <- trMat[,,1]
    colnames(mle$gamma)<-stateNames
    rownames(mle$gamma)<-stateNames
  }

  # conditions of the fit
  conditions <- list(dist=dist,zeroInflation=zeroInflation,oneInflation=oneInflation,
                     estAngleMean=inputs$estAngleMean,circularAngleMean=inputs$circularAngleMean,stationary=stationary,formula=formula,cons=DMinputs$cons,userBounds=userBounds,bounds=p$bounds,Bndind=p$Bndind,DM=DM,fullDM=fullDM,DMind=DMind,workcons=DMinputs$workcons,fixPar=fixPar,wparIndex=wparIndex)

  mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames,knownStates=knownStates)
  
  #compute SEs and CIs on natural and working scale
  CIreal<-tryCatch(CIreal(momentuHMM(mh)),error=function(e) e)
  CIbeta<-tryCatch(CIbeta(momentuHMM(mh)),error=function(e) e)
  
  mh <- list(data=data,mle=mle,CIreal=CIreal,CIbeta=CIbeta,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames,knownStates=knownStates)
  #mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames,knownStates=knownStates)
  
  if(fit) message("DONE")
  
  return(momentuHMM(mh))
}
