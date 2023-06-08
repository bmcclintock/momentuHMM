
#' Fit a continuous-time multivariate HMM to the data
#'
#' Fit an approximate continuous-time (multivariate) hidden Markov model to the data provided, using numerical optimization of the log-likelihood
#' function. The discrete-time approximation of the continuous-time model improves as the time between observations decreases. Note that any time-varying covariates are assumed piece-wise constant between observations.
#' 
#' @param data A \code{\link{momentuHMMData}} (as returned by \code{\link{prepData}}, \code{\link{simData}}, \code{\link{prepCTDS}}, or \code{\link{simCTDS}}) or a \code{\link{momentuHierHMMData}} (as returned by \code{\link{prepData}} or \code{\link{simHierData}}) object.
#' @param ... further arguments passed to or from other methods
#' @export
fitCTHMM <- function(data, ...) {
  UseMethod("fitCTHMM")
}

#' @rdname fitCTHMM
#' @method fitCTHMM momentuHMMData
#' @param Time.name Character string indicating name of the time column. Default: 'time'. Ignored if \code{data} is a \code{ctds} object returned by \code{\link{prepCTDS}}.
#' @param Time.unit Character string indicating units for time difference between observations (e.g. 'auto', 'secs', 'mins', 'hours', 'days', 'weeks'). Ignored unless \code{data[[Time.name]]} is of class \code{\link[base]{date-time}} or \code{\link[base]{date}}. Default: 'auto'.  Ignored if \code{data} is a \code{ctds} object returned by \code{\link{prepCTDS}}.
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'bern', 'beta', 'cat', 'ctds', 'exp', 'gamma', 'lnorm', 'logis', 'negbinom', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 't', and 'weibull'. See \code{\link{fitHMM}}.
#' 
#' For continuous-time HMMs, the (multivariate) normal random walk ('rw_norm', 'rw_mvnorm2', 'rw_mvnorm3') and Poisson (`pois`) distributions are modeled as a function of the time interval between observations \eqn{(\Delta_t)}. All
#' other data stream distributions assume observations do not depend on \eqn{\Delta_t}, i.e., they
#' are ``instantaneous'' and only depend on the state active at time \eqn{t}. See details.
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
#' @param delta0 Initial value for the initial distribution of the HMM. Default: \code{rep(1/nbStates,nbStates)}. If \code{formulaDelta} includes a formula, then \code{delta0} must be specified
#' as a k x (\code{nbStates}-1) matrix, where k is the number of covariates and the columns correspond to states 2:\code{nbStates}. See details below.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). For example, \code{estAngleMean=list(angle=TRUE)} indicates the angle mean is to be 
#' estimated for 'angle'.  Default is \code{NULL}, which assumes any angle means are fixed to zero and are not to be estimated. 
#' Any \code{estAngleMean} elements corresponding to data streams that do not have angular distributions are ignored.
#' \code{estAngleMean} is also ignored for any 'vmConsensus' data streams (because the angle mean must be estimated in consensus models).
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  For example, 
#' \code{circularAngleMean=list(angle=TRUE)} indicates the angle mean is be estimated for 'angle' using circular-circular 
#' regression.  Whenever circular-circular regression is used for an angular data stream, a corresponding design matrix (\code{DM}) 
#' must be specified for the data stream, and the previous movement direction (i.e., a turning angle of zero) is automatically used 
#' as the reference angle (i.e., the intercept). Any circular-circular regression covariates in \code{data} should therefore be relative to the previous 
#' direction of movement (instead of standard directions relative to the x-axis; see \code{\link{prepData}} and \code{\link{circAngles}}).  See Duchesne et al. (2015) for specifics on the circular-circular regression model 
#' using previous movement direction as the reference angle. Default is \code{NULL}, which assumes circular-linear regression is 
#' used for any angular distributions for which the mean angle is to be estimated. \code{circularAngleMean} elements corresponding to angular data 
#' streams are ignored unless the corresponding element of \code{estAngleMean} is \code{TRUE}. Any \code{circularAngleMean} elements 
#' corresponding to data streams that do not have angular distributions are ignored. \code{circularAngleMean} is also ignored for any 'vmConsensus' data streams (because the consensus model is a circular-circular regression model).
#' 
#' Alternatively, \code{circularAngleMean} can be specified as a numeric scalar, where the value specifies the coefficient for the reference angle (i.e., directional persistence) term in the circular-circular regression model. For example, setting \code{circularAngleMean} to \code{0} specifies a 
#' circular-circular regression model with no directional persistence term (thus specifying a biased random walk instead of a biased correlated random walk). Setting \code{circularAngleMean} to 1 is equivalent to setting it to TRUE, i.e., a circular-circular regression model with a coefficient of 1 for the directional persistence reference angle.
#' @param formula Regression formula for the transition probability covariates. Default: \code{~1} (no covariate effect). In addition to allowing standard functions in R formulas
#' (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}), special functions include \code{cosinor(cov,period)} for modeling cyclical patterns, spline functions 
#' (\code{\link[splines]{bs}}, \code{\link[splines]{ns}}, \code{\link[splines2]{bSpline}}, \code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}, and \code{\link[splines2]{mSpline}}),
#'  and state- or parameter-specific formulas (see details).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities.
#' @param formulaDelta Regression formula for the initial distribution. Default for \code{fitCTHMM.momentuHMMData}: \code{NULL} (no covariate effects; both \code{delta0} and \code{fixPar$delta} are specified on the real scale). 
#' Default for \code{fitCTHMM.momentuHierHMMData}: \code{~1} (both \code{delta0} and \code{fixPar$delta} are specified on the working scale).
#' Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{delta0} and \code{fixPar$delta} are specified on the working scale.
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. Default: \code{NULL} (no covariate effects; both \code{beta0$pi} and \code{fixPar$pi} are specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{beta0$pi} and \code{fixPar$pi} are specified on the working scale.
#' Note that only the covariate values from the first row for each individual ID in \code{data} are used (i.e. time-varying covariates cannot be used for the mixture probabilties).
#' @param nlmPar List of parameters to pass to the optimization function \code{\link[stats]{nlm}} (which should be either
#' \code{print.level}, \code{gradtol}, \code{stepmax}, \code{steptol}, \code{iterlim}, or \code{hessian} -- see \code{nlm}'s documentation
#' for more detail). For \code{print.level}, the default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed. Ignored unless \code{optMethod="nlm"}.
#' @param fit \code{TRUE} if an HMM should be fitted to the data, \code{FALSE} otherwise.
#' If fit=\code{FALSE}, a model is returned with the MLE replaced by the initial parameters given in
#' input. This option can be used to assess the initial parameters, parameter bounds, etc. Default: \code{TRUE}.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a ``pseudo'' design matrix.  For example, for a 2-state 
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
#' distribution parameters for each data stream. For each matrix, the first column pertains to the lower bound and the second column the upper bound. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2,dimnames=list(c("mean_1",
#' "mean_2","concentration_1","concentration_2"))))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. For each matrix, the first column pertains to the lower bound and the second column the upper bound.
#' For data streams, each element of \code{workBounds} should be a k x 2 matrix with the same name of the corresponding element of 
#' \code{Par0}, where k is the number of parameters. For transition probability parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``beta'', where k=\code{length(beta0)}. For initial distribution parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``delta'', where k=\code{length(delta0)}.
#' \code{workBounds} is ignored for any given data stream unless \code{DM} is also specified.
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters. See details.
#' @param deltaCons Matrix of the same dimension as \code{delta0} composed of integers identifying any equality constraints among the initial distribution working scale parameters. Ignored unless a formula is provided in \code{formulaDelta}. See details.
#' @param mvnCoords Character string indicating the name of location data that are to be modeled using a multivariate normal distribution. For example, if \code{mu="mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are location data, 
#' then \code{mvnCoords="mu"} needs to be specified in order for these data to be properly treated as locations in functions such as \code{\link{plot.momentuHMM}}, \code{\link{plot.miSum}}, \code{\link{plot.miHMM}}, \code{\link{plotSpatialCov}}, and \code{\link{MIpool}}.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any). Default: NULL (states are not known). This should be a vector with length the number
#' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
#' the state is not known.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. Default: NULL 
#' (no parameters are fixed). For data streams, each element of \code{fixPar} should be a vector of the same name and length as the corresponding element of \code{Par0}. 
#' For transition probability parameters, the corresponding element of \code{fixPar} must be named ``beta'' and have the same dimensions as \code{beta0}. 
#' For initial distribution parameters, the corresponding element of \code{fixPar} must be named ``delta'' and have the same dimensions as \code{delta0}. 
#' Each parameter should either be numeric (the fixed value of the parameter) or NA if the parameter is to be estimated. Corresponding \code{fixPar} parameters must be on the same scale as \code{Par0} (e.g. if \code{DM} is specified for a given data stream, any fixed parameters for this data stream must be on the working scale), \code{beta0}, and \code{delta0}.
#'
#' If \code{optMethod='TMB'} then \code{fixPar} is specified differently. In this case, \code{NA} instead indicates parameters that are fixed (to the value specified by \code{Par0}, \code{beta0}, or \code{delta0}) and the parameters that are not fixed are sequentially indexed (similar to \code{betaCons}) by parameter number (such that parameters can be constrained to be equal). For example, \code{fixPar=list(angle=c(NA,NA,1,2))} indicates the first two parameters of the 'angle' distribution are fixed and the other two parameters are freely estimated. \code{fixPar=list(angle=c(1,NA,2,2,3))} indicates the first angle parameter is freely estimated, the second parameter is fixed, the third and fourth parameters are constrained to be equal, and the fifth parameter is freely estimated.
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization. Normal(0,\code{retrySD}^2) perturbations are used on the working scale parameters. Default: 0.  When \code{retryFits>0}, the model with the largest log likelihood 
#' value is returned. Ignored if \code{fit=FALSE}.
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization. Normal(0,\code{retrySD}^2) perturbations are used on the working scale parameters. Default: 0.  When \code{retryFits>0}, the model with the largest log likelihood 
#' value is returned. Ignored if \code{fit=FALSE}.
#' @param retrySD An optional list of scalars or vectors indicating the standard deviation to use for normal perturbations of each working scale parameter when \code{retryFits>0}. For data streams, each element of \code{retrySD} should be a vector of the same name and length as the corresponding element of \code{Par0} (if a scalar is provided, then this value will be used for all working parameters of the data stream). 
#' For transition probability parameters, the corresponding element of \code{retrySD} must be named ``beta'' and have the same dimensions as \code{beta0}. 
#' For initial distribution parameters, the corresponding element of \code{retrySD} must be named ``delta'' and have the same dimensions as \code{delta0} (if \code{delta0} is on the working scale) or be of length \code{nbStates-1} (if \code{delta0} is on the natural scale).
#' Alternatively \code{retrySD} can be a scalar, in which case this value is used for all parameters. Instead of numeric scalars, \code{retrySD} can also be specified as \code{"adapt"} in which case the standard deviation is adapted as \code{10^(ceiling(log10(abs(x))))} (where \code{x} is the current value of the corresponding working parameter)
#' Default: NULL (in which case \code{retrySD}=1 for data stream parameters and \code{retrySD}=10 for initial distribution and state transition probabilities). Ignored unless \code{retryFits>0}.
#' @param ncores Number of cores to use for parallel processing when \code{retryFits>0}. Default: 1 (no parallel processing). When \code{retryFits} attempts are performed in parallel, each attempt uses perturbations of the initial model fit (i.e. they are not iteratively updated as when \code{ncores=1}).
#' @param optMethod The optimization method to be used.  Can be ``nlm'' (the default; see \code{\link[stats]{nlm}}), ``TMB'' (using Template Model Builder; see \code{\link[optimx]{optimx}} for control parameters), ``Nelder-Mead'' (see \code{\link[stats]{optim}}), or ``SANN'' (see \code{\link[stats]{optim}}).
#' @param control A list of control parameters to be passed to \code{\link[stats]{optim}} (ignored unless \code{optMethod="TMB"}, \code{optMethod="Nelder-Mead"}, or \code{optMethod="SANN"}).
#' @param prior A function that returns the log-density of the working scale parameter prior distribution(s). See 'Details'.
#' @param modelName An optional character string providing a name for the fitted model. If provided, \code{modelName} will be returned in \code{\link{print.momentuHMM}}, \code{\link{AIC.momentuHMM}}, \code{\link{AICweights}}, and other functions. 
#' @param kappa maximum allowed value for the row sums of the off-diagonal elements in the state transition rate matrix, such that the minimum value for the diagonal elements is \code{-kappa}. Default: \code{Inf}. Setting less than \code{Inf} can help avoid numerical issues during optimization, in which case the transition rate parameters \code{beta} are on the logit scale (instead of the log scale).
#'
#' @return A \code{\link{momentuHMM}} or \code{\link{momentuHierHMM}} object, i.e. a list of:
#' \item{mle}{A named list of the maximum likelihood estimates of the parameters of the model (if the numerical algorithm
#' has indeed identified the global maximum of the likelihood function). Elements are included for the parameters of each
#' data strea, as well as \code{beta} (transition probabilities regression coefficients - more information
#' in 'Details'), \code{gamma} (transition probabilities on real scale, based on mean covariate values if \code{formula}
#' includes covariates), and \code{delta} (initial distribution).}
#' \item{CIreal}{Standard errors and 95\% confidence intervals on the real (i.e., natural) scale of parameters}
#' \item{CIbeta}{Standard errors and 95\% confidence intervals on the beta (i.e., working) scale of parameters}
#' \item{data}{The momentuHMMData or momentuHierHMMData object}
#' \item{mod}{List object returned by the numerical optimizer \code{nlm} or \code{optim}. Items in \code{mod} include the best set of free working parameters found (\code{wpar}), 
#' the best full set of working parameters including any fixed parameters (\code{estimate}), the value of the likelihood at \code{estimate} (\code{minimum}), 
#' the estimated variance-covariance matrix at \code{estimate} (\code{Sigma}), and the elapsed time in seconds for the optimization (\code{elapsedTime}).}
#' \item{conditions}{Conditions used to fit the model, e.g., \code{bounds} (parameter bounds), distributions, \code{zeroInflation},
#' \code{estAngleMean}, \code{stationary}, \code{formula}, \code{DM}, \code{fullDM} (full design matrix), etc.}
#' \item{rawCovs}{Raw covariate values for transition probabilities, as found in the data (if any). Used in \code{\link{plot.momentuHMM}}.}
#' \item{stateNames}{The names of the states.}
#' \item{knownStates}{Vector of values of the state process which are known.}
#' \item{covsDelta}{Design matrix for initial distribution.}
#' 
#' @details 
#' \code{fitCTHMM} assumes the snapshot property applies to all data stream distributions (i.e. observations are ``instantaneous'') except for the (multivariate) normal random walk (\code{rw_norm}, \code{rw_mvnorm2}, \code{rw_mvnorm3}) and Poisson (\code{pois}) distributions. For these particular distributions, the observed data are not ``instantaneous''; they depend on the time interval between observations \eqn{(\Delta_t)} and, hence, the state sequence during the entire interval.
#' When fitting with \code{\link{fitCTHMM}} (or \code{\link{MIfitCTHMM}}), it is critical that the frequency of observations is high relative to the serial correlation in the hidden state process in order for the discrete-time approximation of \code{\link{fitCTHMM}} to be reasonably accurate for these distributions.
#' 
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats model.matrix get_all_vars nlm optim terms terms.formula
#' @importFrom CircStats dwrpcauchy dvm pvm
#' @importFrom MASS ginv
#'
#' @useDynLib momentuHMM

fitCTHMM.momentuHMMData <- function(data,Time.name="time",Time.unit="auto",nbStates,dist,
                                  Par0,beta0=NULL,delta0=NULL,
                                  estAngleMean=NULL,circularAngleMean=NULL,
                                  formula=~1,formulaDelta=NULL,stationary=FALSE,mixtures=1,formulaPi=NULL,
                                  nlmPar=list(),fit=TRUE,
                                  DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,deltaCons=NULL,
                                  mvnCoords=NULL,stateNames=NULL,knownStates=NULL,fixPar=NULL,retryFits=0,retrySD=NULL,ncores=1,optMethod="nlm",control=list(),prior=NULL,modelName=NULL, kappa=Inf, ...)
{
  
  #####################
  ## Check arguments ##
  #####################
  
  # check that the data is a momentuHMMData object
  if(!is.momentuHMMData(data))
    stop("'data' must be a momentuHMMData object (as output by prepData or simData)")
  
  if(length(data)<1 | any(dim(data)<1))
    stop("The data input is empty.")
  
  if("dt" %in% names(data)) stop("'dt' is reserved and cannot be used as a field name in 'data'; please choose a different name") 
  data$dt <- 0
  
  if(inherits(data,"ctds")){
    Time.name <- "time"
    Time.unit <- attr(data,"Time.unit")
    if(is.null(data$tau)) stop("'tau' not found in 'data'")
    data$dt <- data$tau
    data$tau <- NULL
    if(!any(data[[attr(data,"ctdsData")]]==(attr(data,"directions")+1))) stop("cell moves cannot occur at every time step; some data$",attr(data,"ctdsData"),"'s must be equal to ",attr(data,"directions")+1)
  } else {
    if(!Time.name %in% names(data)) stop("Time.name not found in 'data'") 
    for(i in unique(data$ID)){
      if(inherits(data,"hierarchical")){
        for(iLevel in levels(data$level)){
          iDat <- which(data$ID==i & data$level==iLevel)
          if(!grepl("i",iLevel)){
            if(inherits(data[[Time.name]],"POSIXt")) {
              tmpTime <- difftime(data[[Time.name]][iDat][-1],data[[Time.name]][iDat][-length(iDat)],units=Time.unit)
              Time.unit <- units(tmpTime)
              data$dt[iDat] <- c(tmpTime,0)
            } else data$dt[iDat] <- c(diff(data[[Time.name]][iDat]),0)
          } else data$dt[iDat] <- 1
        }
      } else {
        iDat <- which(data$ID==i)
        if(inherits(data[[Time.name]],"POSIXt")) {
          tmpTime <- difftime(data[[Time.name]][iDat][-1],data[[Time.name]][iDat][-length(iDat)],units=Time.unit)
          Time.unit <- units(tmpTime)
          data$dt[iDat] <- c(tmpTime,0)
        } else data$dt[iDat] <- c(diff(data[[Time.name]][iDat]),0)
      }
    }
  }
  
  for(i in names(dist)){
    if(!dist[[i]] %in% CTHMMdists) stop("Sorry, currently fitCTHMM only supports the following distributions: ",paste0(CTHMMdists,sep=", "))
    if(dist[[i]]=="ctds" & !inherits(data,"ctds")) stop("'ctds' models require data to be a ctds object (see ?prepCTDS")
  }
  
  if(!is.numeric(kappa) || length(kappa)>1 || kappa<0) stop('kappa must be a non-negative numeric scalar')
  
  chkDots(...)
  
  withCallingHandlers(mfit <- fitHMM.momentuHMMData(data=data,nbStates=nbStates,dist=dist,
                               Par0=Par0,beta0=beta0,delta0=delta0,
                               estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                               formula=formula,formulaDelta=formulaDelta,stationary=stationary,mixtures=mixtures,formulaPi=formulaPi,
                               nlmPar=nlmPar,fit=fit,
                               DM=DM,userBounds=userBounds,workBounds=workBounds,betaCons=betaCons,betaRef=NULL,deltaCons=deltaCons,
                               mvnCoords=mvnCoords,stateNames=stateNames,knownStates=knownStates,fixPar=fixPar,retryFits=retryFits,retrySD=retrySD,ncores=ncores,optMethod=optMethod,control=control,prior=prior,modelName=modelName, CT=TRUE, Time.name=Time.name, kappa=kappa),warning=muffleCTwarning)
  
  class(mfit) <- append(class(mfit),"CTHMM")
  mfit$conditions$CT <- TRUE
  attr(mfit$data,"CT") <- TRUE
  attr(mfit$data,"Time.name") <- Time.name
  attr(mfit$data,"Time.unit") <- Time.unit
  if(inherits(data,"ctds")){
    class(mfit) <- append(class(mfit),"ctds")
    attr(mfit$data,"directions") <- attr(data,"directions")
    attr(mfit$data,"ctdsData") <- attr(data,"ctdsData")
    attr(mfit$data,"normalize.gradients") <- attr(data,"normalize.gradients")
    attr(mfit$data,"grad.point.decreasing") <- attr(data,"grad.point.decreasing")
  }

  mfit$CIreal <- tryCatch(CIreal(mfit),error=function(e) e)
  if(inherits(mfit$CIreal,"error") & fit==TRUE) warning("Failed to compute SEs and confidence intervals on the natural scale -- ",mfit$CIreal)
  
  if(fit) {
    # check for numerical underflow in state-dependent observation distributions
    probs <- tryCatch(allProbs(mfit),warning=function(w) w)
    if(inherits(probs,"warning")) warning(probs)
    
    # check for numerical overflow in state transition rates
    trProbs <- tryCatch(getTrProbs(momentuHMM(mfit)),error=function(e) e)
  }
  mfit$conditions$Time.name <- Time.name
  mfit$conditions$Time.unit <- Time.unit
  return(mfit)
  

}

#' @rdname fitCTHMM
#' @method fitCTHMM momentuHierHMMData
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See details.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See details.
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See details.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (only hierarchical-level effects, with no covariate effects).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy. See details.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' 
#' @details
#' \code{fitCTHMM.momentuHierHMMData} is very similar to \code{\link{fitCTHMM.momentuHMMData}} except that instead of simply specifying the number of states (\code{nbStates}), distributions (\code{dist}), and a single t.p.m. formula (\code{formula}), the \code{hierStates} argument specifies the hierarchical nature of the states,
#' the \code{hierDist} argument specifies the hierarchical nature of the data streams, and the \code{hierFormula} argument specifies a t.p.m. formula for each level of the hierarchy.  All are specified as 
#' \code{\link[data.tree]{Node}} objects from the \code{\link[data.tree]{data.tree}} package.
#' 
#' @seealso \code{\link{simHierData}}
#' 
#' @export
fitCTHMM.momentuHierHMMData <- function(data,Time.name="time",Time.unit="auto",hierStates,hierDist,
                                      Par0,hierBeta=NULL,hierDelta=NULL,
                                      estAngleMean=NULL,circularAngleMean=NULL,
                                      hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,
                                      nlmPar=list(),fit=TRUE,
                                      DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,deltaCons=NULL,
                                      mvnCoords=NULL,knownStates=NULL,fixPar=NULL,retryFits=0,retrySD=NULL,ncores=1,optMethod="nlm",control=list(),prior=NULL,modelName=NULL, kappa=Inf, ...)
{
  
  if(!inherits(data,"momentuHierHMMData")) stop("data must be a momentuHierHMMData object (as returned by prepData or simHierData)")
  
  if(length(data)<1 | any(dim(data)<1))
    stop("The data input is empty.")
  
  inputHierHMM <- formatHierHMM(data,hierStates,hierDist,hierBeta,hierDelta,hierFormula,hierFormulaDelta,mixtures,workBounds,betaCons,deltaCons,fixPar,checkData = TRUE, CT=TRUE)

  chkDots(...)
  
  withCallingHandlers(hfit <- fitCTHMM.momentuHMMData(momentuHMMData(data),Time.name=Time.name,Time.unit=Time.unit,inputHierHMM$nbStates,inputHierHMM$dist,Par0,inputHierHMM$beta,inputHierHMM$delta,
                                estAngleMean,circularAngleMean,
                                formula=inputHierHMM$formula,inputHierHMM$formulaDelta,stationary=FALSE,mixtures,formulaPi,
                                nlmPar,fit,
                                DM,userBounds,workBounds=inputHierHMM$workBounds,inputHierHMM$betaCons,inputHierHMM$betaRef,deltaCons=inputHierHMM$deltaCons,
                                mvnCoords,inputHierHMM$stateNames,knownStates,inputHierHMM$fixPar,retryFits,retrySD,ncores,optMethod,control,prior,modelName, CT=TRUE, kappa=kappa),warning=muffleCTwarning)
  
  # replace initial values with estimates in hierBeta and hierDelta (if provided)
  if(fit){
    par <- getPar(hfit)
    if(is.list(par$beta)){
      beta <- par$beta$beta
      Pi <- par$beta$pi
      g0 <- par$beta$g0
      names(g0) <- names(hfit$mle$g0)
      theta <- par$beta$theta
      names(theta) <- names(hfit$mle$theta)
    } else {
      beta <- par$beta
      Pi <- g0 <- theta <- NULL
    }
    hier <- mapHier(list(beta=beta,g0=g0,theta=theta),Pi,par$delta,hierBeta,hierDelta,inputHierHMM$hFixPar,inputHierHMM$hBetaCons,inputHierHMM$hDeltaCons,hierStates,inputHierHMM$newformula,hfit$conditions$formulaDelta,inputHierHMM$data,hfit$conditions$mixtures,inputHierHMM$recharge,fill=TRUE)
    if(!is.null(inputHierHMM$recharge)){
      for(j in names(inputHierHMM$recharge)){
        hier$hierBeta[[j]]$g0 <- g0[names(hier$hierBeta[[j]]$g0)]
        hier$hierBeta[[j]]$theta <- theta[names(hier$hierBeta[[j]]$theta)]
      }
    }
    hfit$conditions$hierBeta <- hier$hierBeta
    hfit$conditions$hierDelta <- hier$hierDelta
  } else {
    hfit$conditions$hierBeta <- hierBeta
    hfit$conditions$hierDelta <- hierDelta    
  }
  hfit$conditions$hierStates <- hierStates
  hfit$conditions$hierDist <- hierDist
  hfit$conditions$hierFormula <- hierFormula
  hfit$conditions$hierFormulaDelta <- hierFormulaDelta
  class(hfit$data) <- class(data)
  
  hfit <- momentuHierHMM(hfit)
  
  if(fit){
    hfit$CIreal <- tryCatch(CIreal.hierarchical(hfit),error=function(e) e)
    if(inherits(hfit$CIreal,"error") & fit==TRUE) warning("Failed to compute SEs and confidence intervals on the natural scale -- ",hfit$CIreal)
    trProbs <- tryCatch(getTrProbs(momentuHMM(hfit)),error=function(e) e)
  }
  return(hfit)
}
