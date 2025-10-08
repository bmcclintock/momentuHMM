
#' Fit a multivariate HMM to the data
#'
#' Fit a (multivariate) hidden Markov model to the data provided, using numerical optimization of the log-likelihood
#' function.
#' 
#' @param data A \code{\link{momentuHMMData}} (as returned by \code{\link{prepData}} or \code{\link{simData}}) or a \code{\link{momentuHierHMMData}} (as returned by \code{\link{prepData}} or \code{\link{simHierData}}) object.
#' @param ... further arguments passed to or from other methods
#' @export
fitHMM <- function(data, ...) {
  UseMethod("fitHMM")
}

#' @rdname fitHMM
#' @method fitHMM momentuHMMData
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'bern', 'beta', 'cat', exp', 'gamma', 'lnorm', 'logis', 'negbinom', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. For example,
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
#' @param formulaDelta Regression formula for the initial distribution. Default for \code{fitHMM.momentuHMMData}: \code{NULL} (no covariate effects; both \code{delta0} and \code{fixPar$delta} are specified on the real scale). 
#' Default for \code{fitHMM.momentuHierHMMData}: \code{~1} (both \code{delta0} and \code{fixPar$delta} are specified on the working scale).
#' Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{delta0} and \code{fixPar$delta} are specified on the working scale.
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. Default: \code{FALSE}.
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. Default: \code{NULL} (no covariate effects; both \code{beta0$pi} and \code{fixPar$pi} are specified on the real scale). Standard functions in R formulas are allowed (e.g., \code{cos(cov)}, \code{cov1*cov2}, \code{I(cov^2)}). When any formula is provided, then both \code{beta0$pi} and \code{fixPar$pi} are specified on the working scale.
#' Note that only the covariate values from the first row for each individual ID in \code{data} are used (i.e. time-varying covariates cannot be used for the mixture probabilities).
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
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the t.p.m. multinomial logit link. Default: NULL, in which case
#' the diagonal elements of the t.p.m. are the reference. See details.
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
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization. Normal(0,\code{retrySD}^2) perturbations are used on the working scale parameters. Default: 0.  When \code{retryFits>0}, the model with the largest log likelihood 
#' value is returned. Ignored if \code{fit=FALSE}.
#' @param retrySD An optional list of scalars or vectors indicating the standard deviation to use for normal perturbations of each working scale parameter when \code{retryFits>0}. For data streams, each element of \code{retrySD} should be a vector of the same name and length as the corresponding element of \code{Par0} (if a scalar is provided, then this value will be used for all working parameters of the data stream). 
#' For transition probability parameters, the corresponding element of \code{retrySD} must be named ``beta'' and have the same dimensions as \code{beta0}. 
#' For initial distribution parameters, the corresponding element of \code{retrySD} must be named ``delta'' and have the same dimensions as \code{delta0} (if \code{delta0} is on the working scale) or be of length \code{nbStates-1} (if \code{delta0} is on the natural scale).
#' Alternatively \code{retrySD} can be a scalar, in which case this value is used for all parameters.
#' Default: NULL (in which case \code{retrySD}=1 for data stream parameters and \code{retrySD}=10 for initial distribution and state transition probabilities). Ignored unless \code{retryFits>0}.
#' @param optMethod The optimization method to be used.  Can be ``nlm'' (the default; see \code{\link[stats]{nlm}}), ``Nelder-Mead'' (see \code{\link[stats]{optim}}), or ``SANN'' (see \code{\link[stats]{optim}}).
#' @param control A list of control parameters to be passed to \code{\link[stats]{optim}} (ignored unless \code{optMethod="Nelder-Mead"} or \code{optMethod="SANN"}).
#' @param prior A function that returns the log-density of the working scale parameter prior distribution(s). See 'Details'.
#' @param modelName An optional character string providing a name for the fitted model. If provided, \code{modelName} will be returned in \code{\link{print.momentuHMM}}, \code{\link{AIC.momentuHMM}}, \code{\link{AICweights}}, and other functions. 
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
#' \itemize{
#' \item By default the matrix \code{beta0} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 \code{formula} covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix -
#' filled in row-wise).
#' In a covariate-free model (default), \code{beta0} has one row, for the intercept. While the diagonal elements are by default the reference elements, other elements can serve as the reference
#' using the \code{betaRef} argument. For example, in a 3-state model, setting \code{betaRef=c(3,2,3)} changes the reference elements to state transition 1 -> 3 for state 1 (instead of 1 -> 1), state
#' transition 2 -> 2 for state 2 (same as default), and state transition 3 -> 3 for state 3 (same as default).
#' 
#' \item When covariates are not included in \code{formulaDelta} (i.e. \code{formulaDelta=NULL}), then \code{delta0} (and \code{fixPar$delta}) are specified as a vector of length \code{nbStates} that 
#' sums to 1.  When any formula is specified for \code{formulaDelta} (e.g. \code{formulaDelta=~1}, \code{formulaDelta=~cov1}), then \code{delta0}  (and \code{fixPar$delta}) must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta0} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta0)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#'
#' \item The choice of initial parameters (particularly \code{Par0} and \code{beta0}) is crucial to fit a model. The algorithm might not find
#' the global optimum of the likelihood function if the initial parameters are poorly chosen.
#' 
#' \item If \code{DM} is specified for a particular data stream, then the initial values are specified on 
#' the working (i.e., beta) scale of the parameters. The working scale of each parameter is determined by the link function used.
#' If a parameter P is bound by (0,Inf) then the working scale is the log(P) scale.  If the parameter bounds are (-pi,pi) then the working 
#' scale is tan(P/2) unless circular-circular regression is used. Otherwise if the parameter bounds are finite then logit(P) is the working scale. However, when both 
#' zero- and one-inflation are included, then a multinomial logit link is used because the sum of the zeromass and onemass probability parameters cannot exceed 1.
#' The function \code{\link{getParDM}} is intended to help with obtaining initial values on the working scale when specifying a design matrix and other 
#' parameter constraints (see example below). When circular-circular regression is specified using \code{circularAngleMean}, the working scale 
#' for the mean turning angle is not as easily interpretable, but the 
#' link function is atan2(sin(X)*B,1+cos(X)*B), where X are the angle covariates and B the angle coefficients (see Duchesne et al. 2015). 
#' Under this formulation, the reference turning angle is 0 (i.e., movement in the same direction as the previous time step). 
#' In other words, the mean turning angle is zero when the coefficient(s) B=0.
#' 
#' \item Circular-circular regression in \code{momentuHMM} is designed for turning angles (not bearings) as computed by \code{\link{simData}} and \code{\link{prepData}}. 
#' Any circular-circular regression angle covariates for time step t should therefore be relative to the previous 
#' direction of movement for time step t-1.  In other words, circular-circular regression covariates for time step t should be the turning angle
#' between the direction of movement for time step t-1 and the standard direction of the covariate relative to the x-axis for time step t.  If provided standard directions in radians relative to the x-axis 
#' (where 0 = east, pi/2 = north, pi = west, and -pi/2 = south), \code{\link{circAngles}} or \code{\link{prepData}} can perform this calculation for you.  
#'
#' When the circular-circular regression model is used, the special function \code{angleFormula(cov,strength,by)} can be used in \code{DM} for the mean of angular distributions (i.e. 'vm', 'vmConsensus', and 'wrpcauchy'),
#' where \code{cov} is an angle covariate (e.g. wind direction), \code{strength} is an optional positive real covariate (e.g. wind speed), and \code{by} is an optional factor variable for individual- or group-level effects (e.g. ID, sex). The \code{strength} argument allows angle covariates to be weighted based on their relative strength or importance at time step t as in
#' Rivest et al. (2016).  In this case, the link function for the mean angle is atan2((Z * sin(X)) \%*\% B,1+(Z * cos(X)) \%*\% B), where X are the angle covariates, Z the strength covariates, and B the angle coefficients (see Rivest et al. 2016). 
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
#' \item \code{betaCons} can be used to impose equality constraints among the t.p.m. parameters.  It must be a matrix of the same dimension as \code{beta0} and be composed of integers, where each beta parameter is sequentially indexed in a column-wise fashion (see \code{\link{checkPar0}}). Parameter indices in \code{betaCons} must therefore be integers between \code{1} and \code{nbStates*(nbStates-1)}. 
#' 
#' Use of \code{betaCons} is perhaps best demonstrated by example.  If no constraints are imposed (the default), then \code{betaCons=matrix(1:length(beta0),nrow(beta0),ncol(beta0))} such that
#' each beta parameter is (column-wise) sequentially identified by a unique integer.  Suppose we wish to fit a model with \code{nbStates=3} states and a covariate (`cov1') on the t.p.m. With no constraints on the t.p.m., we would have
#' \code{betaCons=matrix(1:(2*(nbStates*(nbStates-1))),nrow=2,ncol=nbStates*(nbStates-1),dimnames=list(c("(Intercept)","cov1"),c("1 -> 2","1 -> 3","2 -> 1","2 -> 3","3 -> 1","3 -> 2")))}.  If we then wanted to constrain the t.p.m. such that the covariate effect is identical for transitions from state 1 to states 2 and 3 (and vice versa), we have
#' \code{betaCons=matrix(c(1,2,3,2,5,6,7,8,9,6,11,12),nrow=2,ncol=nbStates*(nbStates-1),dimnames=list(c("(Intercept)","cov1"),c("1 -> 2","1 -> 3","2 -> 1","2 -> 3","3 -> 1","3 -> 2")))}; this results in 10 estimated beta parameters (instead of 12), the ``cov1'' effects indexed by a ``2'' (``1 -> 2'' and ``1 -> 3'') constrained to be equal, and 
#' the ``cov1'' effects indexed by a ``6'' (``2 -> 1'' and ``3 -> 1'') constrained to be equal. 
#' 
#' Now suppose we instead wish to constrain these sets of state transition probabilities to be equal, i.e., Pr(1 -> 2) = Pr(1 -> 3) and Pr(2 -> 1) = Pr(3 -> 1); then we have \code{betaCons=matrix(c(1,2,1,2,5,6,7,8,5,6,11,12),nrow=2,ncol=nbStates*(nbStates-1),dimnames=list(c("(Intercept)","cov1"),c("1 -> 2","1 -> 3","2 -> 1","2 -> 3","3 -> 1","3 -> 2")))}
#' 
#' \item Cyclical relationships (e.g., hourly, monthly) may be modeled in \code{DM} or \code{formula} using the \code{cosinor(x,period)} special formula function for covariate \code{x}
#' and sine curve period of time length \code{period}. For example, if 
#' the data are hourly, a 24-hour cycle can be modeled using \code{~cosinor(cov1,24)}, where the covariate \code{cov1} is a repeating sequential series
#' of integers indicating the hour of day (\code{0,1,...,23,0,1,...,23,0,1,...}) (note that \code{fitHMM} will not do this for you, the appropriate covariate must be included in \code{data}; see example below). 
#' The \code{cosinor(x,period)} function converts \code{x} to 2 covariates \code{cosinorCos(x)=cos(2*pi*x/period)} and \code{cosinorSin(x)=sin(2*pi*x/period} for inclusion in the model (i.e., 2 additional parameters per state). The amplitude of the sine wave
#' is thus \code{sqrt(B_cos^2 + B_sin^2)}, where \code{B_cos} and \code{B_sin} are the working parameters correponding to \code{cosinorCos(x)} and \code{cosinorSin(x)}, respectively (e.g., see Cornelissen 2014).
#' 
#' \item Similar to that used in \code{\link{crawlWrap}}, the \code{prior} argument is a user-specified function that returns the log-density of the working scale parameter prior distribution(s). In addition to including prior information about parameters, one area where priors can be particularly useful is for handling numerical issues that can arise when parameters are near a boundary. 
#' When parameters are near boundaries, they can wander into the ``nether regions'' of the parameter space during optimization. For example, setting \code{prior=function(par) {sum(dnorm(par,0,sd,log=TRUE))}} with a reasonably large \code{sd} (e.g. 100 or 1000) can help prevent working parameters 
#' from straying too far along the real line.  Here \code{par} is the vector of working scale parameters (as returned by \code{fitHMM}, e.g., see \code{example$m$mod$estimate}) in the following order: data stream working parameters (in order \code{names(dist)}), beta working parameters, and delta working parameters. Instead of specifying the same prior on all parameters, different priors could be specified on different parameters (and not all parameters must have user-specified priors).  For example,
#' \code{prior=function(par){dnorm(par[3],0,100,log=TRUE)}} would only specify a prior for the third working parameter. Note that the \code{prior} function must return a scalar on the log scale. See 'harbourSealExample.R' in the ``vignettes'' source directory for an example using the \code{prior} argument.
#' 
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
#' data <- example$m$data
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
#' 
#' ### 8. Cosinor and state-dependent formulas
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
#'                       beta=beta,formula=formula,covs=covs)    
#' 
#' m.cosinor<-fitHMM(data.cos,nbStates=nbStates,dist=dist,Par0=Par,formula=formula)
#' m.cosinor    
#' 
#' ### 9. Piecewise constant B-spline on step length mean and angle concentration
#' nObs <- 1000 # length of simulated track
#' cov <- data.frame(time=1:nObs) # time covariate for splines
#' dist <- list(step="gamma",angle="vm")
#' stepDM <- list(mean=~splines2::bSpline(time,df=2,degree=0),sd=~1)
#' angleDM <- list(mean=~1,concentration=~splines2::bSpline(time,df=2,degree=0))
#' DM <- list(step=stepDM,angle=angleDM)
#' Par <- list(step=c(log(1000),1,-1,log(100)),angle=c(0,log(10),2,-5))
#' 
#' data.spline<-simData(obsPerAnimal=nObs,nbStates=1,dist=dist,Par=Par,DM=DM,covs=cov) 
#' 
#' Par0 <- list(step=Par$step,angle=Par$angle[-1])
#' m.spline<-fitHMM(data.spline,nbStates=1,dist=dist,Par0=Par0,
#'                  DM=list(step=stepDM,
#'                          angle=angleDM["concentration"]))  
#' 
#' ### 10. Initial state (delta) based on covariate                       
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
#' Par0 <- list(step=Par$step, angle=Par$angle[3:4])   
#' m.delta <- fitHMM(data.delta, nbStates=2, dist=dist, Par0 = Par0, 
#'                   formulaDelta=formulaDelta)
#'                   
#' ### 11. Two mixtures based on covariate                       
#' nObs <- 100
#' nbAnimals <- 20
#' dist <- list(step="gamma",angle="vm")
#' Par <- list(step=c(100,1000,50,100),angle=c(0,0,0.1,2))
#' 
#' # create sex covariate
#' cov <- data.frame(sex=factor(rep(c("F","M"),each=nObs*nbAnimals/2)))
#' formulaPi <- ~ sex + 0
#' 
#' # Females more likely in mixture 1, males more likely in mixture 2
#' beta <- list(beta=matrix(c(-1.5,-0.5,-1.5,-3),2,2),
#'              pi=matrix(c(-2,2),2,1,dimnames=list(c("sexF","sexM"),"mix2"))) 
#' 
#' data.mix<-simData(nbAnimals=nbAnimals,obsPerAnimal=nObs,nbStates=2,dist=dist,Par=Par,
#'                   beta=beta,formulaPi=formulaPi,mixtures=2,covs=cov) 
#' 
#' Par0 <- list(step=Par$step, angle=Par$angle[3:4])   
#' m.mix <- fitHMM(data.mix, nbStates=2, dist=dist, Par0 = Par0, 
#'                 beta0=beta,formulaPi=formulaPi,mixtures=2)
#' }
#' 
#' @references
#' 
#' Cornelissen, G. 2014. Cosinor-based rhythmometry. Theoretical Biology and Medical Modelling 11:16.
#' 
#' Duchesne, T., Fortin, D., Rivest L-P. 2015. Equivalence between step selection functions and 
#' biased correlated random walks for statistical inference on animal movement. PLoS ONE 10 (4):
#' e0122947.
#' 
#' Langrock R., King R., Matthiopoulos J., Thomas L., Fortin D., Morales J.M. 2012.
#' Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions.
#' Ecology, 93 (11), 2336-2342.
#' 
#' Leos-Barajas, V., Gangloff, E.J., Adam, T., Langrock, R., van Beest, F.M., Nabe-Nielsen, J. and Morales, J.M. 2017. 
#' Multi-scale modeling of animal movement and general behavior data using hidden Markov models with hierarchical structures. 
#' Journal of Agricultural, Biological and Environmental Statistics, 22 (3), 232-248.
#' 
#' Maruotti, A., and T. Ryden. 2009. A semiparametric approach to hidden Markov models under longitudinal 
#' observations. Statistics and Computing 19: 381-393.
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
#' Rivest, LP, Duchesne, T, Nicosia, A, Fortin, D, 2016. A general angular regression model for the analysis of data on animal movement in ecology. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 65(3):445-463.
#' 
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats model.matrix get_all_vars nlm optim terms terms.formula
#' @importFrom CircStats dwrpcauchy dvm pvm
#' @importFrom MASS ginv
#'
#' @useDynLib momentuHMM

fitHMM.momentuHMMData <- function(data,nbStates,dist,
                   Par0,beta0=NULL,delta0=NULL,
                   estAngleMean=NULL,circularAngleMean=NULL,
                   formula=~1,formulaDelta=NULL,stationary=FALSE,mixtures=1,formulaPi=NULL,
                   nlmPar=list(),fit=TRUE,
                   DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,betaRef=NULL,deltaCons=NULL,
                   mvnCoords=NULL,stateNames=NULL,knownStates=NULL,fixPar=NULL,retryFits=0,retrySD=NULL,optMethod="nlm",control=list(),prior=NULL,modelName=NULL, ...)
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
  
  if(!is.null(formulaDelta) && !is.formula(formulaDelta))
    stop("Check the argument 'formulaDelta'.")
  
  if(nbStates<1) stop('nbStates must be >0')
  
  chkDots(...)

  # check that there is no response variable in the formula
  if(attr(stats::terms(formula),"response")!=0)
    stop("The response variable should not be specified in formula.")
  if(!is.null(formulaDelta) && attr(stats::terms(formulaDelta),"response")!=0)
    stop("The response variable should not be specified in formulaDelta.")
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par0) | is.null(names(Par0))) stop("'Par0' must be a named list")
  distnames<-tmpdistnames<-names(dist)
  if(any(is.na(match(distnames,names(data))))){
    for(i in which(is.na(match(distnames,names(data))))){
      if(dist[[distnames[i]]] %in% mvndists){
        if(dist[[distnames[i]]] %in% c("mvnorm2","rw_mvnorm2")){
          tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
        } else if(dist[[distnames[i]]] %in% c("mvnorm3","rw_mvnorm3")){
          tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
        }
      }
    }
    if(any(is.na(match(tmpdistnames,names(data))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(data)))],collapse=", ")," not found in data")
  } 
  for(i in tmpdistnames){
    if(all(is.na(data[[i]]))) warning("data$",i," consists entirely of NAs")
  }
  
  if(!all(distnames %in% names(Par0))) stop(paste0(distnames[which(!(distnames %in% names(Par0)))],collapse=", ")," missing in 'Par0'")
  Par0 <- Par0[distnames]
  
  match.arg(optMethod,fitMethods)
  
  if(!is.null(mvnCoords)){
    if(length(mvnCoords)>1 | !is.character(mvnCoords)) stop("mvnCoords must be a character string")
    if(!(mvnCoords %in% distnames)) stop("mvnCoords not found. Permitted values are: ",paste0(distnames,collapse=", "))
    if(!(dist[[mvnCoords]] %in% mvndists)) stop("mvnCoords must correspond to a multivariate normal data stream")
  }
  
  # check knownStates
  if(length(knownStates) > 0){
    if(length(knownStates) != nrow(data)) 
      stop("'knownStates' should be of same length as the data, i.e. ",nrow(data))
    if(!all(is.na(knownStates))) {
      if(max(knownStates, na.rm = TRUE) > nbStates | min(knownStates, na.rm = TRUE) < 1 | !isTRUE(all.equal(knownStates,as.integer(knownStates)))) 
        stop("'knownStates' should only contain integers between 1 and ", nbStates, " (or NAs)")
    }
  }
  
  # convert RW data
  if(any(unlist(dist) %in% rwdists)){
    data <- RWdata(dist,data,knownStates)
    knownStates <- data$knownStates
    data$knownStates <- NULL
  }
  
  newForm <- newFormulas(formula,nbStates,betaRef,hierarchical = TRUE)
  formulaStates <- newForm$formulaStates
  newformula <- newForm$newformula
  recharge <- newForm$recharge
  
  # build design matrix for t.p.m.
  covsCol <- tryCatch(stats::get_all_vars(newformula,data=data),error=function(e) e)#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
  if(inherits(covsCol,"error")){
    if(any(grepl("MIfitHMM",unlist(lapply(sys.calls(),function(x) deparse(x)[1]))))){
      stop(covsCol$message,"\n     -- has MIfitHMM 'covNames' argument been correctly specified?")
    } else stop(covsCol)
  }
  if(!all(names(covsCol) %in% names(data))){
    for(j in names(covsCol)[which(!names(covsCol) %in% names(data))]){
      if(exists(j)) stop("'",j,"' covariate not found in data, but a variable named '",j,"' is present in the environment (this can cause major problems!)",
                         ifelse(any(grepl("MIfitHMM",unlist(lapply(sys.calls(),function(x) deparse(x)[1])))),
                                " \n       -- has MIfitHMM 'covNames' argument been correctly specified?",
                                ""))
    }
    covsCol <- covsCol[,names(covsCol) %in% names(data),drop=FALSE]
  }
  covs <- stats::model.matrix(newformula,data)
  if(nrow(covs)!=nrow(data)) stop("covariates cannot contain missing values")
  nbCovs <- ncol(covs)-1 # substract intercept column
  
  # aInd = list of indices of first observation for each animal
  aInd <- NULL
  nbAnimals <- length(unique(data$ID))
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])
  
  if(mixtures>1){
    if(!is.null(beta0)){
      if(!is.list(beta0)) stop("beta0 must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
  }
  
  # build design matrix for recharge model
  if(!is.null(recharge)){
    reForm <- formatRecharge(nbStates,formula,betaRef,data=data)
    formulaStates <- reForm$formulaStates
    newformula <- reForm$newformula
    recharge <- reForm$recharge
    hierRecharge <- reForm$hierRecharge
    newdata <- reForm$newdata
    g0covs <- reForm$g0covs
    nbG0covs <- ncol(g0covs)-1
    recovs <- reForm$recovs
    nbRecovs <- ncol(recovs)-1
    if(!nbRecovs | (!inherits(data,"hierarchical") && !attributes(stats::terms(recharge$theta))$intercept)) stop("invalid recharge model -- theta must include an intercept and at least 1 covariate")
    covs <- reForm$covs
    nbCovs <- ncol(covs)-1
    if(!is.null(beta0)){
      if(!is.list(beta0)) stop("beta0 must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
    }
    recovsCol <- get_all_vars(recharge$theta,data)#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
    if(!all(names(recovsCol) %in% names(data))){
      recovsCol <- recovsCol[,names(recovsCol) %in% names(data),drop=FALSE]
    }
    g0covsCol <- get_all_vars(recharge$g0,data)#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
    if(!all(names(g0covsCol) %in% names(data))){
      g0covsCol <- g0covsCol[,names(g0covsCol) %in% names(data),drop=FALSE]
    }
  } else {
    if(mixtures==1) beta0 <- list(beta=beta0)
    nbRecovs <- 0
    nbG0covs <- 0
    g0covs <- g0covsCol <- NULL
    recovs <- recovsCol <- NULL
    newdata <- NULL
  }
  
  # build design matrix for initial distribution
  if(is.null(formulaDelta)){
    formDelta <- ~1
  } else formDelta <- formulaDelta
  covsDelta <- stats::model.matrix(formDelta,data[aInd,,drop=FALSE])
  if(nrow(covsDelta)!=nrow(data[aInd,,drop=FALSE])) stop("covariates cannot contain missing values")
  nbCovsDelta <- ncol(covsDelta)-1 # substract intercept column
  
  # build design matrix for mixture probabilties
  if(is.null(formulaPi)){
    formPi <- ~1
  } else formPi <- formulaPi
  covsPi <- stats::model.matrix(formPi,data[aInd,,drop=FALSE])
  if(nrow(covsPi)!=nrow(data[aInd,,drop=FALSE])) stop("covariates cannot contain missing values")
  nbCovsPi <- ncol(covsPi)-1 # substract intercept column
  
  # check that stationary==FALSE if there are covariates
  if(nbCovs>0 & stationary==TRUE)
    if(nbAnimals!=sum(mapply(function(x) nrow(unique(covs[which(data$ID==x),,drop=FALSE])),unique(data$ID)))) stop("stationary can't be set to TRUE if there are time-varying covariates in formula.")
  if(nbCovsDelta>0 & stationary==TRUE)
    stop("stationary can't be set to TRUE if there are covariates in formulaDelta.")

  if(length(covsCol)>0 | !is.null(recharge)) {
    if(!is.null(recharge)){
      covsCol <- cbind(covsCol,get_all_vars(recharge$g0,data),get_all_vars(recharge$theta,data))#rownames(attr(stats::terms(formula),"factors"))#attr(stats::terms(formula),"term.labels")#seq(1,ncol(data))[-match(c("ID","x","y",distnames),names(data),nomatch=0)]
    }  
    if(!all(names(covsCol) %in% names(data))){
        covsCol <- covsCol[,names(covsCol) %in% names(data),drop=FALSE]
    }
    rawCovs <- covsCol[,unique(names(covsCol)),drop=FALSE]
  } else {
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
  
  for(i in which(unlist(lapply(dist,function(x) x %in% "bern"))))
    if(any(data[[distnames[[i]]]]!=0 & data[[distnames[[i]]]]!=1,na.rm=TRUE))
      stop(distnames[[i]]," data should be binary (0 or 1)")
  
  # determine whether zero-inflation or one-inflation should be included
  inflation <- get_inflation(data,dist,DM,Par0,nbStates)
  Par0 <- inflation$Par0
  zeroInflation <- inflation$zeroInflation
  oneInflation <- inflation$oneInflation

  mHind <- (requireNamespace("moveHMM", quietly = TRUE) && is.null(DM) & is.null(userBounds) & is.null(workBounds) & ("step" %in% distnames) & is.null(fixPar) & !length(attr(stats::terms.formula(newformula),"term.labels")) & !length(attr(stats::terms.formula(formDelta),"term.labels")) & stationary & optMethod=="nlm" & is.null(prior) & is.null(betaCons) & is.null(betaRef) & is.null(deltaCons) & is.null(mvnCoords) & mixtures==1) # indicator for moveHMMwrap below
  
  inputs <- checkInputs(nbStates,dist,Par0,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,stateNames,checkInflation = TRUE)
  p <- inputs$p
  
  DMinputs<-getDM(data,inputs$DM,inputs$dist,nbStates,p$parNames,p$bounds,Par0,zeroInflation,oneInflation,inputs$circularAngleMean)
  fullDM <- DMinputs$fullDM
  DMind <- DMinputs$DMind

  # check elements of nlmPar
  lsPars <- c("print.level","gradtol","stepmax","steptol","iterlim",'hessian')
  if(length(which(!(names(nlmPar) %in% lsPars)))>0)
    stop("Check the names of the elements of 'nlmPar'; they should be in
         ('print.level','gradtol','stepmax','steptol','iterlim','hessian')")


  ####################################
  ## Prepare initial values for nlm ##
  ####################################
  fixParIndex <- get_fixParIndex(Par0,beta0,delta0,fixPar,distnames,inputs,p,nbStates,DMinputs,recharge,nbG0covs,nbRecovs,workBounds,mixtures,newformula,formulaStates,nbCovs,betaCons,betaRef,deltaCons,stationary,nbCovsDelta,formulaDelta,formulaPi,nbCovsPi,data,newdata)
  
  parCount<- lapply(fullDM,ncol)
  for(i in distnames[!unlist(lapply(inputs$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount))[-length(fullDM)])
  names(parindex) <- distnames
  
  if(is.null(workBounds)) {
    workBounds <- vector('list',length(distnames))
    names(workBounds) <- distnames
  }
  workBounds <- getWorkBounds(workBounds,distnames,fixParIndex$Par0,parindex,parCount,inputs$DM,fixParIndex$beta0,fixParIndex$delta0)
  
  wpar <- n2w(fixParIndex$Par0,p$bounds,fixParIndex$beta0,fixParIndex$delta0,nbStates,inputs$estAngleMean,inputs$DM,p$Bndind,inputs$dist)
  if(any(!is.finite(wpar))) stop("Scaling error. Check initial parameter values and bounds.")

  parmInd <- length(wpar)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures-(nbCovsPi+1)*(mixtures-1)-ncol(covsDelta)*(nbStates-1)*(!stationary)*mixtures-ifelse(is.null(recharge),0,(nbRecovs+1+nbG0covs+1))

  retrySD <- get_retrySD(retryFits,retrySD,wpar,parmInd,distnames,parCount,nbStates,fixParIndex$beta0,mixtures,recharge,fixParIndex$delta0)
  
  optInd <- sort(c(fixParIndex$wparIndex,parmInd+which(duplicated(c(betaCons))),parmInd+length(fixParIndex$beta0$beta)+length(fixParIndex$fixPar[["pi"]])+which(duplicated(c(deltaCons)))))
  
  if(!is.null(prior)){
    if(!is.function(prior)) stop("prior must be a function")
    environment(prior) <- environment()
    pr <- tryCatch(prior(wpar),error=function(e) e)
    if(inherits(pr,"error")){
      stop("Invalid prior: ",pr)
    } else {
      if(length(pr)>1) stop("prior function must return a scalar")
      if(!is.finite(pr)) stop("prior is not finite")
    }
  }
  
  ##################
  ## Optimization ##
  ##################
  printMessage(nbStates,dist,p,DM,formula,formDelta,formPi,mixtures,stationary=stationary,hierarchical=inherits(data,"hierarchical"))
  
  ncmean <- get_ncmean(distnames,fullDM,inputs$circularAngleMean,nbStates)
  nc <- ncmean$nc
  meanind <- ncmean$meanind
  
  optPar <- wpar
  optInd <- sort(c(fixParIndex$wparIndex,parmInd+which(duplicated(c(betaCons))),parmInd+length(fixParIndex$beta0$beta)+length(fixParIndex$fixPar[["pi"]])+which(duplicated(c(deltaCons)))))
  if(length(optInd)){
    optPar <- wpar[-optInd]
  }
  
  # aInd = list of indices of first observation for each animal
  aInd <- NULL
  for(i in 1:nbAnimals){
    idInd <- which(data$ID==unique(data$ID)[i])
    aInd <- c(aInd,idInd[1])
  }
  
  if(fit) {
    
    hessian <- TRUE
    if(!is.null(control$hessian)){
      hessian <- control$hessian
      control$hessian <- NULL
    }
    
    fitCount<-0
    
    while(fitCount<=retryFits){
      
      # just use moveHMM if simpler models are specified
      if(all(distnames %in% c("step","angle")) & all(unlist(dist) %in% moveHMMdists) & mHind){
        fullPar<-w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary,fullDM,DMind,nrow(data),dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
        Par<-lapply(fullPar[distnames],function(x) x[,1])
        for(i in distnames){
          if(dist[[i]] %in% angledists & !inputs$estAngleMean[[i]])
            Par[[i]]<-Par[[i]][-(1:nbStates)]
        }
        startTime <- proc.time()
        withCallingHandlers(curmod<-tryCatch(moveHMMwrap(data,nbStates,dist,Par,fullPar$beta,fullPar$delta[1,],inputs$estAngleMean,newformula,stationary,nlmPar,fit,nbAnimals,knownStates)$mod,error=function(e) e),warning=muffleOPTwarning)
        endTime <- proc.time()-startTime
        curmod$wpar <- curmod$estimate
      } else {

        # check additional parameters for nlm
        print.level <- ifelse(is.null(nlmPar$print.level),0,nlmPar$print.level)
        gradtol <- ifelse(is.null(nlmPar$gradtol),1e-6,nlmPar$gradtol)
        typsize <- rep(1, length(wpar))
        defStepmax <- max(1000 * sqrt(sum((wpar/typsize)^2)),1000)
        stepmax <- ifelse(is.null(nlmPar$stepmax),defStepmax,nlmPar$stepmax)
        steptol <- ifelse(is.null(nlmPar$steptol),1e-6,nlmPar$steptol)
        iterlim <- ifelse(is.null(nlmPar$iterlim),1000,nlmPar$iterlim)
        
        optPar <- wpar
        optInd <- sort(c(fixParIndex$wparIndex,parmInd+which(duplicated(c(betaCons))),parmInd+length(fixParIndex$beta0$beta)+length(fixParIndex$fixPar[["pi"]])+which(duplicated(c(deltaCons)))))
        if(length(optInd)){
          optPar <- wpar[-optInd]
        }
        
        startTime <- proc.time()
  
        # call to optimizer nlm
        if(!length(optPar)){
          curmod <- list()
          
          curmod$minimum <- nLogLike(optPar,nbStates,newformula,p$bounds,p$parSize,data,inputs$dist,covs,
                                    inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,zeroInflation,oneInflation,
                                    stationary,fullDM,DMind,p$Bndind,knownStates,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,
                                    nc,meanind,covsDelta,workBounds,prior,betaCons,fixParIndex$betaRef,deltaCons,optInd,recovs,g0covs,mixtures,covsPi,hierRecharge,aInd)
          curmod$estimate <- numeric()
          
        } else if(optMethod=="nlm"){
          withCallingHandlers(curmod <- tryCatch(nlm(nLogLike,optPar,nbStates,newformula,p$bounds,p$parSize,data,inputs$dist,covs,
                                               inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,zeroInflation,oneInflation,
                                               stationary,fullDM,DMind,p$Bndind,knownStates,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,
                                               nc,meanind,covsDelta,workBounds,prior,betaCons,fixParIndex$betaRef,deltaCons,optInd,recovs,g0covs,mixtures,covsPi,hierRecharge,aInd,
                                               print.level=print.level,gradtol=gradtol,
                                               stepmax=stepmax,steptol=steptol,
                                               iterlim=iterlim,hessian=ifelse(is.null(nlmPar$hessian),TRUE,nlmPar$hessian)),error=function(e) e),warning=muffleOPTwarning)
        } else {
          withCallingHandlers(curmod <- tryCatch(optim(optPar,nLogLike,gr=NULL,nbStates,newformula,p$bounds,p$parSize,data,inputs$dist,covs,
                                                     inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,zeroInflation,oneInflation,
                                                     stationary,fullDM,DMind,p$Bndind,knownStates,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,
                                                     nc,meanind,covsDelta,workBounds,prior,betaCons,fixParIndex$betaRef,deltaCons,optInd,recovs,g0covs,mixtures,covsPi,hierRecharge,aInd,
                                                     method=optMethod,control=control,hessian=hessian),error=function(e) e),warning=muffleOPTwarning)
        }
        endTime <- proc.time()-startTime
      }
      
      if(fitCount==0){
        if(inherits(curmod,"error")) stop(curmod)
        else {
          names(curmod)[which(names(curmod)=="par")] <- "estimate"
          names(curmod)[which(names(curmod)=="value")] <- "minimum"
          curmod$elapsedTime <- endTime[3]
          mod <- curmod
          if(retryFits>=1) cat("Attempting to improve fit using random perturbation. Press 'esc' to force exit from 'fitHMM'\n")
        }
      }
      
      wpar <- expandPar(mod$estimate,optInd,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,betaCons,deltaCons,nbStates,nbCovsDelta,stationary,nbCovs,nbRecovs+nbG0covs,mixtures,nbCovsPi)
      
      if((fitCount+1)<=retryFits){
        cat("\r    Attempt ",fitCount+1," of ",retryFits," -- current log-likelihood value: ",-mod$minimum,"         ",sep="")
        if(!inherits(curmod,"error")){
          names(curmod)[which(names(curmod)=="par")] <- "estimate"
          names(curmod)[which(names(curmod)=="value")] <- "minimum"
          curmod$elapsedTime <- endTime[3]
          if(curmod$minimum < mod$minimum) mod <- curmod
          wpar <- expandPar(mod$estimate,optInd,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,betaCons,deltaCons,nbStates,nbCovsDelta,stationary,nbCovs,nbRecovs+nbG0covs,mixtures,nbCovsPi)
        }
        wpar[1:parmInd] <- wpar[1:parmInd]+rnorm(parmInd,0,retrySD[1:parmInd])
        if(nbStates>1)
          wpar[-(1:parmInd)] <- wpar[-(1:parmInd)]+rnorm(length(wpar)-parmInd,0,retrySD[-(1:parmInd)])
        if(length(fixParIndex$wparIndex)) wpar[fixParIndex$wparIndex] <- unlist(fixParIndex$fixPar)[fixParIndex$wparIndex]
        if(!is.null(betaCons) & nbStates>1){
          wpar[parmInd+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures)] <- wpar[parmInd+1:((nbCovs+1)*nbStates*(nbStates-1)*mixtures)][betaCons]
        }
      }
      fitCount<-fitCount+1
    }
    
    # convert the parameters back to their natural scale
    mod$wpar <- mod$estimate
    wpar <- mod$estimate <- expandPar(mod$wpar,optInd,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,betaCons,deltaCons,nbStates,nbCovsDelta,stationary,nbCovs,nbRecovs+nbG0covs,mixtures,nbCovsPi)
    mle <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary,fullDM,DMind,nrow(data),inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
    
    if(!is.null(mod$hessian)){
      Sigma <- tryCatch(MASS::ginv(mod$hessian),error=function(e) e)
      mod$Sigma <- Sigma
      if(!inherits(Sigma,"error")){
        if(length(optInd)){
          mod$Sigma <- matrix(0,length(mod$estimate),length(mod$estimate))
          mod$Sigma[(1:length(mod$estimate))[-optInd],(1:length(mod$estimate))[-optInd]] <- Sigma
        }
      } else {
        warning("ginv of the hessian failed -- ",Sigma)
      }
    }
  }
  else {
    mod <- list()
    mod$minimum <- nLogLike(optPar,nbStates,newformula,p$bounds,p$parSize,data,inputs$dist,covs,
                               inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,zeroInflation,oneInflation,
                               stationary,fullDM,DMind,p$Bndind,knownStates,unlist(fixParIndex$fixPar),fixParIndex$wparIndex,
                               nc,meanind,covsDelta,workBounds,prior,betaCons,fixParIndex$betaRef,deltaCons,optInd,recovs,g0covs,mixtures,covsPi,hierRecharge,aInd)
    mod$estimate <- wpar
    mod$wpar <- optPar
    mle <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,inputs$estAngleMean,inputs$circularAngleMean,inputs$consensus,stationary,fullDM,DMind,nrow(data),inputs$dist,p$Bndind,nc,meanind,covsDelta,workBounds,covsPi)
  }
  

  ####################
  ## Prepare output ##
  ####################
  if(is.null(stateNames)){
    for(i in 1:nbStates)
      stateNames[i] <- paste("state",i)
  }
  
  # name columns and rows of MLEs
  for(i in distnames){
    if(dist[[i]] %in% angledists)
      if(!inputs$estAngleMean[[i]]){
        p$parNames[[i]] <- c("mean",p$parNames[[i]])
      }
    if(DMind[[i]]){
      mle[[i]]<-matrix(mle[[i]][1:(length(p$parNames[[i]])*nbStates),1],nrow=length(p$parNames[[i]]),ncol=nbStates,byrow=TRUE)
      rownames(mle[[i]]) <- p$parNames[[i]]
      colnames(mle[[i]]) <- stateNames
    } else {
      mle[[i]]<-matrix(wpar[parindex[[i]]+1:parCount[[i]]],1)
      rownames(mle[[i]])<-"[1,]"
      if(!isFALSE(inputs$circularAngleMean[[i]])){
        colnames(mle[[i]]) <- unique(gsub("cos","",gsub("sin","",colnames(fullDM[[i]]))))
      } else colnames(mle[[i]])<-colnames(fullDM[[i]])
      #if(is.null(names(mle[[i]]))) warning("No names for the regression coeffs were provided in DM$",i)
    }
  }

  if(!is.null(mle$beta)) {
    rownames(mle$beta) <- rep(colnames(covs),mixtures)#c("(Intercept)",attr(stats::terms(formula),"term.labels"))
    columns <- NULL
    for(i in 1:nbStates)
      for(j in 1:nbStates) {
        if(fixParIndex$betaRef[i]<j)
          columns[(i-1)*nbStates+j-i] <- paste(i,"->",j)
        if(j<fixParIndex$betaRef[i])
            columns[(i-1)*(nbStates-1)+j] <- paste(i,"->",j)
      }
    colnames(mle$beta) <- columns
    if(!is.null(recharge)){
      names(mle$g0)<-colnames(g0covs)
      names(mle$theta)<-colnames(recovs)
    }
    if(mixtures>1){
      rownames(mle$beta) <- paste0(rownames(mle$beta),"_mix",rep(1:mixtures,each=length(colnames(covs))))
      colnames(mle[["pi"]]) <- paste0("mix",1:mixtures)
      rownames(mle[["pi"]]) <- paste0("ID:",data$ID[aInd])
    } else mle[["pi"]] <- NULL
  }

  # compute stationary distribution
  if(stationary) {
    mle$delta <- matrix(0,nbAnimals*mixtures,nbStates)
    
    for(mix in 1:mixtures){
      
      for(i in 1:nbAnimals){
        
        gamma <- trMatrix_rcpp(nbStates,mle$beta[(nbCovs+1)*(mix-1)+1:(nbCovs+1),,drop=FALSE],covs,fixParIndex$betaRef)[,,aInd[i]]
    
        # error if singular system
        tryCatch(
          mle$delta[nbAnimals*(mix-1)+i,] <- solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)),
          error = function(e) {
            stop(paste("A problem occurred in the calculation of the stationary",
                       "distribution. You may want to try different initial values",
                       "and/or the option stationary=FALSE."))
          }
        )
      }
    }
  }

  if(nbStates==1)
    mle$delta <- matrix(1,nrow(covsDelta),1)
  colnames(mle$delta) <- stateNames
  rownames(mle$delta) <- paste0("ID:",rep(data$ID[aInd],mixtures))
  if(mixtures>1) rownames(mle$delta) <- paste0("ID:",rep(data$ID[aInd],mixtures),"_mix",rep(1:mixtures,each=nbAnimals))
  
  # compute t.p.m. if no covariates
  if(nbCovs==0 & nbStates>1) {
    mle$gamma <- matrix(0,nbStates*mixtures,nbStates)
    for(mix in 1:mixtures){
      trMat <- trMatrix_rcpp(nbStates,mle$beta[(nbCovs+1)*(mix-1)+1:(nbCovs+1),,drop=FALSE],covs,fixParIndex$betaRef)
      mle$gamma[nbStates*(mix-1)+1:nbStates,] <- trMat[,,1]
    }
    colnames(mle$gamma)<-stateNames
    rownames(mle$gamma)<-rep(stateNames,mixtures)
    if(mixtures>1) rownames(mle$gamma) <- paste0(rep(stateNames,mixtures),"_mix",rep(1:mixtures,each=nbStates))
  }
  
  if(is.null(betaCons) & nbStates>1) betaCons <- matrix(1:length(mle$beta),nrow(mle$beta),ncol(mle$beta))

  # conditions of the fit
  conditions <- list(dist=dist,zeroInflation=zeroInflation,oneInflation=oneInflation,
                     estAngleMean=inputs$estAngleMean,circularAngleMean=inputs$circularAngleMean,stationary=stationary,formula=formula,userBounds=userBounds,workBounds=workBounds,bounds=p$bounds,Bndind=p$Bndind,DM=DM,fullDM=fullDM,DMind=DMind,fixPar=fixParIndex$ofixPar,wparIndex=fixParIndex$wparIndex,formulaDelta=formulaDelta,betaCons=betaCons,betaRef=fixParIndex$betaRef,deltaCons=deltaCons,optInd=optInd,recharge=recharge,mvnCoords=mvnCoords,mixtures=mixtures,formulaPi=formulaPi,fit=fit)

  mh <- list(data=data,mle=mle,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames,knownStates=knownStates,covsDelta=covsDelta,prior=prior,modelName=modelName,reCovs=recovsCol,g0covs=g0covsCol,covsPi=covsPi)
  
  #compute SEs and CIs on natural and working scale
  CIreal<-tryCatch(CIreal(momentuHMM(mh)),error=function(e) e)
  if(inherits(CIreal,"error") & fit==TRUE) warning("Failed to compute SEs and confidence intervals on the natural scale -- ",CIreal)
  CIbeta<-tryCatch(CIbeta(momentuHMM(mh)),error=function(e) e)
  if(inherits(CIbeta,"error") & fit==TRUE) warning("Failed to compute SEs confidence intervals on the working scale -- ",CIbeta)
  
  mh <- list(data=data,mle=mle,CIreal=CIreal,CIbeta=CIbeta,mod=mod,conditions=conditions,rawCovs=rawCovs,stateNames=stateNames,knownStates=knownStates,covsDelta=covsDelta,prior=prior,modelName=modelName,reCovs=recovsCol,g0covs=g0covsCol,covsPi=covsPi)
  
  if(!is.null(mvnCoords)) attr(mh$data,'coords') <- paste0(mvnCoords,c(".x",".y"))
  
  if(fit) message(ifelse(retryFits>=1,"\n",""),"DONE")
  
  return(momentuHMM(mh))
}

#' @rdname fitHMM
#' @method fitHMM momentuHierHMMData
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). See details.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See details.
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See details.
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (only hierarchical-level effects, with no covariate effects).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy. See details.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' 
#' @details
#' \itemize{
#' \item \code{fitHMM.momentuHierHMMData} is very similar to \code{\link{fitHMM.momentuHMMData}} except that instead of simply specifying the number of states (\code{nbStates}), distributions (\code{dist}), and a single t.p.m. formula (\code{formula}), the \code{hierStates} argument specifies the hierarchical nature of the states,
#' the \code{hierDist} argument specifies the hierarchical nature of the data streams, and the \code{hierFormula} argument specifies a t.p.m. formula for each level of the hierarchy.  All are specified as 
#' \code{\link[data.tree]{Node}} objects from the \code{\link[data.tree]{data.tree}} package.
#' }
#' 
#' @seealso \code{\link{simHierData}}
#' 
#' @export
fitHMM.momentuHierHMMData <- function(data,hierStates,hierDist,
                       Par0,hierBeta=NULL,hierDelta=NULL,
                       estAngleMean=NULL,circularAngleMean=NULL,
                       hierFormula=NULL,hierFormulaDelta=NULL,mixtures=1,formulaPi=NULL,
                       nlmPar=list(),fit=TRUE,
                       DM=NULL,userBounds=NULL,workBounds=NULL,betaCons=NULL,deltaCons=NULL,
                       mvnCoords=NULL,knownStates=NULL,fixPar=NULL,retryFits=0,retrySD=NULL,optMethod="nlm",control=list(),prior=NULL,modelName=NULL, ...)
{
  
  if(!inherits(data,"momentuHierHMMData")) stop("data must be a momentuHierHMMData object (as returned by prepData or simHierData)")
  
  installDataTree()
  
  inputHierHMM <- formatHierHMM(data,hierStates,hierDist,hierBeta,hierDelta,hierFormula,hierFormulaDelta,mixtures,workBounds,betaCons,deltaCons,fixPar)
  
  chkDots(...)
  
  hfit <- fitHMM.momentuHMMData(momentuHMMData(data),inputHierHMM$nbStates,inputHierHMM$dist,Par0,inputHierHMM$beta,inputHierHMM$delta,
                estAngleMean,circularAngleMean,
                formula=inputHierHMM$formula,inputHierHMM$formulaDelta,stationary=FALSE,mixtures,formulaPi,
                nlmPar,fit,
                DM,userBounds,workBounds=inputHierHMM$workBounds,inputHierHMM$betaCons,inputHierHMM$betaRef,deltaCons=inputHierHMM$deltaCons,
                mvnCoords,inputHierHMM$stateNames,knownStates,inputHierHMM$fixPar,retryFits,retrySD,optMethod,control,prior,modelName)
  
  # replace initial values with estimates in hierBeta and hierDelta (if provided)
  if(fit){
    par <- getPar(hfit)
    if(is.list(par$beta)){
      beta <- par$beta$beta
      Pi <- par$beta[["pi"]]
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
    hfit$CIreal <- CIreal.hierarchical(hfit)
  }
  hfit
}
