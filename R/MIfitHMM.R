
#' Fit HMMs to multiple imputation data
#' 
#' Fit a (multivariate) hidden Markov model to multiple imputation data. Multiple imputation is a method for accommodating 
#' missing data, temporal-irregularity, or location measurement error in hidden Markov models, where pooled parameter estimates reflect uncertainty
#' attributable to observation error.
#' 
#' \code{miData} can either be a \code{\link{crwData}} or \code{\link{crwHierData}} object (as returned by \code{\link{crawlWrap}}), a \code{\link{crwSim}} or \code{\link{crwHierSim}} object (as returned by \code{MIfitHMM} when \code{fit=FALSE}), 
#' or a list of \code{\link{momentuHMMData}} or \code{\link{momentuHierHMMData}} objects (e.g., each element of the list as returned by \code{\link{prepData}}). 
#' 
#' If \code{miData} is a \code{crwData} (or \code{crwHierData}) object, \code{MIfitHMM} uses a combination of 
#' \code{\link[crawl]{crwSimulator}}, \code{\link[crawl]{crwPostIS}}, \code{\link{prepData}}, and \code{\link{fitHMM}} to draw \code{nSims} realizations of the position process
#' and fit the specified HMM to each imputation of the data. The vast majority of \code{MIfitHMM} arguments are identical to the corresponding arguments from these functions.
#' 
#' If \code{miData} is a \code{\link{crwData}} or \code{\link{crwHierData}} object, \code{nSims} determines both the number of realizations of the position process to draw 
#' (using \code{\link[crawl]{crwSimulator}} and \code{\link{crwPostIS}}) as well as the number of HMM fits.
#' 
#' If \code{miData} is a \code{\link{crwSim}} (or \code{\link{crwHierSim}}) object or a list of \code{\link{momentuHMMData}} (or \code{\link{momentuHierHMMData}}) object(s), the specified HMM will simply be fitted to each of the \code{momentuHMMData} (or \code{momentuHierHMMData}) objects
#' and all arguments related to \code{\link[crawl]{crwSimulator}}, \code{\link[crawl]{crwPostIS}}, or \code{\link{prepData}} are ignored.
#' 
#' @param miData A \code{\link{crwData}} object, a \code{\link{crwHierData}} object, a \code{\link{crwSim}} object, a \code{\link{crwHierSim}} object, a list of \code{\link{momentuHMMData}} objects, or a list of \code{\link{momentuHierHMMData}} objects.
#' @param ... further arguments passed to or from other methods
#' @export
MIfitHMM <- function(miData, ...) {
  UseMethod("MIfitHMM")
}

#' @rdname MIfitHMM
#' @method MIfitHMM default
#' @param nSims Number of imputations in which to fit the HMM using \code{\link{fitHMM}}. If \code{miData} is a list of \code{momentuHMMData} 
#' objects, \code{nSims} cannot exceed the length of \code{miData}.
#' @param ncores Number of cores to use for parallel processing. Default: 1 (no parallel processing).
#' @param poolEstimates Logical indicating whether or not to calculate pooled parameter estimates across the \code{nSims} imputations using \code{\link{MIpool}}. Default: \code{TRUE}.
#' @param alpha Significance level for calculating confidence intervals of pooled estimates when \code{poolEstimates=TRUE} (see \code{\link{MIpool}}). Default: 0.95.
#' @param na.rm Logical indicating whether or not to exclude model fits with \code{NA} parameter estimates or standard errors from pooling when \code{poolEstimates=TRUE} (see \code{\link{MIpool}}). Default: FALSE.
#' @param nbStates Number of states of the HMM. See \code{\link{fitHMM}}.
#' @param dist A named list indicating the probability distributions of the data streams. See \code{\link{fitHMM}}.
#' @param Par0 A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. See \code{\link{fitHMM}}.  \code{Par0} may also be a list of length \code{nSims}, where each element is a named list containing vectors
#' of initial state-dependent probability distribution parameters for each imputation. Note that if \code{useInitial=TRUE} then \code{Par0} is ignored after the first imputation.
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities. See \code{\link{fitHMM}}. \code{beta0} may also be a list of length \code{nSims}, where each element 
#' is an initial matrix of regression coefficients for the transition probabilities for each imputation.
#' @param delta0 Initial values for the initial distribution of the HMM. See \code{\link{fitHMM}}. \code{delta0} may also be a list of length \code{nSims}, where each element 
#' is the initial values for the initial distribution of the HMM for each imputation.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). See \code{\link{fitHMM}}.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear or circular-circular
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. See \code{\link{fitHMM}}.
#' @param formula Regression formula for the transition probability covariates. See \code{\link{fitHMM}}.
#' @param formulaDelta Regression formula for the initial distribution. See \code{\link{fitHMM}}.
#' @param stationary \code{FALSE} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. See \code{\link{fitHMM}}.
#' @param mixtures Number of mixtures for the state transition probabilities  (i.e. discrete random effects *sensu* DeRuiter et al. 2017). Default: \code{mixtures=1}.
#' @param formulaPi Regression formula for the mixture distribution probabilities. See \code{\link{fitHMM}}.
#' @param nlmPar List of parameters to pass to the optimization function \code{\link[stats]{nlm}} (which should be either
#' \code{print.level}, \code{gradtol}, \code{stepmax}, \code{steptol}, \code{iterlim}, or \code{hessian} -- see \code{nlm}'s documentation
#' for more detail). For \code{print.level}, the default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed. Ignored unless \code{optMethod="nlm"}.
#' @param fit \code{TRUE} if the HMM should be fitted to the data, \code{FALSE} otherwise. See \code{\link{fitHMM}}. If \code{fit=FALSE} and \code{miData} is a \code{\link{crwData}}
#' object, then \code{MIfitHMM} returns a list containing a \code{\link{momentuHMMData}} object (if \code{nSims=1}) or, if \code{nSims>1}, a \code{\link{crwSim}} object.
#' @param useInitial Logical indicating whether or not to use parameter estimates for the first model fit as initial values for all subsequent model fits.  
#' If \code{ncores>1} then the first model is fit on a single core and then used as the initial values for all subsequent model fits on each core 
#' (in this case, the progress of the initial model fit can be followed using the \code{print.level} option in \code{nlmPar}). Default: FALSE. Ignored if \code{nSims<2}.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. See \code{\link{fitHMM}}.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. See \code{\link{fitHMM}}.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. See \code{\link{fitHMM}}.
#' @param betaCons Matrix of the same dimension as \code{beta0} composed of integers identifying any equality constraints among the t.p.m. parameters. See \code{\link{fitHMM}}.
#' @param betaRef Numeric vector of length \code{nbStates} indicating the reference elements for the t.p.m. multinomial logit link. See \code{\link{fitHMM}}. 
#' @param deltaCons Matrix of the same dimension as \code{delta0} composed of integers identifying any equality constraints among the initial distribution working scale parameters. Ignored unless a formula is provided in \code{formulaDelta}. See \code{\link{fitHMM}}.
#' @param mvnCoords Character string indicating the name of location data that are to be modeled using a multivariate normal distribution. For example, if \code{mu="mvnorm2"} was included in \code{dist} and (mu.x, mu.y) are location data, 
#' then \code{mvnCoords="mu"} needs to be specified in order for these data to be properly treated as locations in functions such as \code{\link{plot.momentuHMM}}, \code{\link{plot.miSum}}, \code{\link{plot.miHMM}}, \code{\link{plotSpatialCov}}, and \code{\link{MIpool}}.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any). See \code{\link{fitHMM}}. If \code{miData} is a list of \code{\link{momentuHMMData}} objects, then \code{knownStates} can alternatively
#' be a list of vectors containing the known values for the state process for each element of \code{miData}.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. See \code{\link{fitHMM}}. 
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization.  See \code{\link{fitHMM}}. 
#' @param retrySD An optional list of scalars or vectors indicating the standard deviation to use for normal perturbations of each working scale parameter when \code{retryFits>0}. See \code{\link{fitHMM}}.
#' @param optMethod The optimization method to be used.  Can be ``nlm'' (the default; see \code{\link[stats]{nlm}}), ``Nelder-Mead'' (see \code{\link[stats]{optim}}), or ``SANN'' (see \code{\link[stats]{optim}}).
#' @param control A list of control parameters to be passed to \code{\link[stats]{optim}} (ignored unless \code{optMethod="Nelder-Mead"} or \code{optMethod="SANN"}).
#' @param prior A function that returns the log-density of the working scale parameter prior distribution(s).  See \code{\link{fitHMM}}.
#' @param modelName An optional character string providing a name for the fitted model. If provided, \code{modelName} will be returned in \code{\link{print.momentuHMM}}, \code{\link{AIC.momentuHMM}}, \code{\link{AICweights}}, and other functions. 
#' @param covNames Names of any covariates in \code{miData$crwPredict} (if \code{miData} is a \code{\link{crwData}} object; otherwise 
#' \code{covNames} is ignored). See \code{\link{prepData}}. If \code{miData} is a \code{\link{crwData}} object, any covariate in \code{miData$crwPredict} that is used in \code{formula}, \code{formulaDelta}, \code{formulaPi}, or \code{DM} must be included in \code{covNames}.
#' @param spatialCovs List of raster layer(s) for any spatial covariates. See \code{\link{prepData}}. 
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on realizations of the position process. 
#' See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param centroids List where each element is a data frame containing the x-coordinates ('x'), y-coordinates ('y'), and times (with user-specified name, e.g., `time') for centroids (i.e., dynamic activity centers where the coordinates can change over time)
#' from which distance and angle covariates will be calculated based on the location data. See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{miData$crwPredict} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} or \code{\link{crwHierData}} object.
#' @param altCoordNames Character string indicating an alternative name for the returned location data. See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} or \code{\link{crwHierData}} object.
#' @param method Method for obtaining weights for movement parameter samples. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param parIS Size of the parameter importance sample. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param dfSim Degrees of freedom for the t approximation to the parameter posterior. See 'df' argument in \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param grid.eps Grid size for \code{method="quadrature"}. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param crit Criterion for deciding "significance" of quadrature points
#' (difference in log-likelihood). See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param scaleSim Scale multiplier for the covariance matrix of the t approximation. See 'scale' argument in \code{\link[crawl]{crwSimulator}}. 
#' Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param quad.ask Logical, for method='quadrature'. Whether or not the sampler should ask if quadrature sampling should take place. It is used to stop the sampling if the number of likelihood evaluations would be extreme. Default: FALSE. Ignored if \code{ncores>1}.
#' @param force.quad A logical indicating whether or not to force the execution 
#' of the quadrature method for large parameter vectors. See \code{\link[crawl]{crwSimulator}}. Default: TRUE. Ignored unless \code{miData} is a \code{\link{crwData}} object and \code{method=``quadrature''}.
#' @param fullPost Logical indicating whether to draw parameter values as well to simulate full posterior. See \code{\link[crawl]{crwPostIS}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param dfPostIS Degrees of freedom for multivariate t distribution approximation to parameter posterior. See 'df' argument in \code{\link[crawl]{crwPostIS}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param scalePostIS Extra scaling factor for t distribution approximation. See 'scale' argument in \code{\link[crawl]{crwPostIS}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param thetaSamp If multiple parameter samples are available in \code{\link[crawl]{crwSimulator}} objects,
#' setting \code{thetaSamp=n} will use the nth sample. Defaults to the last. See \code{\link[crawl]{crwSimulator}} and \code{\link[crawl]{crwPostIS}}. 
#' Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' 
#' @return  If \code{nSims>1}, \code{poolEstimates=TRUE}, and \code{fit=TRUE}, a \code{\link{miHMM}} object, i.e., a list consisting of:
#' \item{miSum}{\code{\link{miSum}} object returned by \code{\link{MIpool}}.}
#' \item{HMMfits}{List of length \code{nSims} comprised of \code{\link{momentuHMM}} objects.}
#' If \code{poolEstimates=FALSE} and \code{fit=TRUE}, a list of length \code{nSims} consisting of \code{\link{momentuHMM}} objects is returned. 
#' 
#' However, if \code{fit=FALSE} and \code{miData} is a \code{\link{crwData}} 
#' object, then \code{MIfitHMM} returns a \code{\link{crwSim}} object, i.e., a list containing \code{miData} (a list of 
#' \code{\link{momentuHMMData}} objects) and \code{crwSimulator} (a list of \code{\link[crawl]{crwSimulator}} objects),and most other arguments related to \code{\link{fitHMM}} are ignored.
#' 
#' @seealso \code{\link{crawlWrap}}, \code{\link[crawl]{crwPostIS}}, \code{\link[crawl]{crwSimulator}}, \code{\link{fitHMM}}, \code{\link{getParDM}}, \code{\link{MIpool}}, \code{\link{prepData}} 
#' 
#' @examples
#' \dontshow{
#' set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
#' }
#' # Don't run because it takes too long on a single core
#' \dontrun{
#' # extract simulated obsData from example data
#' obsData <- miExample$obsData
#' 
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#'
#' # create crwData object by fitting crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut <- crawlWrap(obsData=obsData,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
#'          err.model=err.model)
#' 
#' # HMM specifications
#' nbStates <- 2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar0 <- c(mu0,sigma0)
#' anglePar0 <- c(-pi/2,pi/2,kappa0)
#' formula <- ~cov1+cos(cov2)
#' nbCovs <- 2
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' 
#' # first fit HMM to best predicted position process
#' bestData<-prepData(crwOut,covNames=c("cov1","cov2"))
#' bestFit<-fitHMM(bestData,
#'                 nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                 Par0=list(step=stepPar0,angle=anglePar0),beta0=beta0,
#'                 formula=formula,estAngleMean=list(angle=TRUE))
#'             
#' print(bestFit)
#' 
#' # extract estimates from 'bestFit'
#' bPar0 <- getPar(bestFit)
#' 
#' # Fit nSims=5 imputations of the position process
#' miFits<-MIfitHMM(miData=crwOut,nSims=5,
#'                   nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                   Par0=bPar0$Par,beta0=bPar0$beta,delta0=bPar0$delta,
#'                   formula=formula,estAngleMean=list(angle=TRUE),
#'                   covNames=c("cov1","cov2"))
#'
#' # print pooled estimates
#' print(miFits)
#'}
#' 
#' @references
#' 
#' Hooten M.B., Johnson D.S., McClintock B.T., Morales J.M. 2017. Animal Movement: Statistical Models for Telemetry Data. CRC Press, Boca Raton.
#' 
#' McClintock B.T. 2017. Incorporating telemetry error into hidden Markov movement models using multiple imputation. Journal of Agricultural, Biological,
#' and Environmental Statistics.
#' 
#' @export
#' @importFrom crawl crwPostIS crwSimulator
# @importFrom raster getZ
#' @importFrom stats terms.formula
MIfitHMM.default<-function(miData,nSims, ncores = 1, poolEstimates = TRUE, alpha = 0.95, na.rm = FALSE,
                   nbStates, dist, 
                   Par0, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, circularAngleMean = NULL,
                   formula = ~1, formulaDelta = NULL, stationary = FALSE, mixtures = 1, formulaPi = NULL,
                   nlmPar = NULL, fit = TRUE, useInitial = FALSE,
                   DM = NULL, userBounds = NULL, workBounds = NULL, betaCons = NULL, betaRef = NULL, deltaCons = NULL,
                   mvnCoords = NULL, stateNames = NULL, knownStates = NULL, fixPar = NULL, retryFits = 0, retrySD = NULL, optMethod = "nlm", control = list(), prior = NULL, modelName = NULL,
                   covNames = NULL, spatialCovs = NULL, centers = NULL, centroids = NULL, angleCovs = NULL, altCoordNames = NULL, 
                   method = "IS", parIS = 1000, dfSim = Inf, grid.eps = 1, crit = 2.5, scaleSim = 1, quad.ask = FALSE, force.quad = TRUE,
                   fullPost = TRUE, dfPostIS = Inf, scalePostIS = 1,thetaSamp = NULL, ...)
{
  
  # for directing list of momentuHierHMMData objects to MIfitHMM.hierarchical
  hierArgs <- list(...)
  argNames <- names(hierArgs)[which(names(hierArgs) %in% c("hierStates","hierDist","hierBeta","hierDelta","hierFormula","hierFormulaDelta"))]
  
  if(is.list(miData) & !is.crwData(miData) & !is.crwSim(miData)){
    if(missing(nbStates) & missing(dist)){
        if(all(c("hierStates","hierDist") %in% argNames)){
          if(is.null(formulaDelta)) formulaDelta <- ~1
          if(length(attr(stats::terms.formula(formula),"term.labels"))>0 && is.null(hierArgs$hierFormula)) stop("hierFormula should be specified instead of formula")
          if(length(attr(stats::terms.formula(formulaDelta),"term.labels"))>0 && is.null(hierArgs$hierFormulaDelta)) stop("hierFormulaDelta should be specified instead of formulaDelta")
          return(MIfitHMM.hierarchical(miData,nSims, ncores, poolEstimates, alpha,
                                       hierArgs$hierStates, hierArgs$hierDist, 
                                       Par0, hierArgs$hierBeta, hierArgs$hierDelta,
                                       estAngleMean, circularAngleMean,
                                       hierArgs$hierFormula, hierArgs$hierFormulaDelta, mixtures, formulaPi,
                                       nlmPar, fit, useInitial,
                                       DM, userBounds, workBounds, betaCons, deltaCons,
                                       mvnCoords, knownStates, fixPar, retryFits, retrySD, optMethod, control, prior, modelName,
                                       covNames, spatialCovs, centers, centroids, angleCovs, altCoordNames,
                                       method, parIS, dfSim, grid.eps, crit, scaleSim, quad.ask, force.quad,
                                       fullPost, dfPostIS, scalePostIS,thetaSamp))
        }
    }
  }
  if(!missing(nbStates) | !missing(dist)){
    if(any(c("hierStates","hierDist") %in% argNames))
      stop("Either nbStates and dist must be specified (for a regular HMM) or hierStates and hierDist must be specified (for a hierarchical HMM)")
  }
  
  chkDots(...)
  
  j <- mf <- mD <- NULL #gets rid of no visible binding for global variable NOTE in R cmd check
  
  if(poolEstimates){
    if(optMethod=="nlm" & !is.null(nlmPar$hessian)){
      if(!nlmPar$hessian) stop("estimates cannot be pooled unless hessian is calculated")
    } else if(optMethod %in% fitMethods[-1] & !is.null(control$hessian)){
      if(!control$hessian) stop("estimates cannot be pooled unless hessian is calculated")
    }
    if (!requireNamespace("mitools", quietly = TRUE)) {
      stop("Package \"mitools\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  
  if(ncores>1){
    for(pkg in c("doFuture","future")){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
    progBar <- progBar
    pkgs <- c("momentuHMM")
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- c("momentuHMM")
  }
  
  if(is.crwData(miData)){
    
    if(nSims>=1) {
      
      model_fits <- miData$crwFits
      predData <- miData$crwPredict
  
      Time.name<-attr(predData,"Time.name")
      ids = as.character(unique(predData$ID))
      
      if(!is.null(covNames) | !is.null(angleCovs)){
        covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
        if(!length(covNames)) covNames <- NULL
      }
      znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
      if(length(znames))
        if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(miData)),"$crwPredict")
      #coordNames <- attr(predData,"coord")
      
      if(fit | !missing("dist")) {
        if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
        distnames <- tmpdistnames <- names(dist)[which(!(names(dist) %in% c("step","angle")))]
        if(any(is.na(match(distnames,names(predData))))){
          for(i in which(is.na(match(distnames,names(predData))))){
            if(dist[[distnames[i]]] %in% mvndists){
              if(dist[[distnames[i]]] %in% c("mvnorm2","rw_mvnorm2")){
                tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
              } else if(dist[[distnames[i]]] %in% c("mvnorm3","rw_mvnorm3")){
                tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
              }
            }
          }
          if(any(is.na(match(tmpdistnames,names(predData))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(predData)))],collapse=", ")," not found in miData")
          tmpdistnames <- tmpdistnames[which(!(tmpdistnames %in% c("mu.x","mu.y")))]
        }
      } else {
        distnames <- tmpdistnames <- names(predData)[which(!(names(predData) %in% c("ID",Time.name,c("mu.x","mu.y"),covNames,znames)))]
      }
      cat('Drawing ',nSims,' realizations from the position process using crawl... \n',sep="")
      if(ncores>1) message("Running simulator in parallel... ")
      nbAnimals <- length(ids)
      predTimes <- lapply(ids,function(x) predData[[Time.name]][which(predData$ID==x & predData$locType=="p")])
      names(predTimes) <- ids
      withCallingHandlers(crwSim <- foreach(mf = model_fits, i = ids, .export="crwSimulator", .packages=pkgs) %dorng% {
          if(ncores==1) message("Running simulator for individual ",i,"... ",sep="")
          else progBar(which(ids==i), length(ids))
          crSim <- suppressMessages(crawl::crwSimulator(mf,predTime=predTimes[[i]], method = method, parIS = parIS,
                                  df = dfSim, grid.eps = grid.eps, crit = crit, scale = scaleSim, quad.ask = ifelse(ncores>1, FALSE, quad.ask), force.quad = force.quad))
          return(crSim)
      },warning=muffleRNGwarning)
      if(ncores==1) cat("DONE\n")
      names(crwSim) <- ids
      
      
      if(ncores>1) message("Drawing imputations in parallel... ",sep="")
      withCallingHandlers(miData<-
        foreach(j = 1:nSims, .export=c("crwPostIS","prepData"), .errorhandling="pass", .packages=pkgs) %dorng% {
          if(ncores==1) message("\rDrawing imputation ",j,"... ",sep="")
          else progBar(j, nSims)
          locs<-data.frame()
          for(i in 1:length(ids)){
              tmp<-tryCatch({crawl::crwPostIS(crwSim[[i]], fullPost = fullPost, df = dfPostIS, scale = scalePostIS, thetaSamp = thetaSamp)},error=function(e) e)
              if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',ids[i],'; ',tmp,'  Check crwPostIS arguments, crawl::crwMLE model fits, and/or consult crawl documentation.')
              locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
          }
          df<-data.frame(x=locs$mu.x,y=locs$mu.y,predData[,c("ID",Time.name,tmpdistnames,covNames,znames),drop=FALSE])[which(predData$locType=="p"),]
          pD <- tryCatch(prepData(df,covNames=covNames,spatialCovs=spatialCovs,centers=centers,centroids=centroids,angleCovs=angleCovs,altCoordNames=altCoordNames),error=function(e) e)
          return(pD)
        }
      ,warning=muffleRNGwarning)
      if(ncores==1) cat("DONE\n")
      
      for(i in which(unlist(lapply(miData,function(x) inherits(x,"error"))))){
        warning('prepData failed for imputation ',i,"; ",miData[[i]])
      }
      ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHMMData"))))
      if(fit) cat('Fitting',length(ind),'realizations of the position process using fitHMM... \n')
      else {
        if(ncores==1) doParallel::stopImplicitCluster()
        else future::plan(future::sequential)
        return(crwSim(list(miData=miData,crwSimulator=crwSim)))
      }
    } else stop("nSims must be >0")
    
  } else {
    if(!is.list(miData)) stop("miData must either be a crwData object (as returned by crawlWrap) or a list of momentuHMMData objects as returned by simData, prepData, or MIfitHMM (when fit=FALSE)")
    if(is.crwSim(miData)) miData <- miData$miData
    ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHMMData"))))
    if(!length(ind)) stop("miData must either be a crwData object (as returned by crawlWrap) or a list of momentuHMMData objects as returned by simData, prepData, or MIfitHMM (when fit=FALSE)")
    if(missing(nSims)) nSims <- length(miData)
    if(nSims>length(miData)) stop("nSims is greater than the length of miData. nSims must be <=",length(miData))
    if(nSims<1) stop("nSims must be >0")
    cat('Fitting',min(nSims,length(ind)),'imputation(s) using fitHMM... \n')
  }
  
  if(!is.list(knownStates)){
    tmpStates<-knownStates
    knownStates<-vector('list',nSims)
    if(!is.null(tmpStates))
      knownStates[1:nSims]<-list(tmpStates)
  } else if(length(knownStates)<nSims) stop("knownStates must be a list of length >=",nSims)
  
  if(all(names(dist) %in% names(Par0))){
    tmpPar0<-Par0
    Par0<-vector('list',nSims)
    Par0[1:nSims]<-list(tmpPar0)
  } else if(length(Par0)<nSims) stop("Par0 must be a list of length >=",nSims)
  
  newForm <- newFormulas(formula,nbStates,betaRef)
  recharge <- newForm$recharge
  
  if(is.null(recharge) & mixtures==1){
    if(!is.list(beta0)){
      tmpbeta0<-beta0
      beta0<-vector('list',nSims)
      if(!is.null(tmpbeta0))
          beta0[1:nSims]<-list(tmpbeta0)
    } else if(length(beta0)<nSims) stop("beta0 must be a list of length >=",nSims)
  } else {
    if(!is.null(beta0) && !is.list(beta0)){
      if(!is.null(recharge)) stop("beta0 must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
      if(mixtures>1) stop("beta0 must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
    if(!is.list(beta0[[1]])){
      tmpbeta0<-beta0
      beta0<-vector('list',nSims)
      if(!is.null(tmpbeta0))
        beta0[1:nSims]<-list(tmpbeta0)
    } else if(length(beta0)<nSims) stop("beta0 must be a list of length >=",nSims)    
  }
  if(!is.list(delta0)){
    tmpdelta0<-delta0
    delta0<-vector('list',nSims)
    if(!is.null(tmpdelta0))
      delta0[1:nSims]<-list(tmpdelta0)
  } else if(length(delta0)<nSims) stop("delta0 must be a list of length >=",nSims)
  
  #check HMM inputs and print model message
  test<-fitHMM.momentuHMMData(miData[[ind[1]]],nbStates, dist, Par0[[ind[1]]], beta0[[ind[1]]], delta0[[ind[1]]],
           estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
           nlmPar, fit = FALSE, DM,
           userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates[[ind[1]]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName)
  
  # fit HMM(s)
  fits <- list()
  parallelStart <- 1
  useInitial <- ifelse(nSims<2,FALSE,useInitial)
  if(useInitial){
    parallelStart <- 2

    cat("\rImputation ",1,"... ",sep="")
    if(retryFits>=1) cat("\n")
    
    fits[[1]]<-suppressMessages(fitHMM.momentuHMMData(miData[[1]],nbStates, dist, Par0[[1]], beta0[[1]], delta0[[1]],
                                    estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
                                    nlmPar, fit, DM,
                                    userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates[[1]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName))
    if(retryFits>=1){
      cat("\n")
    }
    cat("DONE\n")

    tmpPar <- getPar0(fits[[1]])
    Par0[parallelStart:nSims] <- list(tmpPar$Par)
    beta0[parallelStart:nSims] <- list(tmpPar$beta)
    delta0[parallelStart:nSims] <- list(tmpPar$delta)
  }
  
  # suppress printing of optimization progress if in parallel
  if(nSims>1 & ncores>1){
    if(!is.null(nlmPar)) nlmPar$print.level <- 0
    if(!is.null(control)) control$trace <- 0
  }
  
  if(useInitial | ncores>1) cat("Fitting ",ifelse(useInitial,"remaining ",""),length(parallelStart:nSims)," imputation",ifelse(length(parallelStart:nSims)>1,"s",""),ifelse(ncores>1," in parallel",""),"... \n",sep="")
  
  withCallingHandlers(fits[parallelStart:nSims] <-
    foreach(mD = miData[parallelStart:nSims], j = parallelStart:nSims, .export=c("fitHMM"), .errorhandling="pass", .packages=pkgs) %dorng% {
      
      if(nSims==1 | ncores==1 | retryFits>=1) {
        cat("     \rImputation ",j,"... ",sep="")
        if(retryFits>=1) cat("\n")
      } else {
        if(retryFits<1) progBar(j,nSims)
      }
      tmpFit<-suppressMessages(fitHMM(mD,nbStates, dist, Par0[[j]], beta0[[j]], delta0[[j]],
                                      estAngleMean, circularAngleMean, formula, formulaDelta, stationary, mixtures, formulaPi,
                                      nlmPar, fit, DM,
                                      userBounds, workBounds, betaCons, betaRef, deltaCons, mvnCoords, stateNames, knownStates[[j]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName))
      if(retryFits>=1) cat("\n")
      tmpFit
    } 
  ,warning=muffleRNGwarning)
  if((nSims>1 & ncores==1) | retryFits>=1) cat("DONE\n")
  
  if(ncores==1) doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  
  for(i in which(!unlist(lapply(fits,function(x) inherits(x,"momentuHMM"))))){
    warning('Fit #',i,' failed; ',fits[[i]])
  }
  
  fits <- HMMfits(fits)
  
  if(poolEstimates & nSims>1){ 
    mipool <- tryCatch(MIpool(fits,alpha=alpha,ncores=ncores,na.rm=na.rm),error=function(e) e)
    if(!inherits(mipool,"error")) fits <- miHMM(list(miSum=mipool,HMMfits=fits))
    else warning("MIpool failed: ",mipool)
  }
  
  return(fits)
}

#' @rdname MIfitHMM
#' @method MIfitHMM hierarchical
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states.  See \code{\link{fitHMM}}.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams. See \code{\link{fitHMM}}.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy. See \code{\link{fitHMM}}.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). See \code{\link{fitHMM}}.
#' 
#' @export
#' @importFrom crawl crwPostIS crwSimulator
# @importFrom raster getZ
# #' @importFrom data.tree Clone
MIfitHMM.hierarchical<-function(miData,nSims, ncores = 1, poolEstimates = TRUE, alpha = 0.95, na.rm = FALSE,
                       hierStates, hierDist, 
                       Par0, hierBeta = NULL, hierDelta = NULL,
                       estAngleMean = NULL, circularAngleMean = NULL,
                       hierFormula = NULL, hierFormulaDelta = NULL, mixtures = 1, formulaPi = NULL,
                       nlmPar = NULL, fit = TRUE, useInitial = FALSE,
                       DM = NULL, userBounds = NULL, workBounds = NULL, betaCons = NULL, deltaCons = NULL,
                       mvnCoords = NULL, knownStates = NULL, fixPar = NULL, retryFits = 0, retrySD = NULL, optMethod = "nlm", control = list(), prior = NULL, modelName = NULL,
                       covNames = NULL, spatialCovs = NULL, centers = NULL, centroids = NULL, angleCovs = NULL, altCoordNames = NULL,
                       method = "IS", parIS = 1000, dfSim = Inf, grid.eps = 1, crit = 2.5, scaleSim = 1, quad.ask = FALSE, force.quad = TRUE,
                       fullPost = TRUE, dfPostIS = Inf, scalePostIS = 1,thetaSamp = NULL, ...)
{
  
  j <- mf <- mD <- NULL #gets rid of no visible binding for global variable NOTE in R cmd check
  
  if(poolEstimates){
    if(optMethod=="nlm" & !is.null(nlmPar$hessian)){
      if(!nlmPar$hessian) stop("estimates cannot be pooled unless hessian is calculated")
    } else if(optMethod %in% fitMethods[-1] & !is.null(control$hessian)){
      if(!control$hessian) stop("estimates cannot be pooled unless hessian is calculated")
    }
    if (!requireNamespace("mitools", quietly = TRUE)) {
      stop("Package \"mitools\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  
  chkDots(...)
  
  installDataTree()
  
  if(ncores>1){
    for(pkg in c("doFuture","future")){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
    progBar <- progBar
    pkgs <- c("momentuHMM","data.tree")
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- NULL
  }
  
  if(is.crwHierData(miData)){
    
    if(nSims>=1) {
      
      model_fits <- miData$crwFits
      predData <- miData$crwPredict
      
      Time.name<-attr(predData,"Time.name")
      ids = unique(predData$ID)
      
      if(!is.null(covNames) | !is.null(angleCovs)){
        covNames <- unique(c(covNames,angleCovs[!(angleCovs %in% names(spatialCovs))]))
        if(!length(covNames)) covNames <- NULL
      }
      znames <- unique(unlist(lapply(spatialCovs,function(x) names(attributes(x)$z))))
      if(length(znames))
        if(!all(znames %in% names(predData))) stop("z-values for spatialCovs raster stack or brick not found in ",deparse(substitute(miData)),"$crwPredict")
      #coordNames <- attr(predData,"coord")
      
      if(fit | !missing("hierDist")) {
        dist <- getHierDist(hierDist,data=NULL,checkData=FALSE)
        if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
        distnames <- tmpdistnames <- names(dist)[which(!(names(dist) %in% c("step","angle")))]
        if(any(is.na(match(distnames,names(predData))))){
          for(i in which(is.na(match(distnames,names(predData))))){
            if(dist[[distnames[i]]] %in% mvndists){
              if(dist[[distnames[i]]] %in% c("mvnorm2","rw_mvnorm2")){
                tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
              } else if(dist[[distnames[i]]] %in% c("mvnorm3","rw_mvnorm3")){
                tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
              }
            }
          }
          if(any(is.na(match(tmpdistnames,names(predData))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(predData)))],collapse=", ")," not found in miData")
          tmpdistnames <- tmpdistnames[which(!(tmpdistnames %in% c("mu.x","mu.y")))]
        }
      } else {
        distnames <- tmpdistnames <- names(predData)[which(!(names(predData) %in% c("ID",Time.name,"level","locType",c("mu.x","mu.y"),covNames,znames)))]
      }
      cat('Drawing ',nSims,' realizations from the position process using crawl... \n',sep="")
      if(ncores>1) message("Running simulator in parallel... ")
      nbAnimals <- length(ids)
      predTimes <- lapply(ids,function(x) predData[[Time.name]][which(predData$ID==x & predData$locType=="p")])
      names(predTimes) <- ids
      withCallingHandlers(crwHierSim <- foreach(mf=model_fits, i=ids, .export="crwSimulator", .packages = pkgs) %dorng% {
        if(ncores==1) message("Running simulator for individual ",i,"... ",sep="")
        else progBar(which(ids==i), length(ids))
        crSim <- suppressMessages(crawl::crwSimulator(mf,predTime=predTimes[[i]], method = method, parIS = parIS,
                            df = dfSim, grid.eps = grid.eps, crit = crit, scale = scaleSim, quad.ask = ifelse(ncores>1, FALSE, quad.ask), force.quad = force.quad))
        return(crSim)
      },warning=muffleRNGwarning)
      if(ncores==1) cat("DONE\n")
      names(crwHierSim) <- ids
      
      if(ncores>1) message("Drawing imputations in parallel... ",sep="")
      withCallingHandlers(miData<-
                            foreach(j = 1:nSims, .export=c("crwPostIS","prepData"), .errorhandling="pass", .packages=pkgs) %dorng% {
                              if(ncores==1) message("\rDrawing imputation ",j,"... ",sep="")
                              else progBar(j, nSims)
                              locs<-data.frame()
                              for(i in 1:length(ids)){
                                #if(!is.null(model_fits[[i]]$err.model)){
                                tmp<-tryCatch({crawl::crwPostIS(crwHierSim[[i]], fullPost = fullPost, df = dfPostIS, scale = scalePostIS, thetaSamp = thetaSamp)},error=function(e) e)
                                if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',ids[i],'; ',tmp,'  Check crwPostIS arguments, crawl::crwMLE model fits, and/or consult crawl documentation.')
                                tlocs <- data.frame(x=tmp$alpha.sim[,"mu.x"],y=tmp$alpha.sim[,"mu.y"])
                                #tlocs[[Time.name]] <- predData[[Time.name]][which(predData$ID==ids[i] & predData$locType %in% c("o","p"))]
                                #tlocs$level <- predData$level[which(predData$ID==ids[i] & predData$locType %in% c("o","p"))]
                                #tlocs$ID <- ids[i]
                                tlocs$locType <- predData$locType[which(predData$ID==ids[i] & predData$locType %in% c("o","p"))]
                                locs<-rbind(locs,tlocs)
                              }
                              df<-predData[,c("ID",Time.name,"level","locType",tmpdistnames,covNames,znames),drop=FALSE][which(is.na(predData$locType) | predData$locType!="o"),]
                              df[which(df$locType=="p"),"x"] <- locs[which(locs$locType=="p"),"x"]
                              df[which(df$locType=="p"),"y"] <- locs[which(locs$locType=="p"),"y"]
                              df$locType <- NULL
                              pD <- tryCatch(prepData(df,covNames=covNames,spatialCovs=spatialCovs,centers=centers,centroids=centroids,angleCovs=angleCovs,altCoordNames=altCoordNames,hierLevels=levels(predData$level),coordLevel=attr(predData,"coordLevel")),error=function(e) e)
                              return(pD)
                            }
                          ,warning=muffleRNGwarning)
      if(ncores==1) cat("DONE\n")
      
      for(i in which(unlist(lapply(miData,function(x) inherits(x,"error"))))){
        warning('prepData failed for imputation ',i,"; ",miData[[i]])
      }
      ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHierHMMData"))))
      if(fit) cat('Fitting',length(ind),'realizations of the position process using fitHMM... \n')
      else {
        if(ncores==1) doParallel::stopImplicitCluster()
        else future::plan(future::sequential)
        return(crwHierSim(list(miData=miData,crwSimulator=crwHierSim)))
      }
    } else stop("nSims must be >0")
    
  } else {
    if(!is.list(miData)) stop("miData must either be a crwHierData object (as returned by crawlWrap) or a list of momentuHierHMMData objects as returned by simHierData, prepData, or MIfitHMM (when fit=FALSE)")
    if(is.crwHierSim(miData)) miData <- miData$miData
    ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHierHMMData"))))
    if(!length(ind)) stop("miData must either be a crwHierData object (as returned by crawlWrap) or a list of momentuHierHMMData objects as returned by simHierData, prepData, or MIfitHMM (when fit=FALSE)")
    if(missing(nSims)) nSims <- length(miData)
    if(nSims>length(miData)) stop("nSims is greater than the length of miData. nSims must be <=",length(miData))
    if(nSims<1) stop("nSims must be >0")
    cat('Fitting',min(nSims,length(ind)),'imputation(s) using fitHMM... \n')
  }
  
  nbStates <- nbHierStates(hierStates)$nbStates
  dist <- getHierDist(hierDist,data=NULL,checkData=FALSE)
  
  if(!is.list(knownStates)){
    tmpStates<-knownStates
    knownStates<-vector('list',nSims)
    if(!is.null(tmpStates))
      knownStates[1:nSims]<-list(tmpStates)
  } else if(length(knownStates)<nSims) stop("knownStates must be a list of length >=",nSims)
  
  if(all(names(dist) %in% names(Par0))){
    tmpPar0<-Par0
    Par0<-vector('list',nSims)
    Par0[1:nSims]<-list(tmpPar0)
  } else if(length(Par0)<nSims) stop("Par0 must be a list of length >=",nSims)
  
  #newForm <- newFormulas(formula,nbStates,betaRef)
  #recharge <- newForm$recharge
  
  if(mixtures==1){
    if(!is.list(hierBeta)){
      if(!inherits(hierBeta,"Node"))
        tmphierBeta<-hierBeta
      else tmphierBeta <- data.tree::Clone(hierBeta)
      hierBeta<-vector('list',nSims)
      if(!is.null(tmphierBeta))
        hierBeta[1:nSims]<-list(tmphierBeta)
    } else if(length(hierBeta)<nSims) stop("hierBeta must be a list of length >=",nSims)
  } else {
    if(!is.null(hierBeta) && !is.list(hierBeta)){
      #if(!is.null(recharge)) stop("hierBeta must be a list with elements named 'beta', 'g0', and/or 'theta' when a recharge model is specified")
      if(mixtures>1) stop("hierBeta must be a list with elements named 'beta' and/or 'pi' when mixtures>1")
    }
    if(!is.list(hierBeta[[1]])){
      if(inherits(hierBeta$beta,"Node")){
        tmphierBeta <- vector('list',length(hierBeta))
        names(tmphierBeta) <- names(hierBeta)
        tmphierBeta[names(hierBeta)!="beta"] <- hierBeta[names(hierBeta)!="beta"]
        tmphierBeta$beta <- data.tree::Clone(hierBeta$beta)
        hierBeta<-vector('list',nSims)
        hierBeta[1:nSims]<-list(tmphierBeta)
      } else {
        tmphierBeta<-hierBeta
        hierBeta<-vector('list',nSims)
        if(!is.null(tmphierBeta)){
          hierBeta[1:nSims]<-list(tmphierBeta)
        }
      }
    } else if(length(hierBeta)<nSims) stop("hierBeta must be a list of length >=",nSims)     
  }
  if(!is.list(hierDelta)){
    if(inherits(hierDelta,"Node")){
      tmphierDelta<-data.tree::Clone(hierDelta)
      hierDelta<-vector('list',nSims)
      hierDelta[1:nSims]<-list(tmphierDelta)
    } else {
      tmphierDelta<-hierDelta
      hierDelta<-vector('list',nSims)
      if(!is.null(tmphierDelta))
        hierDelta[1:nSims]<-list(tmphierDelta)
    }
  } else if(length(hierDelta)<nSims) stop("hierDelta must be a list of length >=",nSims)
  
  #check HMM inputs and print model message
  test<-fitHMM.momentuHierHMMData(miData[[ind[1]]], hierStates, hierDist, Par0[[ind[1]]], hierBeta[[ind[1]]], hierDelta[[ind[1]]],
                   estAngleMean, circularAngleMean, hierFormula, hierFormulaDelta, mixtures, formulaPi, 
                   nlmPar, fit = FALSE, DM,
                   userBounds, workBounds, betaCons, deltaCons, mvnCoords, knownStates[[ind[1]]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName)
  
  # fit HMM(s)
  fits <- list()
  parallelStart <- 1
  useInitial <- ifelse(nSims<2,FALSE,useInitial)
  if(useInitial){
    parallelStart <- 2

    cat("\rImputation ",1,"... ",sep="")

    fits[[1]]<-suppressMessages(fitHMM.momentuHierHMMData(miData[[1]], hierStates, hierDist, Par0[[1]], hierBeta[[1]], hierDelta[[1]],
                                           estAngleMean, circularAngleMean, hierFormula, hierFormulaDelta, mixtures, formulaPi, 
                                           nlmPar, fit, DM,
                                           userBounds, workBounds, betaCons, deltaCons, mvnCoords, knownStates[[1]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName))
    if(retryFits>=1){
      cat("\n")
    }
    cat("DONE\nFitting remaining imputations... \n")
    
    tmpPar <- getPar0(fits[[1]])
    Par0[parallelStart:nSims] <- list(tmpPar$Par)
    beta0[parallelStart:nSims] <- list(tmpPar$beta)
    delta0[parallelStart:nSims] <- list(tmpPar$delta)
  }
  
  # suppress printing of optimization progress if in parallel
  if(nSims>1 & ncores>1){
    if(!is.null(nlmPar)) nlmPar$print.level <- 0
    if(!is.null(control)) control$trace <- 0
    message("        Fitting ",length(parallelStart:nSims)," imputation",ifelse(length(parallelStart:nSims)>1,"s","")," in parallel... ",sep="")
  }
  
  withCallingHandlers(fits[parallelStart:nSims] <-
                        foreach(mD = miData[parallelStart:nSims], j = parallelStart:nSims, .export=c("fitHMM"), .errorhandling="pass", .packages = pkgs) %dorng% {
                          
                          if(nSims==1 & ncores==1) {
                            cat("     \rImputation ",j,"... ",sep="")
                          } else {
                            progBar(j,nSims)
                          }
                          tmpFit<-suppressMessages(fitHMM(mD, hierStates, hierDist, Par0[[j]], hierBeta[[j]], hierDelta[[j]],
                                                              estAngleMean, circularAngleMean, hierFormula, hierFormulaDelta, mixtures, formulaPi, 
                                                              nlmPar, fit, DM, 
                                                              userBounds, workBounds, betaCons, deltaCons, mvnCoords, knownStates[[j]], fixPar, retryFits, retrySD, optMethod, control, prior, modelName))
                          if(retryFits>=1 & ncores==1) cat("\n")
                          tmpFit
                        } 
                      ,warning=muffleRNGwarning)
  if(nSims>1 & ncores==1) cat("DONE\n")
  
  if(ncores==1) doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  
  for(i in which(!unlist(lapply(fits,function(x) inherits(x,"momentuHierHMM"))))){
    warning('Fit #',i,' failed; ',fits[[i]])
  }
  
  fits <- HMMfits(fits)
  class(fits) <- append(class(fits),"hierarchical")
  
  if(poolEstimates & nSims>1){ 
    mipool <- tryCatch(MIpool(fits,alpha=alpha,ncores=ncores,na.rm=na.rm),error=function(e) e)
    if(!inherits(mipool,"error")) {
      fits <- miHMM(list(miSum=mipool,HMMfits=fits))
      class(fits) <- append(class(fits),"hierarchical")
    } else warning("MIpool failed: ",mipool)
  }
  
  return(fits)
}
