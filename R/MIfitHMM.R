
#' Fit HMMs to multiple imputation data
#' 
#' Fit a (multivariate) hidden Markov model to multiple imputation data. Multiple imputation is a method for accommodating 
#' missing data, temporal-irregularity, or location measurement error in hidden Markov models, where pooled parameter estimates reflect uncertainty
#' attributable to observation error.
#' 
#' \code{miData} can either be a \code{\link{crwData}} object (as returned by \code{\link{crawlWrap}}), a \code{\link{crwSim}} object (as returned by \code{MIfitHMM} when \code{fit=FALSE}), 
#' or a list of \code{\link{momentuHMMData}} objects (e.g., each element of the list as returned by \code{\link{prepData}}). 
#' 
#' If \code{miData} is a \code{crwData} object, \code{MIfitHMM} uses a combination of 
#' \code{\link[crawl]{crwSimulator}}, \code{\link[crawl]{crwPostIS}}, \code{\link{prepData}}, and \code{\link{fitHMM}} to draw \code{nSims} realizations of the position process
#' and fit the specified HMM to each imputation of the data. The vast majority of \code{MIfitHMM} arguments are identical to the corresponding arguments from these functions.
#' 
#' If \code{miData} is a \code{\link{crwData}} object, \code{nSims} determines both the number of realizations of the position process to draw 
#' (using \code{\link[crawl]{crwSimulator}} and \code{\link{crwPostIS}}) as well as the number of HMM fits.
#' 
#' If \code{miData} is a \code{\link{crwSim}} object or a list of \code{\link{momentuHMMData}} object(s), the specified HMM will simply be fitted to each of the \code{momentuHMMData} objects
#' and all arguments related to \code{\link[crawl]{crwSimulator}}, \code{\link[crawl]{crwPostIS}}, or \code{\link{prepData}} are ignored.
#' 
#' @param miData A \code{\link{crwData}} object, a \code{\link{crwSim}} object, or a list of \code{\link{momentuHMMData}} objects.
#' @param nSims Number of imputations in which to fit the HMM using \code{\link{fitHMM}}. If \code{miData} is a list of \code{momentuHMMData} 
#' objects, \code{nSims} cannot exceed the length of \code{miData}.
#' @param ncores Number of cores to use for parallel processing.
#' @param poolEstimates Logical indicating whether or not to calculate pooled parameter estimates across the \code{nSims} imputations using \code{\link{MIpool}}. Default: \code{TRUE}.
#' @param alpha Significance level for calculating confidence intervals of pooled estimates when \code{poolEstimates=TRUE} (see \code{\link{MIpool}}). Default: 0.95.
#' @param nbStates Number of states of the HMM. See \code{\link{fitHMM}}.
#' @param dist A named list indicating the probability distributions of the data streams. See \code{\link{fitHMM}}.
#' @param Par0 A named list containing vectors of initial state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. See \code{\link{fitHMM}}.  \code{Par0} may also be a list of length \code{nSims}, where each element is a named list containing vectors
#' of initial state-dependent probability distribution parameters for each imputation.
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities. See \code{\link{fitHMM}}. \code{beta0} may also be a list of length \code{nSims}, where each element 
#' is an initial matrix of regression coefficients for the transition probabilities for each imputation.
#' @param delta0 Initial values for the initial distribution of the HMM. See \code{\link{fitHMM}}. \code{delta0} may also be a list of length \code{nSims}, where each element 
#' is the initial values for the initial distribution of the HMM for each imputation.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). See \code{\link{fitHMM}}.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear (FALSE) or circular-circular (TRUE) 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles. See \code{\link{fitHMM}}.
#' @param formula Regression formula for the transition probability covariates. See \code{\link{fitHMM}}.
#' @param formulaDelta Regression formula for the initial distribution. See \code{\link{fitHMM}}.
#' @param stationary \code{FALSE} if there are covariates. If \code{TRUE}, the initial distribution is considered
#' equal to the stationary distribution. See \code{\link{fitHMM}}.
#' @param verbose Determines the print level of the \code{fitHMM} optimizer. The default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed.
#' @param nlmPar List of parameters to pass to the \code{fitHMM} optimization function \code{nlm} (which should be either
#' '\code{gradtol}', '\code{stepmax}', '\code{steptol}', or '\code{iterlim}' -- see \code{nlm}'s documentation
#' for more detail)
#' @param fit \code{TRUE} if the HMM should be fitted to the data, \code{FALSE} otherwise. See \code{\link{fitHMM}}. If \code{fit=FALSE} and \code{miData} is a \code{\link{crwData}}
#' object, then \code{MIfitHMM} returns a list containing a \code{\link{momentuHMMData}} object (if \code{nSims=1}) or, if \code{nSims>1}, a \code{\link{crwSim}} object.
#' @param DM An optional named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. See \code{\link{fitHMM}}.
#' @param cons An optional named list of vectors specifying a power to raise parameters corresponding to each column of the design matrix 
#' for each data stream. See \code{\link{fitHMM}}.
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. See \code{\link{fitHMM}}.
#' @param workcons An optional named list of vectors specifying constants to add to the regression coefficients on the working scale for 
#' each data stream. See \code{\link{fitHMM}}.
#' @param stateNames Optional character vector of length nbStates indicating state names.
#' @param knownStates Vector of values of the state process which are known prior to fitting the
#' model (if any). See \code{\link{fitHMM}}. If \code{miData} is a list of \code{\link{momentuHMMData}} objects, then \code{knownStates} can alternatively
#' be a list of vectors containing the known values for the state process for each element of \code{miData}.
#' @param fixPar An optional list of vectors indicating parameters which are assumed known prior to fitting the model. See \code{\link{fitHMM}}. 
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization.  See \code{\link{fitHMM}}. 
#' @param covNames Names of any covariates in \code{miData$crwPredict} (if \code{miData} is a \code{\link{crwData}} object; otherwise 
#' \code{covNames} is ignored). See \code{\link{prepData}}. 
#' @param spatialCovs List of raster layer(s) for any spatial covariates not included in \code{miData$crwPredict} (if \code{miData} is 
#' a \code{\link{crwData}} object; otherwise \code{spatialCovs} is ignored). See \code{\link{prepData}}. 
#' @param centers 2-column matrix providing the x-coordinates (column 1) and y-coordinates (column 2) for any activity centers (e.g., potential 
#' centers of attraction or repulsion) from which distance and angle covariates will be calculated based on realizations of the position process. 
#' See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param angleCovs Character vector indicating the names of any circular-circular regression angular covariates in \code{miData$crwPredict} that need conversion from standard direction (in radians relative to the x-axis) to turning angle (relative to previous movement direction) 
#' See \code{\link{prepData}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param method Method for obtaining weights for movement parameter samples. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param parIS Size of the parameter importance sample. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param dfSim Degrees of freedom for the t approximation to the parameter posterior. See 'df' argument in \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param grid.eps Grid size for \code{method="quadrature"}. See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param crit Criterion for deciding "significance" of quadrature points
#' (difference in log-likelihood). See \code{\link[crawl]{crwSimulator}}. Ignored unless \code{miData} is a \code{\link{crwData}} object.
#' @param scaleSim Scale multiplier for the covariance matrix of the t approximation. See 'scale' argument in \code{\link[crawl]{crwSimulator}}. 
#' Ignored unless \code{miData} is a \code{\link{crwData}} object.
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
#' 
#' # Don't run because it takes too long on a single core
#' \dontrun{
#' # extract simulated obsData from example data
#' obsData <- miExample$obsData
#' 
#' # extract crwMLE inputs from example data
#' inits <- miExample$inits # initial state
#' err.model <- miExample$err.model # error ellipse model
#'
#' # create crwData object by fitting crwMLE models to obsData and predict locations 
#' # at default intervals for both individuals
#' crwOut <- crawlWrap(obsData=obsData,ncores=1,
#'          theta=c(4,0),fixPar=c(1,1,NA,NA),
#'          initial.state=inits,
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
#' miFits<-MIfitHMM(miData=crwOut,nSims=5,ncores=1,
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
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom raster getZ
MIfitHMM<-function(miData,nSims, ncores, poolEstimates = TRUE, alpha = 0.95,
                   nbStates, dist, 
                   Par0, beta0 = NULL, delta0 = NULL,
                   estAngleMean = NULL, circularAngleMean = NULL,
                   formula = ~1, formulaDelta = ~1, stationary = FALSE, 
                   verbose = 0, nlmPar = NULL, fit = TRUE, 
                   DM = NULL, cons = NULL, userBounds = NULL, workcons = NULL, 
                   stateNames = NULL, knownStates = NULL, fixPar = NULL, retryFits = 0,
                   covNames=NULL,spatialCovs=NULL,centers=NULL,angleCovs=NULL,
                   method = "IS", parIS = 1000, dfSim = Inf, grid.eps = 1, crit = 2.5, scaleSim = 1, force.quad = TRUE,
                   fullPost = TRUE, dfPostIS = Inf, scalePostIS = 1,thetaSamp = NULL
                   ){

  j <- NULL #gets rid of no visible binding for global variable 'j' NOTE in R cmd check
  
  if(is.crwData(miData)){
    
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
      
      if(fit | !missing("dist")) {
        if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
        distnames <- names(dist)[which(!(names(dist) %in% c("step","angle")))]
        if(any(is.na(match(distnames,names(predData))))) stop(paste0(distnames[is.na(match(distnames,names(predData)))],collapse=", ")," not found in miData")
      } else {
        distnames <- names(predData)[which(!(names(predData) %in% c("ID","x","y",covNames,znames)))]
      }
      
      if(all(unlist(lapply(model_fits,function(x) is.null(x$err.model))))) stop("Multiple realizations of the position process cannot be drawn if there is no location measurement error")
      cat('Drawing',nSims,'realizations from the position process using crawl::crwPostIS... ')
      
      registerDoParallel(cores=ncores)
      crwSim <- foreach(i = 1:length(ids), .export="crwSimulator") %dopar% {
        if(!is.null(model_fits[[i]]$err.model))
          crawl::crwSimulator(model_fits[[i]],predTime=predData[[Time.name]][which(predData$ID==ids[i] & predData$locType=="p")], method = method, parIS = parIS,
                                  df = dfSim, grid.eps = grid.eps, crit = crit, scale = scaleSim, force.quad = force.quad)
      }
      stopImplicitCluster()
      
      registerDoParallel(cores=ncores)
      miData<-
        foreach(j = 1:nSims, .export=c("crwPostIS","prepData"), .errorhandling="pass") %dopar% {
          locs<-data.frame()
          for(i in 1:length(ids)){
            if(!is.null(model_fits[[i]]$err.model)){
              tmp<-tryCatch({crawl::crwPostIS(crwSim[[i]], fullPost = fullPost, df = dfPostIS, scale = scalePostIS, thetaSamp = thetaSamp)},error=function(e) e)
              if(!all(class(tmp) %in% c("crwIS","list"))) stop('crawl::crwPostIS error for individual ',ids[i],'; ',tmp,'  Check crwPostIS arguments, crawl::crwMLE model fits, and/or consult crawl documentation.')
              locs<-rbind(locs,tmp$alpha.sim[,c("mu.x","mu.y")])
            } else {
              locs<-rbind(locs,predData[which(predData$ID==ids[i]),c("mu.x","mu.y")])
            }
          }
          df<-data.frame(x=locs$mu.x,y=locs$mu.y,predData[,c("ID",distnames,covNames,znames),drop=FALSE])[which(predData$locType=="p"),]
          prepData(df,covNames=covNames,spatialCovs=spatialCovs,centers=centers,angleCovs=angleCovs)
        }
      stopImplicitCluster()
      cat("DONE\n")
      for(i in which(unlist(lapply(miData,function(x) inherits(x,"error"))))){
        warning('prepData failed for imputation ',i,"; ",miData[[i]])
      }
      ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHMMData"))))
      if(fit) cat('Fitting',length(ind),'realizations of the position process using fitHMM... ')
      else return(crwSim(list(miData=miData,crwSimulator=crwSim)))
    } else stop("nSims must be >0")
    
  } else {
    if(!is.list(miData)) stop("miData must either be a crwData object (as returned by crawlWrap) or a list of momentuHMMData objects as returned by simData, prepData, or MIfitHMM (when fit=FALSE)")
    if(is.crwSim(miData)) miData <- miData$miData
    ind <- which(unlist(lapply(miData,function(x) inherits(x,"momentuHMMData"))))
    if(!length(ind)) stop("miData must either be a crwData object (as returned by crawlWrap) or a list of momentuHMMData objects as returned by simData, prepData, or MIfitHMM (when fit=FALSE)")
    if(missing(nSims)) nSims <- length(miData)
    if(nSims>length(miData)) stop("nSims is greater than the length of miData. nSims must be <=",length(miData))
    if(nSims<1) stop("nSims must be >0")
    cat('Fitting',min(nSims,length(ind)),'imputation(s) using fitHMM... ')
  }
  
  if(!is.list(knownStates)){
    tmpStates<-knownStates
    knownStates<-vector('list',nSims)
    if(!is.null(tmpStates))
      for(i in 1:nSims)
        knownStates[[i]]<-tmpStates
  } else if(length(knownStates)<nSims) stop("knownStates must be a list of length >=",nSims)
  
  if(all(names(dist) %in% names(Par0))){
    tmpPar0<-Par0
    Par0<-vector('list',nSims)
    for(i in 1:nSims)
      Par0[[i]]<-tmpPar0
  } else if(length(Par0)<nSims) stop("Par0 must be a list of length >=",nSims)
  
  if(!is.list(beta0)){
    tmpbeta0<-beta0
    beta0<-vector('list',nSims)
    if(!is.null(tmpbeta0))
      for(i in 1:nSims)
        beta0[[i]]<-tmpbeta0
  } else if(length(beta0)<nSims) stop("beta0 must be a list of length >=",nSims)
  
  if(!is.list(delta0)){
    tmpdelta0<-delta0
    delta0<-vector('list',nSims)
    if(!is.null(tmpdelta0))
      for(i in 1:nSims)
        delta0[[i]]<-tmpdelta0
  } else if(length(delta0)!=nSims) stop("delta0 must be a list of length ",nSims)
  
  #check HMM inputs and print model message
  test<-fitHMM(miData[[ind[1]]],nbStates, dist, Par0[[ind[1]]], beta0[[ind[1]]], delta0[[ind[1]]],
         estAngleMean, circularAngleMean, formula, formulaDelta, stationary, verbose,
         nlmPar, fit = FALSE, DM, cons,
         userBounds, workcons, stateNames, knownStates[[ind[1]]], fixPar, retryFits)
  
  # fit HMM(s)
  registerDoParallel(cores=ncores)
  fits <-
    foreach(j = 1:nSims, .export=c("fitHMM"), .errorhandling="pass") %dopar% {

      if(nSims>1) cat("\rImputation ",j,"...",sep="")
      tmpFit<-suppressMessages(fitHMM(miData[[j]],nbStates, dist, Par0[[j]], beta0[[j]], delta0[[j]],
                              estAngleMean, circularAngleMean, formula, formulaDelta, stationary, verbose,
                              nlmPar, fit, DM, cons,
                              userBounds, workcons, stateNames, knownStates[[j]], fixPar, retryFits))
      if(retryFits>=1) cat("\n")
      tmpFit
    }  
  stopImplicitCluster()
  cat("DONE\n")
  
  for(i in which(!unlist(lapply(fits,function(x) inherits(x,"momentuHMM"))))){
    warning('Fit #',i,' failed; ',fits[[i]])
  }
  
  if(nSims==1) out<-fits[[1]]
  else {
    if(poolEstimates) out<-miHMM(list(miSum=MIpool(fits,alpha=alpha,ncores=ncores),HMMfits=fits))
    else out <- fits
  }
  out
}
