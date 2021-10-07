
#' @rdname simCTHMM
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states ('state').  See details.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams ('dist'). Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. See details.
#' @param hierBeta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the transition probabilities at each level of the hierarchy ('beta'). See \code{\link{fitHMM}}. 
#' @param hierDelta A hierarchical data structure \code{\link[data.tree]{Node}} for the matrix of initial values for the regression coefficients of the initial distribution at each level of the hierarchy ('delta'). See \code{\link{fitHMM}}. 
#' @param hierFormula A hierarchical formula structure for the transition probability covariates for each level of the hierarchy ('formula'). Default: \code{NULL} (only hierarchical-level effects, with no covariate effects).
#' Any formula terms that are not state- or parameter-specific are included on all of the transition probabilities within a given level of the hierarchy. See details.
#' @param hierFormulaDelta A hierarchical formula structure for the initial distribution covariates for each level of the hierarchy ('formulaDelta'). Default: \code{NULL} (no covariate effects and \code{fixPar$delta} is specified on the working scale). 
#' @param nbHierCovs A hierarchical data structure \code{\link[data.tree]{Node}} for the number of covariates ('nbCovs') to simulate for each level of the hierarchy (0 by default). Does not need to be specified if
#' \code{covs} is specified. Simulated covariates are provided generic names (e.g., 'cov1.1' and 'cov1.2' for \code{nbHierCovs$level1$nbCovs=2}) and can be included in \code{hierFormula} and/or \code{DM}.
#' @param obsPerLevel A hierarchical data structure \code{\link[data.tree]{Node}} indicating the number of observations for each level of the hierarchy ('obs'). For each level, the 'obs' field can either be the number of observations per animal (if single value) or the bounds of the number of observations per animal (if vector of two values). In the latter case, 
#' the numbers of obervations generated per level for each animal are uniformously picked from this interval. Alternatively, \code{obsPerLevel} can be specified as
#' a list of length \code{nbAnimals} with each element providing the hierarchical data structure for the number of observations for each level of the hierarchy for each animal, where the 'obs' field can either be the number of observations (if single value) or the bounds of the number of observations (if vector of two values) for each individual.
#'
#' @details \itemize{
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal' (or 'obsPerLevel' for \code{simHierCTHMM}), then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' 
#' \item For \code{simCTHMM}, when covariates are not included in \code{formulaDelta} (i.e. \code{formulaDelta=NULL}), then \code{delta} is specified as a vector of length \code{nbStates} that 
#' sums to 1.  When covariates are included in \code{formulaDelta}, then \code{delta} must be specified
#' as a k x (\code{nbStates}-1) matrix of working parameters, where k is the number of regression coefficients and the columns correspond to states 2:\code{nbStates}. For example, in a 3-state
#' HMM with \code{formulaDelta=~cov1+cov2}, the matrix \code{delta} has three rows (intercept + two covariates)
#' and 2 columns (corresponding to states 2 and 3). The initial distribution working parameters are transformed to the real scale as \code{exp(covsDelta*Delta)/rowSums(exp(covsDelta*Delta))}, where \code{covsDelta} is the N x k design matrix, \code{Delta=cbind(rep(0,k),delta)} is a k x \code{nbStates} matrix of working parameters,
#' and \code{N=length(unique(data$ID))}.
#' 
#' \item For \code{simHierCTHMM}, \code{delta} must be specified
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
# @importFrom raster cellFromXY getValues
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
                        export=NULL)
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
                                         CT=TRUE),warning=muffleCTwarning)
  
  coordLevel <- attr(out,"coordLevel")
  out<- out[,c("ID","time",colnames(out)[which(!colnames(out) %in% c("ID","time"))])]
  if(!is.null(mvnCoords)){
    attr(out,'coords') <- paste0(mvnCoords,c(".x",".y"))
    attr(out,"coordLevel") <- coordLevel
  }
  attr(out,"CT") <- TRUE
  attr(out,"Time.name") <- "time"
  return(out)
}
