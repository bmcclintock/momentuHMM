
#' AIC
#'
#' Akaike information criterion of momentuHMM model(s).
#'
#' @method AIC momentuHMM
#'
#' @param object A \code{momentuHMM} object.
#' @param ... Optional additional \code{momentuHMM} objects, to compare AICs of the different models. These can be passed as a list using the \code{!!!} operator (see \code{\link[rlang]{splice-operator}} and example in \code{\link{AICweights}}).
#' @param k Penalty per parameter. Default: 2 ; for classical AIC.
#' @param n Optional sample size. If specified, the small sample correction AIC is used (i.e., \code{AICc = AIC + kp(p+1)/(n-p-1)} where p is the number of parameters).
#'
#' @return The AIC of the model(s) provided. If several models are provided, the AICs are output
#' in ascending order.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' AIC(m)
#' 
#' \dontrun{
#' # HMM specifications
#' nbStates <- 2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar0 <- c(mu0,sigma0)
#' anglePar0 <- c(-pi/2,pi/2,kappa0)
#' formula <- ~cov1+cov2
#'                 
#' # example$m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' mod1 <- fitHMM(example$m$data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                 Par0=list(step=stepPar0,angle=anglePar0),
#'                 formula=~1,estAngleMean=list(angle=TRUE))
#' 
#' Par0 <- getPar0(mod1,formula=formula)                 
#' mod2 <- fitHMM(example$m$data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                 Par0=Par0$Par,beta0=Par0$beta,
#'                 formula=formula,estAngleMean=list(angle=TRUE))
#'                 
#' AIC(mod1,mod2)
#' 
#' Par0nA <- getPar0(mod1,estAngleMean=list(angle=FALSE))                 
#' mod3 <- fitHMM(example$m$data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                 Par0=Par0nA$Par,beta0=Par0nA$beta,
#'                 formula=~1)
#'  
#' AIC(mod1,mod2,mod3)
#' 
#' # add'l models provided as a list using the !!! operator                              
#' AIC(mod1, !!!list(mod2,mod3))
#' }
#'
#' @importFrom rlang list2
#' @export

AIC.momentuHMM <- function(object,...,k=2,n=NULL)
{
  models <- rlang::list2(...)

  if(length(models)>0) { # if several models are provided
    modNames <- all.vars(match.call()) # store the names of the models given as arguments

    # include "object" in "models"
    modcopy <- list()
    modcopy[[1]] <- object
    for(i in 1:length(models))
      modcopy[[i+1]] <- models[[i]]
    models <- modcopy
    
    for(i in 1:length(models)){
      if(!is.null(models[[i]]$modelName)) modNames[i] <- models[[i]]$modelName
    }
    
    if(any(!unlist(lapply(models,is.momentuHMM)))) stop("all models must be momentuHMM objects")
    
    pr <- is.null(models[[1]]$prior)
    for(i in 2:length(models)) {
      datNames <- sort(colnames(models[[1]]$data)[colnames(models[[1]]$data) %in% colnames(models[[i]]$data)])
      if(!length(datNames)) stop("data must be the same for each momentuHMM object")
      if(!isTRUE(all.equal(models[[i]]$data[,datNames],models[[1]]$data[,datNames]))) stop("data must be the same for each momentuHMM object")
      if(pr!=is.null(models[[i]]$prior)) stop("AIC is not valid for comparing models with and without priors")
    }
    if(!pr) warning("Please be careful when using AIC to compare models with priors!")

    # compute AICs of models
    aic <- rep(NA,length(models))

    for(i in 1:length(models)) {
      if(!inherits(models[[i]],"randomEffects")) aic[i] <- getAIC(models[[i]],k,n)
      else {
        if(!is.null(models[[i]]$conditions$fit) && !models[[i]]$conditions$fit) stop("The given model hasn't been fitted.")
        nbPar <- sum(unlist(models[[i]]$traceG))+length(models[[i]]$mod$wpar)
        maxLogLike <- -models[[i]]$mod$minimum
        aic[i] <- -2*maxLogLike+k*nbPar
        if(!is.null(n)) aic[i] <- aic[i] + k*nbPar*(nbPar+1)/(n-nbPar-1)
      }
    }

    ord <- order(aic) # order models by increasing AIC
    return(data.frame(Model=modNames[ord],AIC=aic[ord]))
  }
  else { # if only one model is provided
    if(!inherits(object,"randomEffects")) aic <- getAIC(object,k,n)
    else {
      if(!is.null(object$conditions$fit) && !object$conditions$fit) stop("The given model hasn't been fitted.")
      nbPar <- sum(unlist(object$traceG))+length(object$mod$wpar)
      maxLogLike <- -object$mod$minimum
      aic <- -2*maxLogLike+k*nbPar
      if(!is.null(n)) aic <- aic + k*nbPar*(nbPar+1)/(n-nbPar-1)
    }
    return(aic)
  }
}

getAIC <- function(m, k=2, n=NULL){
    #nbPar <- length(m$mod$estimate)-sum(!is.na(unlist(m$conditions$fixPar)))-(length(m$conditions$betaCons[is.na(m$conditions$fixPar$beta)])-length(unique(m$conditions$betaCons[is.na(m$conditions$fixPar$beta)])))
    #if(length(m$conditions$fixPar$delta)==length(m$stateNames) & all(!is.na(m$conditions$fixPar$delta))) nbPar <- nbPar + 1
    if(!is.null(m$conditions$fit) && !m$conditions$fit) stop("The given model hasn't been fitted.")
    else nbPar <- length(m$mod$wpar)
    maxLogLike <- -m$mod$minimum
    aic <- -2*maxLogLike+k*nbPar
    if(!is.null(n)) aic <- aic + k*nbPar*(nbPar+1)/(n-nbPar-1)
    aic
}