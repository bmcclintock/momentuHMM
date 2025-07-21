
#' Calculate Akaike information criterion model weights 
#'
#' @param ... \code{\link{momentuHMM}}, \code{\link{HMMfits}}, or \code{\link{miHMM}} objects, to compare AIC weights of the different models. The first object must be a \code{\link{momentuHMM}}, \code{\link{HMMfits}}, or \code{\link{miHMM}} object, but additional model objects can be passed as a list using the \code{!!!} operator (see \code{\link[rlang]{splice-operator}}).
#' @param k Penalty per parameter. Default: 2 ; for classical AIC.
#' @param n Optional sample size. If specified, the small sample correction AIC is used (i.e., \code{AICc = AIC + kp(p+1)/(n-p-1)} where p is the number of parameters).
#'
#' @return The AIC weights of the models. If multiple imputation objects are provided, 
#' then the mean model weights (and standard deviations) are provided.
#' 
#' @details
#' \itemize{
#' \item Model objects must all be either of class \code{\link{momentuHMM}} or multiple imputation model objects (of class \code{\link{HMMfits}} and/or \code{\link{miHMM}}).
#' 
#' \item AIC is only valid for comparing models fitted to the same data. The data for each model fit must therefore be identical. For multiple imputation model objects, 
#' respective model fits must have identical data.}
#' 
#' @examples
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
#' AICweights(mod1,mod2)
#' 
#' Par0nA <- getPar0(mod1,estAngleMean=list(angle=FALSE))                 
#' mod3 <- fitHMM(example$m$data,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
#'                 Par0=Par0nA$Par,beta0=Par0nA$beta,
#'                 formula=~1)
#'  
#' AICweights(mod1,mod2,mod3)
#' 
#' # add'l models provided as a list using the !!! operator                               
#' AICweights(mod1, !!!list(mod2,mod3))
#' }
#' 
#' @export
#' 

AICweights <- function (..., k=2, n=NULL) {
  UseMethod("AICweights")
}

#' @method AICweights momentuHMM
#' @importFrom rlang list2
#' @export
AICweights.momentuHMM <- function(..., k=2, n=NULL)
{
  models <- rlang::list2(...)
  
  modNames <- all.vars(match.call()) # store the names of the models given as arguments
  
  if(any(!unlist(lapply(models,is.momentuHMM)))) stop("all models must be momentuHMM objects")
  
  if(length(models)<2) stop("at least 2 momentuHMM objects must be provided")
  
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
    aic[i] <- AIC.momentuHMM(models[[i]],k=k,n=n)
  }
  
  for(i in 1:length(models)){
    if(!is.null(models[[i]]$modelName)) modNames[i] <- models[[i]]$modelName
  }
  
  ord <- order(aic) # order models by increasing AIC
  weight <- exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
  
  return(data.frame(Model=modNames[ord],weight=weight[ord]))
  
}

#' @method AICweights miHMM
#' @importFrom rlang list2
#' @export
AICweights.miHMM <- function(...,k=2, n=NULL)
{
  models <- rlang::list2(...)
  
  modNames <- all.vars(match.call()) # store the names of the models given as arguments
  
  if(any(!unlist(lapply(models,is.miHMM)) & !unlist(lapply(models,is.HMMfits)))) stop("all model objects must be of class miHMM or HMMfits")
  
  for(i in which(unlist(lapply(models,is.miHMM)))){
    models[[i]] <- models[[i]]$HMMfits
  }
  if(length(models)<2) stop("at least 2 model objects must be provided")
  
  nSims <- unique(unlist(lapply(models,length)))
  if(length(nSims)>1){
    nSims <- min(nSims)
    warning("HMMfits objects are of different lengths; calculating AIC weights based on first ",nSims," fits for each model")
  }
  
  momObs <- which(!unlist(lapply(models[[1]][1:nSims],is.momentuHMM)))
  pr <- logical(nSims)
  for(i in 2:length(models)) {
    for(j in 1:nSims){
      pr[j] <- (!is.null(models[[1]][[j]]$prior))
      if(!(j %in% momObs)){
        if(is.momentuHMM(models[[i]][[j]])){
          datNames <- sort(colnames(models[[1]][[j]]$data)[colnames(models[[1]][[j]]$data) %in% colnames(models[[i]][[j]]$data)])
          if(!length(datNames)) stop("Imputed data must be the same for each model object")
          if(!isTRUE(all.equal(models[[i]][[j]]$data[,datNames],models[[1]][[j]]$data[,datNames]))) stop("Imputed data must be the same for each model object")
          if(pr[j]!=(!is.null(models[[i]][[j]]$prior))) stop("AIC is not valid for comparing models with and without priors")
        } else {
          momObs <- c(momObs,j)
        }
      }
    }
  }
  if(any(pr)) warning("Please be careful when using AIC to compare models with priors!")
  
  if(length(momObs)) {
    momObs <- sort(momObs)
    iSims <- (1:nSims)[-momObs]
    warning("Some model fits are not momentuHMM objects. Imputations ",paste(momObs,collapse=", ")," will be ignored")
  } else iSims <- 1:nSims
  
  # compute AICs of models
  aic <- matrix(NA,length(models),nSims)
  
  for(i in 1:length(models)) {
    # check modelName
    checkNames <- lapply(models[[i]],function(x) x[match("modelName",names(x))])
    if(any(!unlist(lapply(checkNames,function(x) isTRUE(all.equal(x,checkNames[[1]],use.names=FALSE)))))) stop("'modelName' must be identical for fitted models within each miHMM or HMMfits object")
    for(j in iSims){
      aic[i,j] <- AIC.momentuHMM(models[[i]][[j]],k=k,n=n)
    }
    if(!is.null(models[[i]][[1]]$modelName)) modNames[i] <- models[[i]][[1]]$modelName
  }
  
  weights <- apply(aic,2,function(x) exp(-0.5*(x-min(x)))/sum(exp(-0.5*(x-min(x)))))
  weight <- apply(weights,1,mean,na.rm=TRUE)
  #se <- sqrt((nSims+1)/(nSims*(nSims-1))*rowSums((weights-weight)^2))
  sd <- apply(weights,1,sd,na.rm=TRUE)
  ord <- order(weight,decreasing = TRUE) # order models by increasing AIC
  
  return(data.frame(Model=modNames[ord],weight=weight[ord],sd=sd[ord]))
}

#' @method AICweights HMMfits
#' @export
AICweights.HMMfits <- AICweights.miHMM

