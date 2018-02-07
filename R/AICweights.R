
#' Calculate Akaike information criterion model weights 
#'
#' @param ... \code{\link{momentuHMM}}, \code{\link{HMMfits}}, or \code{\link{miHMM}} objects, to compare AIC weights of the different models.
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
#' }
#' 
#' @export
#' 

AICweights <- function (..., k=2, n=NULL) {
  UseMethod("AICweights")
}

#' @method AICweights momentuHMM
#' @export
AICweights.momentuHMM <- function(..., k=2, n=NULL)
{
  models <- list(...)
  
  modNames <- all.vars(match.call()) # store the names of the models given as arguments
  
  if(any(!unlist(lapply(models,is.momentuHMM)))) stop("all model objects must be of the same class")
  
  if(length(models)<2) stop("at least 2 momentuHMM objects must be provided")
  
  for(i in 2:length(models)) {
    if(!isTRUE(all.equal(models[[i]]$data,models[[1]]$data))) stop("data must be the same for each momentuHMM object")
  }
  
  # compute AICs of models
  aic <- rep(NA,length(models))
  
  for(i in 1:length(models)) {
    aic[i] <- AIC.momentuHMM(models[[i]],k=k,n=n)
  }
  
  ord <- order(aic) # order models by increasing AIC
  weight <- exp(-0.5*(aic-min(aic)))/sum(exp(-0.5*(aic-min(aic))))
  
  return(data.frame(Model=modNames[ord],weight=weight[ord]))
  
}

#' @method AICweights miHMM
#' @export
AICweights.miHMM <- function(...,k=2, n=NULL)
{
  models <- list(...)
  
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
  for(i in 2:length(models)) {
    for(j in 1:nSims){
      if(!(j %in% momObs)){
        if(is.momentuHMM(models[[i]][[j]])){
          datNames1 <- colnames(models[[1]][[j]]$data)[colnames(models[[1]][[j]]$data) %in% colnames(models[[i]][[j]]$data)]
          datNames2 <- colnames(models[[i]][[j]]$data)[colnames(models[[i]][[j]]$data) %in% colnames(models[[1]][[j]]$data)]
          if(!isTRUE(all.equal(models[[i]][[j]]$data[,datNames2],models[[1]][[j]]$data[,datNames1]))) stop("Imputed data must be the same for each model object")
        } else {
          momObs <- c(momObs,j)
        }
      }
    }
  }
  
  if(length(momObs)) {
    momObs <- sort(momObs)
    iSims <- (1:nSims)[-momObs]
    warning("Some model fits are not momentuHMM objects. Imputations ",paste(momObs,collapse=", ")," will be ignored")
  } else iSims <- 1:nSims
  
  # compute AICs of models
  aic <- matrix(NA,length(models),nSims)
  
  for(i in 1:length(models)) {
    for(j in iSims){
      aic[i,j] <- AIC.momentuHMM(models[[i]][[j]],k=k,n=n)
    }
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
