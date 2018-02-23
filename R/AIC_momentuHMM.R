
#' AIC
#'
#' Akaike information criterion of momentuHMM model(s).
#'
#' @method AIC momentuHMM
#'
#' @param object A \code{momentuHMM} object.
#' @param ... Optional additional \code{momentuHMM} objects, to compare AICs of the different models.
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
#' @export

AIC.momentuHMM <- function(object,...,k=2,n=NULL)
{
  models <- list(...)

  if(length(models)>0) { # if several models are provided
    modNames <- all.vars(match.call()) # store the names of the models given as arguments

    # include "object" in "models"
    modcopy <- list()
    modcopy[[1]] <- object
    for(i in 1:length(models))
      modcopy[[i+1]] <- models[[i]]
    models <- modcopy
    
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
      aic[i] <- getAIC(models[[i]],k,n)
    }

    ord <- order(aic) # order models by increasing AIC
    return(data.frame(Model=modNames[ord],AIC=aic[ord]))
  }
  else { # if only one model is provided
    aic <- getAIC(object,k,n)
    return(aic)
  }
}

getAIC <- function(m, k=2, n=NULL){
    nbPar <- length(m$mod$estimate)-sum(!is.na(unlist(m$conditions$fixPar)))
    maxLogLike <- -m$mod$minimum
    aic <- -2*maxLogLike+k*nbPar
    if(!is.null(n)) aic <- aic + k*nbPar*(nbPar+1)/(n-nbPar-1)
    aic
}