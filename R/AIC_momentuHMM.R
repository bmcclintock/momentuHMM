
#' AIC
#'
#' Akaike information criterion of momentuHMM model(s).
#'
#' @method AIC momentuHMM
#'
#' @param object A \code{momentuHMM} object.
#' @param ... Optional additional \code{momentuHMM} objects, to compare AICs of the different models.
#' @param k Penalty per parameter. Default: 2 ; for classical AIC.
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

AIC.momentuHMM <- function(object,...,k=2)
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

    # compute AICs of models
    AIC <- rep(NA,length(models))

    for(i in 1:length(models)) {
      m <- models[[i]]
      nbPar <- length(m$mod$estimate)-sum(!is.na(unlist(m$conditions$fixPar)))
      maxLogLike <- -m$mod$minimum
      AIC[i] <- -2*maxLogLike+k*nbPar
    }

    ord <- order(AIC) # order models by increasing AIC
    return(data.frame(Model=modNames[ord],AIC=AIC[ord]))
  }
  else { # if only one model is provided
    m <- object
    nbPar <- length(m$mod$estimate)-sum(!is.na(unlist(m$conditions$fixPar)))
    maxLogLike <- -m$mod$minimum
    AIC <- -2*maxLogLike+k*nbPar

    return(AIC)
  }
}
