
#' Constructor of \code{momentuHierHMM} objects
#'
#' @param m A list of attributes of the fitted model: \code{mle} (the maximum likelihood estimates of
#' the parameters of the model), \code{data} (the \code{fitHMM} data), \code{mod} (the object
#' returned by the \code{fitHMM} numerical optimizer \code{nlm} or \code{optim}), \code{conditions} (conditions used to fit
#' the model: \code{hierStates}, \code{hierDist}, \code{zeroInflation}, \code{estAngleMean}, \code{circularAngleMean}
#' \code{stationary}, \code{formula}, \code{userBounds}, \code{bounds}, \code{workBounds}, \code{DM}, etc.), 
#' \code{stateNames}, and \code{rawCovs} (optional -- only if there are transition probability matrix covariates in the data).
#'
#' @return An object \code{momentuHierHMM}.

momentuHierHMM <- function(m)
{
  if(!is.momentuHMM(m) | is.null(m$conditions$hierStates) | is.null(m$conditions$hierDist))
    stop("Can't construct momentuHierHMM object: fields are missing")
  
  obj <- m
  
  class(obj) <- append(c("momentuHierHMM","hierarchical"),class(obj))
  return(obj)
}

#' Is momentuHierHMM
#'
#' Check that an object is of class \code{\link{momentuHierHMM}}. Used in \code{\link{CIreal}}, \code{\link{CIbeta}},
#' \code{\link{plotPR}}, \code{\link{plotStates}}, \code{\link{pseudoRes}}, \code{\link{stateProbs}},
#' and \code{\link{viterbi}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{momentuHierHMM}}, \code{FALSE} otherwise.

is.momentuHierHMM <- function(x)
  inherits(x,"momentuHierHMM")
