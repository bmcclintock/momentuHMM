
#' Constructor of \code{momentuHMM} objects
#'
#' @param m A list of attributes of the fitted model: \code{mle} (the maximum likelihood estimates of
#' the parameters of the model), \code{data} (the \code{fitHMM} data), \code{mod} (the object
#' returned by the \code{fitHMM} numerical optimizer \code{nlm}), \code{conditions} (conditions used to fit
#' the model: \code{dist}, \code{zeroInflation}, \code{estAngleMean}, \code{circularAngleMean}
#' \code{stationary}, \code{formula}, \code{cons}, \code{userBounds}, \code{bounds}, \code{workcons}, \code{DM}, etc.), 
#' \code{stateNames}, and \code{rawCovs} (optional -- only if there are transition probability matrix covariates in the data).
#'
#' @return An object \code{momentuHMM}.

momentuHMM <- function(m)
{
  if(is.null(m$data) | is.null(m$mle) | is.null(m$mod) | is.null(m$conditions) | is.null(m$stateNames))
    stop("Can't construct momentuHMM object: fields are missing")

  obj <- m

  class(obj) <- append("momentuHMM",class(obj))
  return(obj)
}

#' Is momentuHMM
#'
#' Check that an object is of class \code{\link{momentuHMM}}. Used in \code{\link{CIreal}}, \code{\link{CIbeta}},
#' \code{\link{plotPR}}, \code{\link{plotStates}}, \code{\link{pseudoRes}}, \code{\link{stateProbs}},
#' and \code{\link{viterbi}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{momentuHMM}}, \code{FALSE} otherwise.

is.momentuHMM <- function(x)
  inherits(x,"momentuHMM")
