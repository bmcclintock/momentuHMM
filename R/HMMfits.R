
#' Constructor of \code{HMMfits} objects
#'
#' @param m A list of \code{\link{momentuHMM}} objects.
#' 
#' \code{HMMfits} objects are returned by \code{\link{MIfitHMM}} when arguments \code{fit=TRUE} and \code{poolEstimates=FALSE}.
#'
#' @return An object \code{HMMfits}.

HMMfits <- function(m)
{
  stopifnot(any(unlist(lapply(m,is.momentuHMM))))
  
  obj <- m
  
  class(obj) <- append("HMMfits",class(obj))
  return(obj)
}

#' Is HMMfits
#'
#' Check that an object is of class \code{\link{HMMfits}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{HMMfits}}, \code{FALSE} otherwise.

is.HMMfits <- function(x)
  inherits(x,"HMMfits")