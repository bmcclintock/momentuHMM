
#' Constructor of \code{miHMM} objects
#'
#' @param m A list of attributes of the multiple imputation model: \code{Par} (the multiple imputation estimates of
#' the parameters of the model), \code{data} (the movement data), \code{conditions} (conditions used to fit
#' the model), \code{bounds} (parameter bounds), \code{rawCovs} (optional -- only if there are covariates
#' in the data).
#'
#' @return An object \code{miHMM}.

miHMM <- function(m)
{
  if(is.null(m$miSum) | is.null(m$HMMfits))
    stop("Can't construct miHMM object: fields are missing")
  
  stopifnot(inherits(m$miSum,"miSum"))
  
  obj <- m
  
  class(obj) <- append("miHMM",class(obj))
  return(obj)
}

#' Is miHMM
#'
#' Check that an object is of class \code{\link{miHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{miHMM}}, \code{FALSE} otherwise.

is.miHMM <- function(x)
  inherits(x,"miHMM")