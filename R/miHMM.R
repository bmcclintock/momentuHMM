
#' Constructor of \code{miHMM} objects
#'
#' @param m A list with attributes \code{miSum} (a \code{\link{miSum}} object) and \code{HMMfits} (a list of \code{\link{momentuHMM}} objects).
#' 
#' \code{miHMM} objects are returned by \code{\link{MIfitHMM}} when arguments \code{fit=TRUE}, \code{nSims>1}, and \code{poolEstimates=TRUE}.
#'
#' @return An object \code{miHMM}.

miHMM <- function(m)
{
  if(is.null(m$miSum) | is.null(m$HMMfits))
    stop("Can't construct miHMM object: fields are missing")
  
  stopifnot(inherits(m$miSum,"miSum"))
  
  stopifnot(any(unlist(lapply(m$HMMfits,is.momentuHMM))))
  
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