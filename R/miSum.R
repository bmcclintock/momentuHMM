
#' Constructor of \code{miSum} objects
#'
#' @param m A list of attributes of the multiple imputation model: \code{Par} (the multiple imputation estimates of
#' the parameters of the model), \code{data} (the movement data), \code{conditions} (conditions used to fit
#' the model)
#'
#' @return An object \code{momentuHMM}.

miSum <- function(m)
{
  if(is.null(m$data) | is.null(m$Par) | is.null(m$conditions) | is.null(m$MIcombine))
    stop("Can't construct momentuHMM object: fields are missing")
  
  obj <- m
  
  class(obj) <- append("miSum",class(obj))
  return(obj)
}

#' Is miSum
#'
#' Check that an object is of class \code{\link{miSum}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{miSum}}, \code{FALSE} otherwise.

is.miSum <- function(x)
  inherits(x,"miSum")
