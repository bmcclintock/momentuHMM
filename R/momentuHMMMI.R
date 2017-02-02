
#' Constructor of \code{momentuHMMMI} objects
#'
#' @param m A list of attributes of the multiple imputation model: \code{Par} (the multiple imputation estimates of
#' the parameters of the model), \code{data} (the movement data), \code{conditions} (conditions used to fit
#' the model)
#'
#' @return An object \code{momentuHMM}.

momentuHMMMI <- function(m)
{
  if(is.null(m$data) | is.null(m$Par) | is.null(m$conditions))
    stop("Can't construct momentuHMM object: fields are missing")
  
  obj <- m
  
  class(obj) <- append("momentuHMMMI",class(obj))
  return(obj)
}

#' Is momentuHMMMI
#'
#' Check that an object is of class \code{\link{momentuHMMMI}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{momentuHMMMI}}, \code{FALSE} otherwise.

is.momentuHMMMI <- function(x)
  inherits(x,"momentuHMMMI")
