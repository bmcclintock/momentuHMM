
#' Constructor of \code{miSum} objects
#'
#' @param m A list of attributes required for multiple imputation summaries: \code{data} (averaged across imputations), \code{Par} (the pooled estimates of
#' the parameters of the model), \code{conditions} (conditions used to fit the model), and \code{MIcombine} (as returned by \code{\link[mitools]{MIcombine}} for
#' the working parameters).
#'
#' @return An object \code{miSum}.

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
