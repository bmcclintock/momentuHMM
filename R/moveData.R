
#' Constructor of \code{moveData} objects
#'
#' @param data A dataframe containing: \code{ID} (the ID(s) of the observed animal(s)), \code{step}
#' (the step lengths), \code{angle} (the turning angles, if any), \code{x} (either easting or longitude),
#' \code{y} (either norting or latitude), and covariates, if any.
#'
#' @return An object \code{moveData}.

moveData <- function(data)
{
  if(is.null(data$ID) | is.null(data$step) | is.null(data$x) | is.null(data$y))
    stop("Can't construct moveData object: fields are missing")

  obj <- data

  class(obj) <- append("moveData",class(obj))
  return(obj)
}

#' Is moveData
#'
#' Check that an object is of class \code{\link{moveData}}. Used in \code{\link{fitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{moveData}}, \code{FALSE} otherwise.

is.moveData <- function(x)
  inherits(x,"moveData")
