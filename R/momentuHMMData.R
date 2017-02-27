
#' Constructor of \code{momentuHMMData} objects
#'
#' @param data A dataframe containing: \code{ID} (the ID(s) of the observed animal(s)) and the data streams such as \code{step}
#' (the step lengths, if any), \code{angle} (the turning angles, if any), \code{x} (either easting or longitude, if any),
#' \code{y} (either norting or latitude, if any), and covariates (if any).
#'
#' @return An object \code{momentuHMMData}.

momentuHMMData <- function(data)
{
  if(is.null(data$ID) | ncol(data)<2)
    stop("Can't construct momentuHMMData object: fields are missing")

  obj <- data

  class(obj) <- append("momentuHMMData",class(obj))
  return(obj)
}

#' Is momentuHMMData
#'
#' Check that an object is of class \code{\link{momentuHMMData}}. Used in \code{\link{fitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{momentuHMMData}}, \code{FALSE} otherwise.

is.momentuHMMData <- function(x)
  inherits(x,"momentuHMMData")
