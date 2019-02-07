
#' Constructor of \code{momentuHierHMMData} objects
#'
#' @param data A dataframe containing: \code{ID} (the ID(s) of the observed animal(s)), \code{level} (the level of the hierarchy for each observation), and the data streams such as \code{step}
#' (the step lengths, if any), \code{angle} (the turning angles, if any), \code{x} (either easting or longitude, if any),
#' \code{y} (either norting or latitude, if any), and covariates (if any).
#'
#' @return An object \code{momentuHierHMMData}.

momentuHierHMMData <- function(data)
{
  if(is.null(data$ID) | is.null(data$level) | ncol(data)<3)
    stop("Can't construct momentuHierHMMData object: fields are missing")
  
  obj <- data
  
  class(obj) <- append(c("momentuHierHMMData","hierarchical"),class(obj))
  return(obj)
}

#' Is momentuHierHMMData
#'
#' Check that an object is of class \code{\link{momentuHierHMMData}}. Used in \code{\link{fitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{momentuHierHMMData}}, \code{FALSE} otherwise.

is.momentuHierHMMData <- function(x)
  inherits(x,"momentuHierHMMData")