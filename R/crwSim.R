
#' Constructor of \code{crwSim} objects
#'
#' @param m A list of attributes required for multiple imputation data generated from a \code{\link{crwData}} object using \code{\link{MIfitHMM}}: \code{miData} (a list of 
#' \code{\link{momentuHMMData}} objects), and \code{crwSimulator} (a list of \code{\link[crawl]{crwSimulator}} objects).
#' 
#' \code{crwSim} objects are returned by \code{\link{MIfitHMM}} when argument \code{miData} is a \code{\link{crwData}} object and argument \code{fit=FALSE}.
#'
#' @return An object \code{crwSim}.

crwSim <- function(m)
{
  if(is.null(m$miData) | is.null(m$crwSimulator))
    stop("Can't construct crwSim object: fields are missing")
  
  stopifnot(any(unlist(lapply(m$miData,is.momentuHMMData))))
  
  stopifnot(any(unlist(lapply(m$crwSimulator,function(x) inherits(x,"crwSimulator")))))
  
  obj <- m
  
  class(obj) <- append("crwSim",class(obj))
  return(obj)
}

#' Is crwSim
#'
#' Check that an object is of class \code{\link{crwSim}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{crwSim}}, \code{FALSE} otherwise.

is.crwSim <- function(x)
  inherits(x,"crwSim")
