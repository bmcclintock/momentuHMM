
#' Constructor of \code{crwHierSim} objects
#'
#' @param m A list of attributes required for multiple imputation data generated from a \code{\link{crwHierData}} object using \code{\link{MIfitHMM}}: \code{miData} (a list of 
#' \code{\link{momentuHMMData}} objects), and \code{crwSimulator} (a list of \code{\link[crawl]{crwSimulator}} objects).
#' 
#' \code{crwHierSim} objects are returned by \code{\link{MIfitHMM}} when argument \code{miData} is a \code{\link{crwHierData}} object and argument \code{fit=FALSE}.
#'
#' @return An object \code{crwHierSim}.

crwHierSim <- function(m)
{
  if(is.null(m$miData) | is.null(m$crwSimulator))
    stop("Can't construct crwHierSim object: fields are missing")
  
  stopifnot(any(unlist(lapply(m$miData,is.momentuHierHMMData))))
  
  stopifnot(any(unlist(lapply(m$crwSimulator,function(x) inherits(x,"crwSimulator")))))
  
  obj <- m
  
  class(obj) <- append(c("crwHierSim","hierarchical"),class(obj))
  return(obj)
}

#' Is crwHierSim
#'
#' Check that an object is of class \code{\link{crwHierSim}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{crwHierSim}}, \code{FALSE} otherwise.

is.crwHierSim <- function(x)
  inherits(x,"crwHierSim")
