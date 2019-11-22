
#' Constructor of \code{crwHierData} objects
#'
#' @param m A list of attributes of crawl output: \code{crwFits} (a list of crwFit objects) and \code{crwPredict} (a crwPredict object)
#'
#'
#' @return An object \code{crwHierData}.
#' 
#' @seealso \code{\link{crawlWrap}}, \code{\link{MIfitHMM}}

crwHierData <- function(m)
{
  if(is.null(m$crwFits) | is.null(m$crwPredict))
    stop("Can't construct crwHierData object: fields are missing")
  
  if(any(!unlist(lapply(m$crwFits,function(x) inherits(x,"crwFit")))))
    stop("Can't construct crwHierData object: crwFits must be a list of crwFit objects")
  
  if(!inherits(m$crwPredict,"crwPredict")) stop("Can't construct crwHierData object: crwPredict must be a crwPredict object")
  
  if(!inherits(m$crwPredict,"hierarchical") | is.null(m$crwPredict$level)) stop("Can't construct crwHierData object: crwPredict must be of class 'hierarchical' and include a 'level' field")
  
  obj <- m
  
  class(obj) <- append(c("crwHierData","hierarchical"),class(obj))
  return(obj)
}

#' Is crwHierData
#'
#' Check that an object is of class \code{\link{crwHierData}}. Used in \code{\link{MIfitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{crwHierData}}, \code{FALSE} otherwise.

is.crwHierData <- function(x)
  inherits(x,"crwHierData")
