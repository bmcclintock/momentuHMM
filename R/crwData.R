
#' Constructor of \code{crwData} objects
#'
#' @param m A list of attributes of crawl output: \code{crwFits} (a list of crwFit objects) and \code{crwPredict} (a crwPredict object)
#'
#'
#' @return An object \code{crwData}.
#' 
#' @seealso \code{\link{crawlWrap}}, \code{\link{MIfitHMM}}

crwData <- function(m)
{
  if(is.null(m$crwFits) | is.null(m$crwPredict))
    stop("Can't construct crwData object: fields are missing")
  
  if(any(!unlist(lapply(m$crwFits,function(x) inherits(x,"crwFit")))))
    stop("Can't construct crwData object: crwFits must be a list of crwFit objects")
  
  if(!inherits(m$crwPredict,"crwPredict")) stop("Can't construct crwData object: crwPredict must be a crwPredict object")
  
  obj <- m
  
  class(obj) <- append("crwData",class(obj))
  return(obj)
}

#' Is crwData
#'
#' Check that an object is of class \code{\link{crwData}}. Used in \code{\link{MIfitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{crwData}}, \code{FALSE} otherwise.

is.crwData <- function(x)
  inherits(x,"crwData")
