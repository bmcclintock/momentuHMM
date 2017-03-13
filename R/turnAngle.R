
#' Turning angle
#'
#' Used in \code{\link{prepData}} and \code{\link{simData}}.
#'
#' @param x First point
#' @param y Second point
#' @param z Third point
#'
#' @return The angle between vectors (x,y) and (y,z)
#'
#' @examples
#' \dontrun{
#' x <- c(0,0)
#' y <- c(4,6)
#' z <- c(10,7)
#' momentuHMM:::turnAngle(x,y,z)
#' }

turnAngle <- function(x,y,z)
{
  v <- c(y[1]-x[1],y[2]-x[2])
  w <- c(z[1]-y[1],z[2]-y[2])
  angle <- atan2(w[2],w[1])-atan2(v[2],v[1])
  while(angle<=(-pi)) angle <- angle + 2*pi
  while(angle>pi) angle <- angle -2*pi
  return(angle)
}
