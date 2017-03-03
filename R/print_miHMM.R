
#' Print \code{miHMM}
#' @method print miHMM
#'
#' @param x A \code{miHMM} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # create miHMM object from example data
#' mi <- momentuHMM:::miHMM(list(miSum=MIpool(miExample$HMMfits,ncores=1),HMMfits=miExample$HMMfits))
#' print(mi)
#'
#' @export

print.miHMM <- function(x,...)
{
  m <- x$miSum
  print(m,...)
}
