
#' Print \code{miHMM}
#' @method print miHMM
#'
#' @param x A \code{miHMM} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # m is a miHMM object (as returned by MIfitHMM), automatically loaded with the package
#' m <- example$m
#'
#' print(m)
#'
#' @export

print.miHMM <- function(x,...)
{
  m <- x$miSum
  print(m,...)
}
