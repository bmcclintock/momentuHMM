
#' Print \code{momentuHMM}
#' @method print momentuHMM
#'
#' @param x A \code{momentuHMM} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' print(m)
#'
#' @export

print.momentuHMM <- function(x,...)
{
  m <- x
  distnames <- names(m$conditions$dist)
  DMind <- lapply(m$conditions$fullDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))

  if(length(m$mod)>1)
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")
  
  for(i in distnames){
    cat("\n")
    if(DMind[[i]]) cat(i,"parameters:\n")
    else cat("Regression coeffs for",i,"parameters:\n")
    cat("----------------------\n")
    print(m$mle[[i]])
  }
  
  if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("--------------------------------------------------\n")
    print(m$mle$beta)
  }

  if(!is.null(m$mle$gamma)) {
    cat("\n")
    cat("Transition probability matrix:\n")
    cat("-----------------------------\n")
    print(m$mle$gamma)
  }

  cat("\n")
  cat("Initial distribution:\n")
  cat("--------------------\n")
  print(m$mle$delta)
}
