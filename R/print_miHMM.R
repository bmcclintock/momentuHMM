
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
  distnames <- names(m$conditions$dist)
  nbStates <- length(m$stateNames)
  DMind <- m$conditions$DMind
  
  for(i in distnames){
    cat("\n")
    if(DMind[[i]]) {
      cat(i,"parameters:\n")      
      cat("---------------------------------------------\n")
      print(m$Par$real[[i]]$est)
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat("---------------------------------------------\n")
      print(m$Par$beta[[i]]$est)
      cat("\n")
      cat(i,"parameters (based on mean covariate values):\n")
      cat("---------------------------------------------\n")
      print(x$Par$real[[i]]$est)
    }
  }
  
  if(!is.null(m$Par$beta$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("--------------------------------------------------\n")
    print(m$Par$beta$beta$est)
  }
  
  if(!is.null(m$Par$real$gamma)) {
    cat("\n")
    cat("Transition probability matrix:\n")
    cat("-----------------------------\n")
    print(m$Par$real$gamma$est)
  }
  
  cat("\n")
  cat("Initial distribution:\n")
  cat("--------------------\n")
  print(m$Par$real$delta$est)
  
  if(nbStates>1 & !is.null(m$Par$timeInStates)){
    cat("\n")
    cat("Proportion of time steps spent in each state:\n")
    cat("--------------------------------------------\n")
    print(m$Par$timeInStates$est)   
  }
}
