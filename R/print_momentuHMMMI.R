
#' Print \code{momentuHMMMI}
#' @method print momentuHMMMI
#'
#' @param x A \code{momentuHMMMI} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @export

print.momentuHMMMI <- function(x,...)
{
  m <- x  
  distnames <- names(m$conditions$dist)
  nbStates <- length(m$stateNames)
  
  for(i in distnames){
    cat("\n")
    cat(i,"parameters:\n")
    cat("----------------------\n")
    print(m$Par[[i]])
  }
  
  if(!is.null(m$Par$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("--------------------------------------------------\n")
    print(m$Par$beta)
  }
  
  if(!is.null(m$Par$gamma)) {
    cat("\n")
    cat("Transition probability matrix:\n")
    cat("-----------------------------\n")
    print(m$Par$gamma)
  }
  
  if(nbStates>1 & !is.null(m$Par$timeInStates)){
    cat("\n")
    cat("Proportion of time steps spent in each state:\n")
    cat("--------------------------------------------\n")
    print(m$Par$timeInStates)   
  }
  
  #cat("\n")
  #cat("Initial distribution:\n")
  #cat("--------------------\n")
  #print(m$Par$delta)
}
