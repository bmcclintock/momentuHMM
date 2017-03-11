
#' Print \code{miSum}
#' @method print miSum
#'
#' @param x A \code{miSum} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # create miSum object from example data
#' mi <- MIpool(miExample$HMMfits,ncores=1)
#' print(mi)
#'
#' @export

print.miSum <- function(x,...)
{
  m <- x
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
  
  cat("\n")
  if(!length(attr(terms.formula(m$conditions$formula),"term.labels"))) {
    cat("Transition probability matrix:\n")
    cat("-----------------------------\n")
  } else {
    cat("Transition probability matrix (based on mean covariate values):\n")
    cat("--------------------------------------------------------------\n")
  }
  print(m$Par$real$gamma$est)
  
  if(!is.null(m$Par$real$delta)){
    cat("\n")
    cat("Initial distribution:\n")
    cat("--------------------\n")
    print(m$Par$real$delta$est)
  }
  
  if(nbStates>1 & !is.null(m$Par$timeInStates)){
    cat("\n")
    cat("Proportion of time steps spent in each state:\n")
    cat("--------------------------------------------\n")
    print(m$Par$timeInStates$est)   
  }
}
