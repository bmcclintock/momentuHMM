
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
  nbStates <- ncol(m$Par$stepPar$est)
  
  cat("Step length parameters:\n")
  cat("----------------------\n")
  print(m$Par$stepPar)
  
  cat("\n")
  if(m$conditions$angleDist!="none") {
    cat("Turning angle parameters:\n")
    cat("------------------------\n")
    print(m$Par$anglePar)
  }
  
  cat("\n")
  if(m$conditions$omegaDist!="none") {
    cat("Dive proportion parameters:\n")
    cat("----------------------\n")
    print(m$Par$omegaPar)
  }
  
  cat("\n")  
  if(m$conditions$dryDist!="none") {
    cat("Dry proportion parameters:\n")
    cat("----------------------\n")
    print(m$Par$dryPar)
  }
  
  cat("\n")    
  if(m$conditions$diveDist!="none") {
    cat("Foraging dive parameters:\n")
    cat("----------------------\n")
    print(m$Par$divePar)
  }
  
  cat("\n") 
  if(m$conditions$iceDist!="none") {
    cat("Ice proportion parameters:\n")
    cat("----------------------\n")
    print(m$Par$icePar)
  }
  
  cat("\n")  
  if(m$conditions$landDist!="none") {
    cat("Land proportion parameters:\n")
    cat("----------------------\n")
    print(m$Par$landPar)
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
