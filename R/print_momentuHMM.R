
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
  nbStates <- ncol(m$mle$stepPar)
  p <- parDef(m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,
              nbStates,TRUE,
              m$conditions$zeroInflation,NULL,m$conditions$stepDM,m$conditions$angleDM,m$conditions$omegaDM,m$conditions$dryDM,m$conditions$diveDM,m$conditions$iceDM,m$conditions$landDM)

  if(length(m$mod)>1)
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")

  cat("Step length parameters:\n")
  cat("----------------------\n")
  print(m$mle$stepPar)

  cat("\n")
  if(m$conditions$angleDist!="none") {
    cat("Turning angle parameters:\n")
    cat("------------------------\n")
    print(m$mle$anglePar)
  }

  if(m$conditions$omegaDist!="none") {
    cat("\n")
    cat("Dive proportion parameters:\n")
    cat("----------------------\n")
    print(m$mle$omegaPar)
  }
  
  if(m$conditions$dryDist!="none") {
    cat("\n")  
    cat("Dry proportion parameters:\n")
    cat("----------------------\n")
    print(m$mle$dryPar)
  }
  
  if(m$conditions$diveDist!="none") {
    cat("\n")    
    cat("Foraging dive parameters:\n")
    cat("----------------------\n")
    print(m$mle$divePar)
  }
  
  if(m$conditions$iceDist!="none") {
    cat("\n")    
    cat("Ice proportion parameters:\n")
    cat("----------------------\n")
    print(m$mle$icePar)
  }
  
  if(m$conditions$landDist!="none") {
    cat("\n")    
    cat("Land proportion parameters:\n")
    cat("----------------------\n")
    print(m$mle$landPar)
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
