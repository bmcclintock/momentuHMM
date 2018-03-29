
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
  DMind <- m$conditions$DMind
  
  if(!is.null(m$modelName)) {
    mess <- paste("Model:",m$modelName)
    cat(rep("-",nchar(mess)),"------------\n",sep="")
    cat(mess,"\n")
    cat(rep("-",nchar(mess)),"------------\n\n",sep="")
  }
  
  if(length(m$mod)>1)
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")
  
  for(i in distnames){
    cat("\n")
    if(DMind[[i]]) {
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      print(m$mle[[i]])
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      print(m$CIbeta[[i]]$est)
    }

    if(!DMind[[i]]){
      cat("\n")
      cat(i,                "parameters (based on mean covariate values):\n")
      cat(rep("-",nchar(i)),"---------------------------------------------\n",sep="")
      print(x$CIreal[[i]]$est)
    }
  }
  
  if(length(m$stateNames)>1){
    #if(!is.null(m$mle$beta)) {
      cat("\n")
      cat("Regression coeffs for the transition probabilities:\n")
      cat("---------------------------------------------------\n")
      print(m$mle$beta)
    #}
  
    if(!is.null(m$mle$gamma)) {
      cat("\n")
      cat("Transition probability matrix:\n")
      cat("------------------------------\n")
      print(m$mle$gamma)
    } else {
      cat("\n")
      cat("Transition probability matrix (based on mean covariate values):\n")
      cat("---------------------------------------------------------------\n")
      print(m$CIreal$gamma$est)    
    }
  
    cat("\n")
    cat("Initial distribution:\n")
    cat("---------------------\n")
    m <- delta_bc(m)
    if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
      tmp <- m$mle$delta[1,]
      rownames(tmp)<-NULL
      print(tmp)
    } else print(m$mle$delta)
  }
}
