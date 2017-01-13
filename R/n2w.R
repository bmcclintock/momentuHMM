
#' Scaling function: natural to working parameters.
#'
#' Scales each parameter from its natural interval to the set of real numbers, to allow for
#' unconstrained optimization. Used during the optimization of the log-likelihood.
#'
#' @param par Vector of state-dependent distributions parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in \code{par}. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param beta Matrix of regression coefficients for the transition probabilities.
#' @param delta Initial distribution. Default: \code{NULL} ; if the initial distribution is not estimated.
#' @param nbStates The number of states of the HMM.
#' @param estAngleMean \code{TRUE} if the angle mean is estimated, \code{FALSE} otherwise.
#' @param cons Logical indicating whether to use parameter constraints (cons=\code{1}) or not (cons=\code{0}).
#' 
#' @return A vector of unconstrained parameters.
#'
#' @examples
#' \dontrun{
#' nbStates <- 3
#' par <- c(0.001,0.999,0.5,0.001,1500.3,7.1)
#' bounds <- matrix(c(0,1, # bounds for first parameter
#'                    0,1, # bounds for second parameter
#'                    0,1, # ...
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf),
#'                  byrow=TRUE,ncol=2)
#' beta <- matrix(rnorm(18),ncol=6,nrow=3)
#' delta <- c(0.6,0.3,0.1)
#'
#' # vector of working parameters
#' wpar <- n2w(par=par,bounds=bounds,beta=beta,delta=delta,nbStates=nbStates,
#'            estAngleMean=FALSE)
#' }
#'
#' @importFrom boot logit

n2w <- function(par,bounds,beta,delta=NULL,nbStates,estAngleMean,DM,cons,logitcons)
{
  wpar <- NULL
  for(i in names(par)){
    if(is.null(DM[[i]])){
      if(length(which(par[[i]]<=bounds[[i]][,1] | par[[i]]>=bounds[[i]][,2]))>0)
        stop(paste0("Check the parameter bounds for ",i," (the initial parameters should be ",
                      "strictly between the bounds of their parameter space)."))
      p<-n2wDM(bounds[[i]],diag(length(par[[i]])),par[[i]],cons[[i]],logitcons[[i]])
    } else p <- par[[i]]
    wpar <- c(wpar,p)
  }

  wbeta <- as.vector(beta) # if beta is NULL, wbeta is NULL as well
  wdelta <- log(delta[-1]/delta[1]) # if delta is NULL, wdelta is NULL as well
  return(c(wpar,wbeta,wdelta))
}
