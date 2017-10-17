
#' Print \code{miSum}
#' @method print miSum
#'
#' @param x A \code{miSum} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' \dontrun{
#' # Extract data and crawl inputs from miExample
#' obsData <- miExample$obsData
#' inits <- miExample$inits
#' err.model <- miExample$err.model
#' 
#' # Fit crawl to obsData
#' crwOut <- crawlWrap(obsData,ncores=1,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     initial.state=inits,err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,ncores=1,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,delta0=bPar$delta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' # Pool estimates
#' miSum <- MIpool(HMMfits,ncores=1)
#' print(miSum)
#' }
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
      cat(i,                "parameters:\n")      
      cat(rep("-",nchar(i)),"------------\n",sep="")
      print(m$Par$real[[i]]$est)
    } else {
      cat("Regression coeffs for",i,"parameters:\n")
      cat(rep("-",nchar(i)),"----------------------------------\n",sep="")
      print(m$Par$beta[[i]]$est)
      cat("\n")
      cat(i,                "parameters (based on mean or specified covariate values):\n")
      cat(rep("-",nchar(i)),"----------------------------------------------------------\n",sep="")
      print(x$Par$real[[i]]$est)
    }
  }
  
  if(nbStates>1){
  
    if(!is.null(m$Par$beta$beta)) {
      cat("\n")
      cat("Regression coeffs for the transition probabilities:\n")
      cat("---------------------------------------------------\n")
      print(m$Par$beta$beta$est)
    }
    
    cat("\n")
    if(!length(attr(terms.formula(m$conditions$formula),"term.labels"))) {
      cat("Transition probability matrix:\n")
      cat("------------------------------\n")
    } else {
      cat("Transition probability matrix (based on mean or specified covariate values):\n")
      cat("----------------------------------------------------------------------------\n")
    }
    print(m$Par$real$gamma$est)
    
    if(!is.null(m$Par$real$delta)){
      cat("\n")
      cat("Initial distribution:\n")
      cat("---------------------\n")
      m <- delta_bc(m)
      if(!length(attr(terms.formula(m$conditions$formulaDelta),"term.labels"))){
        tmp <- m$Par$real$delta$est[1,]
        rownames(tmp)<-NULL
        print(tmp)
      } else print(m$Par$real$delta$est)
    }
  
    if(!is.null(m$Par$timeInStates)){
      cat("\n")
      cat("Proportion of time steps spent in each state:\n")
      cat("---------------------------------------------\n")
      print(m$Par$timeInStates$est)   
    }
  }
}
