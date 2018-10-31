
#' Print \code{miSum}
#' @method print miSum
#'
#' @param x A \code{miSum} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' \dontrun{
#' # Extract data from miExample
#' obsData <- miExample$obsData
#' 
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#' 
#' # Fit crawl to obsData
#' crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' # Pool estimates
#' miSum <- MIpool(HMMfits)
#' print(miSum)
#' }
#' @export

print.miSum <- function(x,...)
{
  m <- x
  distnames <- names(m$conditions$dist)
  nbStates <- length(m$stateNames)
  DMind <- m$conditions$DMind
  
  if(!is.null(m$modelName)) {
    mess <- paste("Model:",m$modelName)
    cat(rep("-",nchar(mess)),"------------\n",sep="")
    cat(mess,"\n")
    cat(rep("-",nchar(mess)),"------------\n",sep="")
  }
  
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
      if(is.null(m$conditions$formulaDelta)) {
        formDelta <- ~1
      } else formDelta <- m$conditions$formulaDelta
      if(!length(attr(terms.formula(formDelta),"term.labels")) & is.null(m$conditions$formulaDelta)){
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
