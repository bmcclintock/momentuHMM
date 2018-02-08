
#' Print \code{miHMM}
#' @method print miHMM
#'
#' @param x A \code{miHMM} object.
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
#' crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     initial.state=inits,err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,delta0=bPar$delta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' miHMM <- momentuHMM:::miHMM(list(miSum=MIpool(HMMfits),HMMfits=HMMfits))
#' print(miHMM)
#' }
#' @export

print.miHMM <- function(x,...)
{
  m <- x$miSum
  print(m,...)
}
