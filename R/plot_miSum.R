
#' Plot \code{miSum}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot miSum
#'
#' @param x Object \code{miSum} (as return by \code{\link{MIpool}})
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in plots. If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor appearing in the data is used).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim Parameter \code{ylim} for the step length histograms.
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values.
#' @param sepAnimals If \code{TRUE}, the data is split by individuals in the histograms.
#' Default: \code{FALSE}.
#' @param sepStates If \code{TRUE}, the data is split by states in the histograms.
#' Default: \code{FALSE}.
#' @param col Vector or colors for the states (one color per state).
#' @param cumul	If TRUE, the sum of weighted densities is plotted (default).
#' @param plotTracks If TRUE, the Viterbi-decoded tracks are plotted (default).
#' @param plotCI Logical indicating whether to include confidence intervals in natural parameter plots (default: FALSE)
#' @param alpha Significance level of the confidence intervals (if \code{plotCI=TRUE}). Default: 0.95 (i.e. 95\% CIs).
#' @param plotEllipse Logical indicating whether to plot error ellipses around imputed location means. Default: TRUE.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @details The state-dependent densities are weighted by the frequency of each state in the most
#' probable state sequence (decoded with the function \code{\link{viterbi}} for each imputation). For example, if the
#' most probable state sequence indicates that one third of observations correspond to the first
#' state, and two thirds to the second state, the plots of the densities in the first state are
#' weighted by a factor 1/3, and in the second state by a factor 2/3.
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
#' # Pool estimates
#' miSum <- MIpool(HMMfits)
#' plot(miSum)
#' }
#' @export

plot.miSum <- function(x,animals=NULL,covs=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                        sepStates=FALSE,col=NULL,cumul=TRUE,plotTracks=TRUE,plotCI=FALSE,alpha=0.95,plotEllipse=TRUE,...)
{
  m <- x # the name "x" is for compatibility with the generic method
  m <- delta_bc(m)
  
  m$mle <- lapply(x$Par$real,function(x) x$est)
  m$mle$beta <- x$Par$beta$beta$est
  m$mle$delta <- x$Par$real$delta$est
  m$mod <- list()
  m$mod$estimate <- x$MIcombine$coefficients
  for(i in names(m$conditions$dist)){
    m$conditions$cons[[i]]<-rep(1,length(m$conditions$cons[[i]]))
    m$conditions$workcons[[i]]<-rep(0,length(m$conditions$workcons[[i]]))
    m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
  }
  if(!is.null(m$mle$beta)) m$conditions$workBounds$beta<-matrix(c(-Inf,Inf),length(m$mle$beta),2,byrow=TRUE)
  if(!is.null(m$Par$beta$delta$est)) m$conditions$workBoundsdelta<-matrix(c(-Inf,Inf),length(m$Par$beta$delta$est),2,byrow=TRUE)
  
  m$plotEllipse <- plotEllipse
  
  m <- momentuHMM(m)
  
  plot.momentuHMM(m,animals,covs,ask,breaks,hist.ylim,sepAnimals,
                  sepStates,col,cumul,plotTracks,plotCI,alpha,...)
}
