
#' Plot pseudo-residuals
#'
#' Plots time series, qq-plots (against the standard normal distribution), and sample
#' ACF functions of the pseudo-residuals for each data stream
#'
#' @param m A \code{\link{momentuHMM}}, \code{\link{miHMM}}, or \code{\link{miSum}} object.
#' @param lag.max maximum lag at which to calculate the acf.  See \code{\link[stats]{acf}}.
#'
#' @details \itemize{
#' \item If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#' \item If some data streams are zero-inflated, the corresponding pseudo-
#' residuals are shown as segments, because pseudo-residuals for discrete data are defined as
#' segments (see Zucchini and MacDonald, 2009, Section 6.2).
#' \item Note that pseudo-residuals for multiple imputation analyses are based on pooled parameter 
#' estimates and the means of the data values across all imputations.
#' }
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plotPR(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats acf na.pass qqnorm

plotPR <- function(m, lag.max = NULL)
{
  if(!is.momentuHMM(m) & !is.miHMM(m) & !is.miSum(m))
    stop("'m' must be a momentuHMM, miHMM, or miSum object (as output by fitHMM, MIfitHMM, or MIpool)")
  
  if(is.miHMM(m)) m <- m$miSum
  
  distnames <- names(m$conditions$dist)
  
  if(is.miSum(m)){
    m$mle<-lapply(m$Par[distnames],function(x) x$est)
    m$mle$beta<-m$Par$beta$est
    m$mle$delta<-m$Par$delta$est
  }

  cat("Computing pseudo-residuals... ")
  pr <- pseudoRes(m)
  cat("DONE\n")

  par(mfrow=c(length(distnames),3))

  # reduce top margin
  par(mar=c(5,4,4,2)-c(0,0,3,0)) # bottom, left, top, right
  
  for(i in distnames){
    # time series
    plot(pr[[paste0(i,"Res")]],type="l",xlab="Observation index",ylab=paste0(i," pseudo-residuals"),
         main="")
    
    # qq-plot
    qq <- qqnorm(pr[[paste0(i,"Res")]],plot=FALSE)
    limInf <- min(min(qq$x,na.rm=T),min(qq$y,na.rm=T))
    limSup <- max(max(qq$x,na.rm=T),max(qq$y,na.rm=T))
    q <- qqnorm(pr[[paste0(i,"Res")]],main="",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
    
    # add segments for zero inflation
    if(m$conditions$zeroInflation[[i]]) {
      ind <- which(m$data[[i]]==0)
      x <- q$x[ind]
      y <- q$y[ind]
      segments(x,rep(limInf-5,length(ind)),x,y,col="red")
    }
    
    # add segments for one inflation
    if(m$conditions$oneInflation[[i]]) {
      ind <- which(m$data[[i]]==1)
      x <- q$x[ind]
      y <- q$y[ind]
      segments(x,rep(limInf-5,length(ind)),x,y,col="red")
    }
    
    abline(0,1,lwd=2)
    
    # ACF functions
    acf(pr[[paste0(i,"Res")]],lag.max=lag.max,na.action=na.pass,main="")
  }

  # back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
}
