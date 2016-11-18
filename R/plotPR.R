
#' Plot pseudo-residuals
#'
#' Plots time series, qq-plots (against the standard normal distribution), and sample
#' ACF functions of the pseudo-residuals
#'
#' @param m A \code{\link{momentuHMM}} object
#'
#' @details \itemize{
#' \item If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#' \item If some steps are of length zero (i.e. if there is zero-inflation), the corresponding pseudo-
#' residuals are shown as segments, because pseudo-residuals for discrete data are defined as
#' segments (see Zucchini and MacDonald, 2009, Section 6.2).
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

plotPR <- function(m)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  cat("Computing pseudo-residuals... ")
  pr <- pseudoRes(m)
  cat("DONE\n")

  par(mfrow=c(3,2))

  # time series
  plot(pr$stepRes,type="l",xlab="Observation index",ylab="Steps pseudo-residuals",
       main="Steps pseudo-residuals")
  plot(pr$angleRes,type="l",xlab="Observation index",ylab="Angles pseudo-residuals",
       main="Angles pseudo-residuals")

  # reduce top margin
  par(mar=c(5,4,4,2)-c(0,0,3,0)) # bottom, left, top, right

  # steps qq-plot
  qqStep <- qqnorm(pr$stepRes,plot=FALSE)
  limInf <- min(min(qqStep$x,na.rm=T),min(qqStep$y,na.rm=T))
  limSup <- max(max(qqStep$x,na.rm=T),max(qqStep$y,na.rm=T))
  q <- qqnorm(pr$stepRes,main="",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))

  # add segments for steps of length zero
  if(m$conditions$zeroInflation) {
    ind <- which(m$data$step==0)
    x <- q$x[ind]
    y <- q$y[ind]
    segments(x,rep(limInf-5,length(ind)),x,y,col="red")
  }

  abline(0,1,lwd=2)

  # angles qq-plot
  qqAngle <- qqnorm(pr$angleRes,plot=FALSE)
  limInf <- min(min(qqAngle$x,na.rm=T),min(qqAngle$y,na.rm=T))
  limSup <- max(max(qqAngle$x,na.rm=T),max(qqAngle$y,na.rm=T))
  qqnorm(pr$angleRes,main="",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
  abline(0,1,lwd=2)

  # ACF functions
  acf(pr$stepRes,na.action=na.pass,main="")
  acf(pr$angleRes,na.action=na.pass,main="")

  # back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
}
