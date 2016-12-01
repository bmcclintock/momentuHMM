
#' Confidence intervals for angle parameters
#'
#' Simulation-based computation of confidence intervals for the parameters of the angle distribution.
#' Used in \code{\link{CI_real}}.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations. Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{lower}{Lower bound of the confidence interval for the parameters of the angle distribution}
#' \item{upper}{Upper bound of the confidence interval for the parameters of the angle distribution}
#'
#' @importFrom MASS ginv
#' @importFrom MASS mvrnorm
#' @importFrom stats quantile

angleCI <- function(m,alpha,nbSims=10^6)
{
  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1")

  nbStates <- ncol(m$mle$anglePar)
  lower <- matrix(NA,ncol=nbStates,nrow=2)
  upper <- matrix(NA,ncol=nbStates,nrow=2)

  pdef <- parDef(m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,nbStates,estAngleMean=TRUE,
                 m$conditions$zeroInflation,m$conditions$bounds)
  parSize <- pdef$parSize
  bounds <- pdef$bounds

  for(state in 1:nbStates) {
    # working MLE
    wpar <- m$mod$estimate
    # only keep the angle parameters
    wpar <- wpar[(parSize[1]*nbStates+1):(sum(parSize)*nbStates)]
    x <- wpar[1:(length(wpar)/2)]
    y <- wpar[(length(wpar)/2+1):length(wpar)]

    # compute cov matrix
    Sigma <- ginv(m$mod$hessian)
    # only keep row/columns for angle parameters
    Sigma <- Sigma[(parSize[1]*nbStates+1):(sum(parSize)*nbStates),
                   (parSize[1]*nbStates+1):(sum(parSize)*nbStates)]
    # only keep row/columns for current state
    ind <- c(state,state+nbStates)
    Sigma <- Sigma[ind,ind]

    # simulated working parameters
    wSims <- mvrnorm(nbSims, mu=c(x[state],y[state]), Sigma=Sigma)

    # check whether some angles are close to -pi and others to pi
    theta <- Arg(wSims[,1]+1i*wSims[,2])
    c1 <- (length(which(theta>pi/2))>0) # are some angles above pi/2?
    c2 <- (length(which(theta<pi/2))>0) # are some angles below -pi/2?
    # are there more angles in [-pi,-pi/2]U[pi/2,pi] than on [-pi/2,pi/2]?
    c3 <- (length(which(theta<(-pi/2) | theta>pi/2))>length(which(theta>(-pi/2) & theta<pi/2)))
    if(c1 & c2 & c3) {
      theta_hat <- Arg(-(wSims[,1]+1i*wSims[,2])) # points are rotated to be around 0
      if(m$mle$anglePar[1,state]<0)
        theta <- theta_hat-pi # points are rotated to be around -pi
      else
        theta <- theta_hat+pi # points are rotated to be around pi
    }

    # simulated natural parameters
    nSims <- cbind(theta, sqrt(wSims[,1]^2+wSims[,2]^2))

    # scale concentration if necessary
    # (assumes that concentration has bounds like ]-Inf,b])
    if(is.finite(bounds[sum(parSize)*nbStates,2]))
      nSims[,2] <- -(exp(-nSims[,2])-bounds[sum(parSize)*nbStates,2])

    # define appropriate quantile
    quantInf <- (1-alpha)/2
    quantSup <- 1-quantInf

    # compute CIs
    lower[1,state] <- quantile(nSims[,1],quantInf)
    lower[2,state] <- quantile(nSims[,2],quantInf)
    upper[1,state] <- quantile(nSims[,1],quantSup)
    upper[2,state] <- quantile(nSims[,2],quantSup)
  }

  return(list(lower=lower,upper=upper))
}
