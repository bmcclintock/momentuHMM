
#' Confidence intervals
#'
#' Computes the confidence intervals of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations in the computation of the CIs for the angle parameters.
#' Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{stepPar}{Confidence intervals for the parameters of the step lengths distribution}
#' \item{anglePar}{Confidence intervals for the parameters of the turning angles distribution}
#' \item{beta}{Confidence intervals for the regression coefficients of the transition probabilities.}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CI(m)
#'
#' @importFrom MASS ginv

CI <- function(m,alpha=0.95,nbSims=10^6)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- ncol(m$mle$stepPar)

  # identify covariates
  covsCol <- which(names(m$data)!="ID" & names(m$data)!="x" & names(m$data)!="y" &
                     names(m$data)!="step" & names(m$data)!="angle")
  nbCovs <- length(covsCol)-1 # substract intercept column

  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)
  var <- diag(Sigma)

  p <- parDef(m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,nbStates,m$conditions$estAngleMean,
              m$conditions$zeroInflation,m$bounds)

  # check if some parameters are close to their bounds
  check <- FALSE
  stepBounds <- p$bounds[1:(p$parSize[1]*nbStates),]
  stepPar <- as.vector(t(m$mle$stepPar))
  # are step parameters close to their lower bounds?
  if(length(which(round(abs(stepPar-stepBounds[,1]),5)==0))>0)
    check <- TRUE
  # are step parameters close to their upper bounds?
  if(length(which(round(abs(stepPar-stepBounds[,2]),5)==0))>0)
    check <- TRUE

  if(m$conditions$angleDist!="none") {
    angleBounds <- p$bounds[(p$parSize[1]*nbStates+1):nrow(p$bounds),]
    anglePar <- as.vector(t(m$mle$anglePar))
    # are angle parameters close to their lower bounds?
    if(length(which(round(abs(anglePar-angleBounds[,1]),5)==0))>0)
      check <- TRUE
    # are angle parameters close to their upper bounds?
    if(length(which(round(abs(anglePar-angleBounds[,2]),5)==0))>0)
      check <- TRUE
  }

  if(check)
    warning(paste("Some of the parameter estimates seem to lie close to the boundaries of",
                  "their parameter space.\n  The associated CIs are probably unreliable",
                  "(or might not be computable)."))

  # identify parameters of interest
  i1 <- p$parSize[1]*nbStates
  i2 <- sum(p$parSize)*nbStates+1
  i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)-1

  if(m$conditions$estAngleMean) {
    if(nbStates>1) {
      # select step parameters and "beta" parameters
      est <- c(m$mod$estimate[1:i1],m$mod$estimate[i2:i3])
      var <- c(var[1:i1],var[i2:i3])
    } else {
      # only select step parameters
      est <- m$mod$estimate[1:i1]
      var <- var[1:i1]
    }
  } else {
    if(nbStates>1) {
      # select step parameters, angle parameters, and "beta" parameters
      est <- c(m$mod$estimate[1:i3])
      var <- c(var[1:i3])
    } else {
      # only select step parameters and angle parameters
      est <- m$mod$estimate[1:(i2-1)]
      var <- var[1:(i2-1)]
    }
  }

  # if negative variance, replace by NA
  var[which(var<0)] <- NA

  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)

  # compute lower and upper for working parameters
  wlower <- est-quantSup*sqrt(var)
  wupper <- est+quantSup*sqrt(var)

  # compute lower and upper on natural scale
  if(m$conditions$estAngleMean) {
    lower <- w2n(wlower,m$conditions$stepDist,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
    upper <- w2n(wupper,m$conditions$stepDist,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  } else {
    lower <- w2n(wlower,m$conditions$stepDist,p$bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
    upper <- w2n(wupper,m$conditions$stepDist,p$bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  }

  # CIs for angle parameters
  if(m$conditions$estAngleMean)
    anglePar <- angleCI(m,alpha,nbSims)
  else {
    low <- rbind(rep(NA,nbStates),lower$anglePar)
    up <- rbind(rep(NA,nbStates),upper$anglePar)
    anglePar <- list(lower=low,upper=up)
  }

  # group CIs for step parameters and t.p. coefficients
  stepPar <- list(lower=lower$stepPar,upper=upper$stepPar)
  beta <- list(lower=lower$beta,upper=upper$beta)

  # name the rows and columns of the CIs
  rownames(stepPar$lower) <- rownames(m$mle$stepPar)
  rownames(stepPar$upper) <- rownames(m$mle$stepPar)
  colnames(stepPar$lower) <- colnames(m$mle$stepPar)
  colnames(stepPar$upper) <- colnames(m$mle$stepPar)
  if(m$conditions$angleDist!="none") {
    rownames(anglePar$lower) <- rownames(m$mle$anglePar)
    rownames(anglePar$upper) <- rownames(m$mle$anglePar)
    colnames(anglePar$lower) <- colnames(m$mle$anglePar)
    colnames(anglePar$upper) <- colnames(m$mle$anglePar)
  }
  if(!is.null(m$mle$beta)) {
    rownames(beta$lower) <- rownames(m$mle$beta)
    rownames(beta$upper) <- rownames(m$mle$beta)
    colnames(beta$lower) <- colnames(m$mle$beta)
    colnames(beta$upper) <- colnames(m$mle$beta)
  }

  if(!is.null(m$mle$beta))
    return(list(stepPar=stepPar,anglePar=anglePar,beta=beta))
  else
    return(list(stepPar=stepPar,anglePar=anglePar))
}
