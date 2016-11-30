
#' Confidence intervals
#'
#' Computes the standard errors and confidence intervals on the beta (i.e., working) scale of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters. Working scale depends on the real (i.e., natural) scale of the parameters:
#' 1) if both lower and upper bounds are finite then logit is the working scale;
#' 2) if lower bound is finite and upper bound is infinite then log is the working scale.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations in the computation of the CIs for the angle parameters.
#' Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{stepPar}{Standard errors and confidence intervals for the working parameters of the step lengths distribution ('gammabeta' or 'weibullbeta' depending on stepDist)}
#' \item{anglePar}{Standard errors and confidence intervals for the working parameters of the turning angles distribution ('concbeta' or 'sdbeta' depending on angleDist)}
#' \item{omegaPar}{Standard errors and confidence intervals for the working parameters of the omega distribution ('omegabeta')}
#' \item{dryPar}{Standard errors and confidence intervals for the working parameters of the dry distribution ('drybeta')}
#' \item{divePar}{Standard errors and confidence intervals for the working parameters of the dive distribution ('divebeta')}
#' \item{icePar}{Standard errors and confidence intervals for the working parameters of the ice distribution ('icebeta')}
#' \item{landPar}{Standard errors and confidence intervals for the working parameters of the land distribution ('landbeta')}
#' \item{betaPar}{Standard errors and confidence intervals for the working parameters of the transition probabilities}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CI_beta(m)
#'
#' @export
#' @importFrom MASS ginv

CI_beta <- function(m,alpha=0.95,nbSims=10^6)
{
  if(!is.momentuHMM(m))
    stop("'m' must be a momentuHMM object (as output by fitHMM)")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- length(m$stateNames)
  
  dist <- m$conditions$dist
  distnames <- names(dist)
  DM <- m$conditions$DM

  # identify covariates
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column


  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)

  p <- parDef(dist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation,m$conditions$bounds,DM)
  bounds <- p$bounds
  #if(!all(unlist(lapply(p$bounds,is.numeric)))){
  #  for(i in distnames){
  #    if(!is.numeric(bounds[[i]])){
  #      bounds[[i]] <- gsub(i,"",bounds[[i]],fixed=TRUE)
  #    }
  #  }
  #}
  
  parindex <- c(0,cumsum(unlist(lapply(DM,ncol)))[-length(DM)])
  names(parindex) <- distnames
  
  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)
  
  wpar <- m$mod$estimate
  
  Par <- list()
  for(i in distnames){
    est <- wpar[parindex[[i]]+1:ncol(DM[[i]])]^m$conditions$cons[[i]]
    var <- diag(Sigma)[parindex[[i]]+1:ncol(DM[[i]])]
    
    # if negative variance, replace by NA
    var[which(var<0)] <- NA
    
    # compute lower and upper for working parameters
    wse <- sqrt(var)
    wlower <- est-quantSup*wse
    wupper <- est+quantSup*wse
    
    if(!(dist[[i]] %in% c("wrpcauchy","vm"))){
      Par[[i]]<-parm_list(matrix(wse,ncol=length(est),byrow=T),matrix(wlower,ncol=length(est),byrow=T),matrix(wupper,ncol=length(est),byrow=T),t(bounds[[i]]))
    } else {
      if(m$conditions$estAngleMean){
        anglePar <- angleCI(m,alpha,nbSims)
      } else {
        Par[[i]]<-parm_list(matrix(wse,ncol=length(est),byrow=T),matrix(wlower,ncol=length(est),byrow=T),matrix(wupper,ncol=length(est),byrow=T),t(bounds[[i]]))
      }
    }
  }

  # group CIs for t.p. coefficients
  if(nbStates>1){
    est <- wpar[tail(cumsum(unlist(lapply(DM,ncol))),1)+1:((nbCovs+1)*nbStates*(nbStates-1))]
    var <- diag(Sigma)[tail(cumsum(unlist(lapply(DM,ncol))),1)+1:((nbCovs+1)*nbStates*(nbStates-1))]
    wse <- sqrt(var)
    wlower <- est-quantSup*wse
    wupper <- est+quantSup*wse
    Par$beta <- list(se=matrix(wse,nrow=1+nbCovs),lower=matrix(wlower,nrow=1+nbCovs),upper=matrix(wupper,nrow=1+nbCovs))
  }

  if(!is.null(m$mle$beta)) {
    rownames(Par$beta$se) <- rownames(m$mle$beta)
    rownames(Par$beta$lower) <- rownames(m$mle$beta)
    rownames(Par$beta$upper) <- rownames(m$mle$beta)
    colnames(Par$beta$se) <- colnames(m$mle$beta)
    colnames(Par$beta$lower) <- colnames(m$mle$beta)
    colnames(Par$beta$upper) <- colnames(m$mle$beta)
  }

  return(Par)
}
