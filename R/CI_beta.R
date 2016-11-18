
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

  nbStates <- ncol(m$mle$stepPar)

  # identify covariates
  covsCol <- which(names(m$data)!="ID" & names(m$data)!="x" & names(m$data)!="y" &
                     names(m$data)!="step" & names(m$data)!="angle" & names(m$data)!="omega" & names(m$data)!="dry" & names(m$data)!="dive" & names(m$data)!="ice" & names(m$data)!="land")
  covs <- model.matrix(m$conditions$formula,m$data)
  nbCovs <- ncol(covs)-1 # substract intercept column

  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)
  var <- diag(Sigma)

  p <- parDef(m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,nbStates,m$conditions$estAngleMean,
              m$conditions$zeroInflation,m$bounds,m$conditions$stepDM,m$conditions$angleDM,m$conditions$omegaDM,m$conditions$dryDM,m$conditions$diveDM,m$conditions$iceDM,m$conditions$landDM)

  stepPar <- as.vector(t(m$mle$stepPar))
  anglePar<-omegaPar<-dryPar<-divePar<-icePar<-landPar<-NULL
  if(m$conditions$angleDist!="none") anglePar <- as.vector(t(m$mle$anglePar))
  if(m$conditions$omegaDist!="none") omegaPar <- as.vector(t(m$mle$omegaPar))
  if(m$conditions$dryDist!="none") dryPar <- as.vector(t(m$mle$dryPar))
  if(m$conditions$diveDist!="none") divePar <- as.vector(t(m$mle$divePar))
  if(m$conditions$iceDist!="none") icePar <- as.vector(t(m$mle$icePar))
  if(m$conditions$landDist!="none") landPar <- as.vector(t(m$mle$landPar))
  
  bounds <- p$bounds
  if(!is.numeric(bounds)){
    bounds<-matrix(sapply(bounds,function(x) eval(parse(text=x))),ncol=2,dimnames=list(rownames(p$bounds)))
  }
  
  stepind <- which(grepl(m$conditions$stepDist,rownames(bounds)))
  angleind <- which(grepl(ifelse(m$conditions$angleDist=="wrpcauchy","conc","sd"),rownames(bounds)))
  omegaind <- which(grepl("omega",rownames(bounds)))
  dryind <- which(grepl("dry",rownames(bounds)))
  diveind <- which(grepl("dive",rownames(bounds)))
  iceind <- which(grepl("ice",rownames(bounds)))
  landind <- which(grepl("land",rownames(bounds)))

  # identify parameters of interest
  i1 <- length(stepind)
  i2 <- nrow(bounds)+1
  i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)-1
  
  estcons <- m$mod$estimate
  estcons[1:(i2-1)] <- estcons[1:(i2-1)]^(unlist(m$conditions$cons)[c(stepind,angleind,omegaind,dryind,diveind,iceind,landind)])

  if(m$conditions$estAngleMean) {
    if(nbStates>1) {
      # select step parameters and "beta" parameters
      est <- c(estcons[1:i1],estcons[(i1+length(angleind)+1):i3])
      var <- c(var[1:i1],var[(i1+length(angleind)+1):i3])
    } else {
      # only select step parameters
      est <- estcons[1:i1,(i1+length(angleind)+1):(i2-1)]
      var <- var[1:i1,(i1+length(angleind)+1):(i2-1)]
    }
  } else {
    if(nbStates>1) {
      # select step parameters, angle parameters, and "beta" parameters
      est <- c(estcons[1:i3])
      var <- c(var[1:i3])
    } else {
      # only select step parameters and angle parameters
      est <- estcons[1:(i2-1)]
      var <- var[1:(i2-1)]
    }
  }

  # if negative variance, replace by NA
  var[which(var<0)] <- NA

  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)

  # compute lower and upper for working parameters
  wse <- sqrt(var)
  wlower <- est-quantSup*wse
  wupper <- est+quantSup*wse
  
  stepPar <- parm_list(matrix(wse[stepind],ncol=length(stepind),byrow=T),matrix(wlower[stepind],ncol=length(stepind),byrow=T),matrix(wupper[stepind],ncol=length(stepind),byrow=T),t(bounds[stepind,]))
  
  anglePar <- omegaPar <- dryPar <- divePar <-  icePar <- landPar <-  NULL
  
  if(m$conditions$angleDist!="none") {
    if(m$conditions$estAngleMean){
      anglePar <- angleCI(m,alpha,nbSims)
    } else {
      se <- matrix(wse[angleind],ncol=length(angleind),byrow=T)
      lower <- matrix(wlower[angleind],ncol=length(angleind),byrow=T)
      upper <- matrix(wupper[angleind],ncol=length(angleind),byrow=T)    
      anglePar <- parm_list(se,lower,upper,t(bounds[angleind,]))
    }
  }
  
  if(m$conditions$omegaDist!="none") omegaPar <- parm_list(matrix(wse[omegaind],ncol=length(omegaind),byrow=T),matrix(wlower[omegaind],ncol=length(omegaind),byrow=T),matrix(wupper[omegaind],ncol=length(omegaind),byrow=T),t(bounds[omegaind,]))
  if(m$conditions$dryDist!="none") dryPar <- parm_list(matrix(wse[dryind],ncol=length(dryind),byrow=T),matrix(wlower[dryind],ncol=length(dryind),byrow=T),matrix(wupper[dryind],ncol=length(dryind),byrow=T),t(bounds[dryind,]))
  if(m$conditions$diveDist!="none") divePar <- parm_list(matrix(wse[diveind],ncol=length(diveind),byrow=T),matrix(wlower[diveind],ncol=length(diveind),byrow=T),matrix(wupper[diveind],ncol=length(diveind),byrow=T),t(bounds[diveind,]))
  if(m$conditions$iceDist!="none") icePar <- parm_list(matrix(wse[iceind],ncol=length(iceind),byrow=T),matrix(wlower[iceind],ncol=length(iceind),byrow=T),matrix(wupper[iceind],ncol=length(iceind),byrow=T),t(bounds[iceind,]))
  if(m$conditions$landDist!="none") landPar <- parm_list(matrix(wse[landind],ncol=length(landind),byrow=T),matrix(wlower[landind],ncol=length(landind),byrow=T),matrix(wupper[landind],ncol=length(landind),byrow=T),t(bounds[landind,]))
  
  # compute lower and upper on natural scale
  #if(m$conditions$estAngleMean) {
  #  lower <- w2n(wlower,m$conditions$stepDist,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #  upper <- w2n(wupper,m$conditions$stepDist,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #} else {
  #  lower <- w2n(wlower,m$conditions$stepDist,p$bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #  upper <- w2n(wupper,m$conditions$stepDist,p$bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #}

  # group CIs for t.p. coefficients
  beta <- list(se=matrix(wse[i2:i3],nrow=1+nbCovs),lower=matrix(wlower[i2:i3],nrow=1+nbCovs),upper=matrix(wupper[i2:i3],nrow=1+nbCovs))


  if(!is.null(m$mle$beta)) {
    rownames(beta$se) <- rownames(m$mle$beta)
    rownames(beta$lower) <- rownames(m$mle$beta)
    rownames(beta$upper) <- rownames(m$mle$beta)
    colnames(beta$se) <- colnames(m$mle$beta)
    colnames(beta$lower) <- colnames(m$mle$beta)
    colnames(beta$upper) <- colnames(m$mle$beta)
  }

  if(!is.null(m$mle$beta))
    return(list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar,betaPar=beta))
  else
    return(list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar))
}
