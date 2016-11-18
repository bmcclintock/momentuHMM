
#' Confidence intervals
#'
#' Computes the standard errors and confidence intervals on the real (i.e., natural) scale of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters.
#'
#' @param m A \code{momentuHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations in the computation of the CIs for the angle parameters.
#' Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{stepPar}{Standard errors and confidence intervals for the parameters of the step lengths distribution}
#' \item{anglePar}{Standard errors and confidence intervals for the parameters of the turning angles distribution}
#' \item{omegaPar}{Standard errors and confidence intervals for the parameters of the omega distribution}
#' \item{dryPar}{Standard errors and confidence intervals for the parameters of the dry distribution}
#' \item{divePar}{Standard errors and confidence intervals for the parameters of the dive distribution}
#' \item{icePar}{Standard errors and confidence intervals for the parameters of the ice distribution}
#' \item{landPar}{Standard errors and confidence intervals for the parameters of the land distribution}
#' \item{gamma}{Standard errors and confidence intervals for the transition probabilities. If covariates are included in TPM formula, then the covariate values for the first observation in the data are used.}
#'
#' @examples
#' # m is a momentuHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' CI_real(m)
#'
#' @export
#' @importFrom MASS ginv

CI_real <- function(m,alpha=0.95,nbSims=10^6)
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
  
  stepBounds <- bounds[stepind,]
  if(m$conditions$angleDist!="none") angleBounds <- bounds[angleind,]
  if(m$conditions$omegaDist!="none") omegaBounds <- bounds[omegaind,]
  if(m$conditions$dryDist!="none") dryBounds <- bounds[dryind,]
  if(m$conditions$diveDist!="none") diveBounds <- bounds[diveind,]
  if(m$conditions$iceDist!="none") iceBounds <- bounds[iceind,]
  if(m$conditions$landDist!="none") landBounds <- bounds[landind,]
  
  lower<-list()
  upper<-list()
  se<-list()
  
  stepPar <- get_CI(as.vector(t(m$mle$stepPar)),m,stepind,m$conditions$stepDM,stepBounds,m$conditions$cons$step,p$boundInd$step,Sigma,nbStates,alpha,m$mle$stepPar)
  
  if(m$conditions$angleDist!="none") {
    if(m$conditions$estAngleMean){
      anglePar <- angleCI(m,alpha,nbSims)
    } else {
      wangle<-m$mod$estimate[angleind]#n2wDM(angleBounds,m$conditions$angleDM,c(unique(angleDM)%*%ginv(angleDM)%*%anglePar[(1-m$conditions$estAngleMean)*nbStates+1:ncol(m$conditions$angleDM)]),m$conditions$cons$angle,logitcons=m$conditions$logitcons)
      anglelower<-angleupper<-anglese<-numeric(nrow(m$conditions$angleDM))
      for(k in 1:nrow(m$conditions$angleDM)){
        dN<-numDeriv::grad(w2nDM,wangle,bounds=angleBounds,DM=m$conditions$angleDM,cons=m$conditions$cons$angle,boundInd=p$boundInd$angle,logitcons=m$conditions$logitcons,k=k)
        anglese[k]<-suppressWarnings(sqrt(dN%*%Sigma[angleind,angleind]%*%dN))
        #cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(anglese[k]/anglePar[(1-m$conditions$estAngleMean)*nbStates+k])^2)))
        anglelower[k]<-m$mle$anglePar[-1,k]-anglese[k]*qnorm(1-(1-alpha)/2)
        angleupper[k]<-m$mle$anglePar[-1,k]+anglese[k]*qnorm(1-(1-alpha)/2)
      }
      lower$anglePar<-matrix(c(rep(NA,nbStates),anglelower),ncol=nbStates,byrow=T)
      upper$anglePar<-matrix(c(rep(NA,nbStates),angleupper),ncol=nbStates,byrow=T)  
      se$anglePar<-matrix(c(rep(NA,nbStates),anglese),ncol=nbStates,byrow=T)
      anglePar <- parm_list(se$anglePar,lower$anglePar,upper$anglePar,m$mle$anglePar)
    }
  }
  #else {
  #  low <- rbind(rep(NA,nbStates),lower$anglePar)
  #  up <- rbind(rep(NA,nbStates),upper$anglePar)
  #  anglePar <- list(lower=low,upper=up)
  #}
  
  if(m$conditions$omegaDist!="none") omegaPar <- get_CI(as.vector(t(m$mle$omegaPar)),m,omegaind,m$conditions$omegaDM,omegaBounds,m$conditions$cons$omega,p$boundInd$omega,Sigma,nbStates,alpha,m$mle$omegaPar)
  if(m$conditions$dryDist!="none") dryPar <- get_CI(as.vector(t(m$mle$dryPar)),m,dryind,m$conditions$dryDM,dryBounds,m$conditions$cons$dry,p$boundInd$dry,Sigma,nbStates,alpha,m$mle$dryPar)
  if(m$conditions$diveDist!="none") divePar <- get_CI(as.vector(t(m$mle$divePar)),m,diveind,m$conditions$diveDM,diveBounds,m$conditions$cons$dive,p$boundInd$dive,Sigma,nbStates,alpha,m$mle$divePar)
  if(m$conditions$iceDist!="none") icePar <- get_CI(as.vector(t(m$mle$icePar)),m,iceind,m$conditions$iceDM,iceBounds,m$conditions$cons$ice,p$boundInd$ice,Sigma,nbStates,alpha,m$mle$icePar)
  if(m$conditions$landDist!="none") landPar <- get_CI(as.vector(t(m$mle$landPar)),m,landind,m$conditions$landDM,landBounds,m$conditions$cons$land,p$boundInd$land,Sigma,nbStates,alpha,m$mle$landPar)
  
  #if(check)
  #  warning(paste("Some of the parameter estimates seem to lie close to the boundaries of",
  #                "their parameter space.\n  The associated CIs are probably unreliable",
  #                "(or might not be computable)."))

  if(nbStates>1) {
    # identify parameters of interest
    i2 <- nrow(bounds)+1
    i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)-1
    gammaPar <- m$mle$gamma
    wgamma <- matrix(m$mod$estimate[i2:i3],nrow=1)
    quantSup <- qnorm(1-(1-alpha)/2)
    gammalower<-gammaupper<-gammase<-matrix(NA,nbStates,nbStates)
    for(i in 1:nbStates){
      for(j in 1:nbStates){
        dN<-numDeriv::grad(get_gamma,wgamma,covs=covs,nbStates=nbStates,i=i,j=j)
        gammase[i,j]<-suppressWarnings(sqrt(dN%*%Sigma[i2:i3,i2:i3]%*%dN))
        gammalower[i,j]<-gammaPar[i,j]-quantSup*gammase[i,j]
        gammaupper[i,j]<-gammaPar[i,j]+quantSup*gammase[i,j]
      }
    }
    lower$gammaPar<-matrix(gammalower,nrow=nbStates,ncol=nbStates,byrow=T)
    upper$gammaPar<-matrix(gammaupper,nrow=nbStates,ncol=nbStates,byrow=T)  
    se$gammaPar<-matrix(gammase,nrow=nbStates,ncol=nbStates,byrow=T)
    gammaPar <- list(se=se$gamma,lower=lower$gamma,upper=upper$gamma)
    rownames(gammaPar$se) <- m$stateNames
    rownames(gammaPar$lower) <- m$stateNames
    rownames(gammaPar$upper) <- m$stateNames
    colnames(gammaPar$se) <- m$stateNames
    colnames(gammaPar$lower) <- m$stateNames
    colnames(gammaPar$upper) <- m$stateNames
  }

  # if negative variance, replace by NA
  #var[which(var<0)] <- NA

  # define appropriate quantile


  ## compute lower and upper for working parameters
  #wlower <- est-quantSup*sqrt(var)
  #wupper <- est+quantSup*sqrt(var)

  ## compute lower and upper on natural scale
  #if(m$conditions$estAngleMean) {
  #  lower <- w2n(wlower,m$conditions$stepDist,bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #  upper <- w2n(wupper,m$conditions$stepDist,bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #} else {
  #  lower <- w2n(wlower,m$conditions$stepDist,bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #  upper <- w2n(wupper,m$conditions$stepDist,bounds[1:(i2-1),],c(p$parSize[1],1),nbStates,nbCovs,FALSE,TRUE,m$cons,m$omegaDM)
  #}

  if(!is.null(m$mle$beta))
    return(list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar,gamma=gammaPar))
  else
    return(list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar))
}

get_gamma <- function(beta,covs,nbStates,i,j){
  gamma <- trMatrix_rcpp(nbStates,beta,covs)[,,1]
  gamma[i,j]
}

parm_list<-function(se,lower,upper,m){
  Par <- list(se=se,lower=lower,upper=upper)
  rownames(Par$se) <- rownames(m)
  rownames(Par$lower) <- rownames(m)
  rownames(Par$upper) <- rownames(m)
  colnames(Par$se) <- colnames(m)
  colnames(Par$lower) <- colnames(m)
  colnames(Par$upper) <- colnames(m)
  Par
}

get_CI<-function(Par,m,ind,DM,Bounds,cons,boundInd,Sigma,nbStates,alpha,mmlePar){

  w<-m$mod$estimate[ind]
  lower<-upper<-se<-numeric(nrow(DM))
  for(k in 1:nrow(DM)){
    dN<-numDeriv::grad(w2nDM,w,bounds=Bounds,DM=DM,cons=cons,boundInd=boundInd,k=k)
    se[k]<-suppressWarnings(sqrt(dN%*%Sigma[ind,ind]%*%dN))
    cn<-exp(qnorm(1-(1-alpha)/2)*sqrt(log(1+(se[k]/Par[k])^2)))
    lower[k]<-Par[k]/cn
    upper[k]<-Par[k]*cn
  }
  l<-matrix(lower,ncol=nbStates,byrow=T)
  u<-matrix(upper,ncol=nbStates,byrow=T)  
  s<-matrix(se,ncol=nbStates,byrow=T)
  out <- parm_list(s,l,u,mmlePar)
  out
}