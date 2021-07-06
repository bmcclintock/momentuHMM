
#' Random effects estimation
#'
#' Approximate individual-level random effects estimation for state transition probabilities based on Burnham & White (2002)
#'
#' @param m A \code{\link{momentuHMM}} object.
#' @param Xformula Formula for the design matrix of the random effects model. The default \code{Xformula=~1} specifies an intercept-only model with no additional individual covariate effects.
#' @param alpha Significance level of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param ncores number of cores to use for parallel processing
#' @param nlmPar List of parameters to pass to the optimization function \code{\link[stats]{nlm}}. See \code{\link{fitHMM}}.
#' @param fit \code{TRUE} if the HMM should be re-fitted at the shrinkage estimates, \code{FALSE} otherwise. 
#' @param retryFits Non-negative integer indicating the number of times to attempt to iteratively fit the model using random perturbations of the current parameter estimates as the 
#' initial values for likelihood optimization. See \code{\link{fitHMM}}.
#' @param retrySD An optional list of scalars or vectors indicating the standard deviation to use for normal perturbations of each working scale parameter when \code{retryFits>0}. See \code{\link{fitHMM}}.
#' @param optMethod The optimization method to be used. See \code{\link{fitHMM}}.
#' @param control A list of control parameters to be passed to \code{\link[stats]{optim}} (ignored unless \code{optMethod="Nelder-Mead"} or \code{optMethod="SANN"}).
#' @param modelName An optional character string providing a name for the fitted model. See \code{\link{fitHMM}}.
#' @param ... further arguments passed to or from other methods. Not currently used.
#' 
#' @return A \code{randomEffects} model similar to a \code{\link{momentuHMM}} object, but including the additional random effect components:
#' \item{varcomp}{A list of length \code{nbStates*(nbStates-1)} with each element containing the random effect mean coefficient(s) (\code{mu}), random effect variance (\code{sigma}), 
#' and logit-scale shrinkage estimates for the state transition probability parameters (\code{ztilde}).}
#' \item{traceG}{The trace of the projection matrix for each random effect.}
#'
#' @examples
#' \dontrun{
#' # simulated data with normal random effects
#' # and binary individual covariate 
#' 
#' nbAnimals <- 5 # should be larger for random effects estimation
#' obsPerAnimal <- 110
#' indCov <- rbinom(nbAnimals,1,0.5) # individual covariate
#' betaCov <- c(-0.5,0.5) # covariate effects
#' mu <- c(-0.1,0.1) # mean for random effects
#' sigma <- c(0.2,0.4) # sigma for random effects
#' beta0 <- cbind(rnorm(nbAnimals,mu[1],sigma[1]),
#'                rnorm(nbAnimals,mu[2],sigma[2]))
#' 
#' reData <- simData(nbAnimals=nbAnimals,obsPerAnimal=obsPerAnimal,nbStates=2,
#'                   dist=list(step="gamma"),formula=~0+ID+indCov,
#'                   Par=list(step=c(1,10,1,2)),
#'                   beta=rbind(beta0,betaCov),
#'                   covs=data.frame(indCov=rep(indCov,each=obsPerAnimal)))
#'
#' # fit null model
#' nullFit <- fitHMM(reData,nbStates=2,
#'                   dist=list(step="gamma"),
#'                   Par0=list(step=c(1,10,1,2)))
#' 
#' # fit covariate model
#' covFit <- fitHMM(reData,nbStates=2,
#'                  dist=list(step="gamma"),formula=~indCov,
#'                  Par0=list(step=c(1,10,1,2)),
#'                  beta0=rbind(mu,betaCov)) 
#'
#' # fit fixed effects model
#' fixFit <- fitHMM(reData,nbStates=2,
#'                  dist=list(step="gamma"),formula=~0+ID,
#'                  Par0=list(step=c(1,10,1,2)),
#'                  beta0=beta0)
#' 
#' # fit random effect model
#' reFit <- randomEffects(fixFit)
#' 
#' # fit random effect model with individual covariate
#' reCovFit <- randomEffects(fixFit, Xformula=~indCov)
#' 
#' # compare by AICc
#' AIC(nullFit,covFit,fixFit,reFit,reCovFit, n=nrow(reData))
#' }
#' @references
#' Burnham, K.P. and White, G.C. 2002. Evaluation of some random effects methodology applicable to bird ringing data. Journal of Applied Statistics 29: 245-264.
#' 
#' McClintock, B.T. 2021. Worth the effort? A practical examination of random effects in hidden Markov models for animal telemetry data. Methods in Ecology and Evolution \doi{10.1111/2041-210X.13619}.
#' @export
#' @importFrom MASS ginv
# @importFrom BB BBsolve 
#' @importFrom stats optimise get_all_vars
# @importFrom expm sqrtm
# @importFrom matrixcalc matrix.trace

randomEffects <- function(m, Xformula = ~1, alpha = 0.95, ncores = 1, nlmPar = list(), fit = TRUE, retryFits = 0, retrySD = NULL, optMethod = "nlm", control = list(), modelName = NULL, ...)
{
  
  for(pkg in c("BB","expm","matrixcalc")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package \"",pkg,"\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }
  
  modelterms <- colnames(stats::get_all_vars(m$conditions$formula,m$data))
  if(!("ID" %in% modelterms)) stop("'ID' must be included in the state transition probability formula of the model")
  if(!all(modelterms=="ID")) stop("Only 'ID' is permitted in the state transition probability formula of the model")

  mName <- deparse(substitute(m))
  m <- delta_bc(m)
  
  chkDots(...)
  
  nbStates <- length(m$stateNames)
  nbRE <- nbStates*(nbStates-1)
  
  VC <- getVC(m)
  
  nbAnimals <- length(unique(m$data$ID))
  if(nbAnimals<2) stop("number of individuals must be >1")
  
  # aInd = list of indices of first observation for each animal
  aInd <- NULL
  for(i in 1:nbAnimals){
    idInd <- which(m$data$ID==unique(m$data$ID)[i])
    aInd <- c(aInd,idInd[1])
  }
  
  if(!is.formula(Xformula)) stop("Xformula must be a formula")
  X <- stats::model.matrix(Xformula,m$data)
  if(any(mapply(function(x) nrow(unique(X[which(m$data$ID==x),,drop=FALSE]))>1,unique(m$data$ID))))
    stop("time-varying covariates are not permitted in Xformula")
  X <- X[aInd,,drop=FALSE]
  
  trProbs <- getTrProbs(m,covIndex=aInd)
  if(any(trProbs<0.01) | any(trProbs>0.99)) warning("estimated state transition probabilites for ",mName," appear to be near a boundary; proceed with caution")
  
  if(ncores>1 & nbRE>1){
    for(pkg in c("doFuture","future")){
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"",pkg,"\" needed for parallel processing to work. Please install it.",
             call. = FALSE)
      }
    }
    oldDoPar <- doFuture::registerDoFuture()
    on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    # hack so that foreach %dorng% can find internal momentuHMM variables without using ::: (forbidden by CRAN)
    getZtilde <- getZtilde
    progBar <- progBar
    pkgs <- c("momentuHMM")
  } else { 
    doParallel::registerDoParallel(cores=ncores)
    pkgs <- NULL
  }
  varcomp <- withCallingHandlers(foreach(j=1:nbRE,.packages = pkgs) %dorng% {
    cat("Random effect for state transition ",colnames(m$mle$beta)[j]," ...\n",sep="")
    getZtilde(W=VC[(j-1)*nbAnimals+1:nbAnimals,(j-1)*nbAnimals+1:nbAnimals],betahat=c(m$mle$beta[,j]),X=X,alpha)
  }
  ,warning=muffleRNGwarning)
  if(!((ncores>1 & nbRE>1))) doParallel::stopImplicitCluster()
  else future::plan(future::sequential)
  cat("DONE\n")
  
  betaFix <- matrix(NA,nrow(m$mle$beta),ncol(m$mle$beta))
  for(j in 1:nbRE){
    betaFix[,j] <- varcomp[[j]]$ztilde$est
  }
  fixPar <- m$conditions$fixPar[names(m$conditions$fixPar)[which(unlist(lapply(m$condition$fixPar,function(x) any(is.numeric(x)))))]]
  fixPar$beta <- betaFix
  
  Par0 <- getPar(m)
  
  if(fit){
    what <- "Re-fitting"
    dist <- m$conditions$dist
    distnames <- names(dist)
    message("=======================================================================")
    if(!inherits(m$data,"hierarchical")) message(what," HMM with ",nbStates," state",ifelse(nbStates>1,"s","")," and ",length(distnames)," data stream",ifelse(length(distnames)>1,"s",""))
    else message(what," hierarchical HMM with ",nbStates," state",ifelse(nbStates>1,"s","")," and ",length(distnames)," data stream",ifelse(length(distnames)>1,"s",""))
    message("-----------------------------------------------------------------------\n")
    for(i in distnames){
      pNames<-rownames(m$mle[[i]])
      if(is.null(m$conditions$DM[[i]])){
        message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=~1",collapse=", "),")")
      } else if(is.list(m$conditions$DM[[i]])){
        message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,"=",m$conditions$DM[[i]],collapse=", "),")")
      } else message(" ",i," ~ ",dist[[i]],"(",paste0(pNames,": custom",collapse=", "),")")
    }
    message("\n State transition probability random effect formula: ",paste0(Xformula,collapse=""))
    if(is.null(m$conditions$formulaDelta)){
      formDelta <- ~1
    } else formDelta <- m$conditions$formulaDelta
    message("\n Initial distribution formula: ",ifelse(m$conditions$stationary,"stationary",paste0(formDelta,collapse="")))
    if(m$conditions$mixtures>1) {
      message("\n Number of mixtures: ",m$conditions$mixtures)
      message(" Mixture probability formula: ",paste0(m$conditions$formulaPi,collapse=""))
    }
    message("=======================================================================")
  }
  mod <- tryCatch(
           suppressMessages(fitHMM(m$data,nbStates=nbStates,dist=m$conditions$dist,
           Par0=Par0$Par,beta0=betaFix,delta0=Par0$delta,
           estAngleMean=m$conditions$estAngleMean,circularAngleMean=m$conditions$circularAngleMean,
           formula=m$conditions$formula,formulaDelta=m$conditions$formulaDelta,stationary=m$conditions$stationary,
           nlmPar=nlmPar,fit=fit,
           DM=m$conditions$DM,userBounds=m$conditions$userBounds,workBounds=m$conditions$workBounds,
           betaCons=m$conditions$betaCons,betaRef=m$conditions$betaRef,deltaCons=m$conditions$deltaCons,
           mvnCoords=m$conditions$mvnCoords,stateNames=m$stateNames,knownStates=m$knownStates,
           fixPar=fixPar,retryFits=retryFits,retrySD=retrySD,optMethod=optMethod,control=control,
           prior=m$prior,modelName=modelName,...))
           ,error=function(e) e)
  
  if(!inherits(mod,"error") & fit){
    for(j in 1:length(varcomp)){
      diag(mod$mod$Sigma)[4+(j-1)*nbAnimals+1:nbAnimals] <- varcomp[[j]]$ztilde$se^2
    }
    mod$CIbeta <- CIbeta(mod)
    for(i in names(mod$CIbeta$beta)){
      for(j in 1:length(varcomp)){
        mod$CIbeta$beta[[i]][,j] <- varcomp[[j]]$ztilde[[i]]
      }
    }
  }
  
  names(varcomp) <- colnames(m$mle$beta)
  
  varc <- list()
  for(j in names(varcomp)){
    if(inherits(varcomp[[j]]$sigma,"badSigma")) warning(paste0("sigma for ",j," does not appear to have converged"))
    varc[[j]] <- list(mu=varcomp[[j]]$muhat,
                    sigma=varcomp[[j]]$sigma,
                    ztilde=varcomp[[j]]$ztilde)
  }
  mod$varcomp <- varc
  mod$traceG <- lapply(varcomp,function(x) x$traceG)
  
  class(mod) <- append(class(mod),"randomEffects")
  if(fit) message(ifelse(retryFits>=1,"\n",""),"DONE")
  
  return(mod)
}

getZtilde <- function(W,betahat,X,alpha){
  
  quantSup <- qnorm(1-(1-alpha)/2)
  
  sigma <- ztilde <- list()
  df <- c(length(betahat)-ncol(X),qchisq(1-(1-alpha)/2,length(betahat)-ncol(X)),qchisq((1-alpha)/2,length(betahat)-ncol(X)))
  names(df) <- c("est","lower","upper")
  
  for(j in 1:length(df)){
    tmp <- list()
    tmp[[1]] <- stats::optimise(opt,c(-10,100),W=W,betahat=betahat,X=X,df=df[j],tol=sqrt(.Machine$double.eps))
    tmp[[2]] <- stats::optimise(opt,c(-100,100),W=W,betahat=betahat,X=X,df=df[j],tol=sqrt(.Machine$double.eps))
    tmp[[3]] <- stats::optimise(opt,c(-1000,1000),W=W,betahat=betahat,X=X,df=df[j],tol=sqrt(.Machine$double.eps))
    tmp[[4]] <- stats::optimise(opt,c(-1,10),W=W,betahat=betahat,X=X,df=df[j],tol=sqrt(.Machine$double.eps))
    sig <- tmp[[which.min(unlist(lapply(tmp,function(x) opt12(x$minimum,W=W,betahat=betahat,X=X,df=df[j]))))]]$minimum
    sigma[[names(df)[j]]] <- BB::BBsolve(sig,opt12,control=list(tol=1.e-12,trace=1),W=W,betahat=betahat,X=X,df=df[j],opt="solve")$par[1]
  }
  
  if(sigma$est < sigma$lower | sigma$est > sigma$upper){
    class(sigma) <- append(class(sigma),"badSigma")
  }
  if(sigma$est < 0) sigma$est <- 0
  if(sigma$lower < 0) sigma$lower <- 0
  if(sigma$upper < 0) sigma$upper <- 0
  sigma$est <- sqrt(sigma$est)
  sigma$lower <- min(sigma$est,sqrt(sigma$lower))
  sigma$upper <- max(sigma$est,sqrt(sigma$upper))
  D <- sigma$est^2 * diag(length(betahat)) + W
  Dinv <- MASS::ginv(D)
  if(sigma$est==0) H <- matrix(0,nrow(W),ncol(W)) #sigma$est * MASS::ginv(expm::sqrtm(D))
  else H <- MASS::ginv(expm::sqrtm(diag(length(betahat))+1/(sigma$est^2)*W))
  if(any(is.complex(c(H)))){
    if(isTRUE(all.equal(Im(H),matrix(0,nrow(H),ncol(H))))) H <- Re(H)
  }
  
  muhat <- list()
  muhatVC <- MASS::ginv(t(X)%*%Dinv%*%X)
  muhat$est <- c(muhatVC%*%t(X)%*%Dinv%*%betahat)
  muhat$se <- c(sqrt(diag(muhatVC)))
  muhat$lower <- c(muhat$est - quantSup*muhat$se)
  muhat$upper <- c(muhat$est + quantSup*muhat$se)
  names(muhat$est) <- names(muhat$se) <- names(muhat$lower) <- names(muhat$upper) <- colnames(X)
  #muhat$VC <- muhatVC
  
  A <- X %*% MASS::ginv(t(X)%*%Dinv%*%X) %*% t(X)
  G <- H + (diag(length(betahat))-H) %*% A %*% Dinv
  ztilde$est <- c(H %*% (betahat-X%*%muhat$est) + X %*% muhat$est )
  #ztilde <- G %*% betahat
  
  ztildeVC <- G %*% W %*% t(G)
  for(k in 1:length(betahat)){
    for(kk in 1:length(betahat)){
      ztildeVC[k,kk] <- ztildeVC[k,kk] + (ztilde$est[k]-betahat[k])*(ztilde$est[kk]-betahat[kk])
    }
  }
  ztilde$se <- c(sqrt(diag(ztildeVC))) #rmse
  ztilde$lower <- c(ztilde$est - quantSup*ztilde$se)
  ztilde$upper <- c(ztilde$est + quantSup*ztilde$se)
  names(ztilde$est) <- names(ztilde$se) <- names(ztilde$lower) <- names(ztilde$upper) <- rownames(ztildeVC) <- colnames(ztildeVC) <- paste0("ID",1:nrow(X))
  #ztilde$VC <- ztildeVC
  
  traceG <- matrixcalc::matrix.trace(G)
  if(is.complex(traceG)){
    if(isTRUE(all.equal(Im(traceG),0))) traceG <- Re(traceG)
  }
  return(list(sigma=sigma, ztilde=ztilde, muhat=muhat, traceG = traceG))
}

getVC <- function(m){
  
  nbStates <- length(m$stateNames)
  mixtures <- m$conditions$mixtures
  reForm <- formatRecharge(nbStates,m$conditions$formula,m$conditions$betaRef,m$data,par=m$mle)
  recharge <- reForm$recharge
  hierRecharge <- reForm$hierRecharge
  newformula <- reForm$newformula
  covs <- reForm$covs
  nbCovs <- reForm$nbCovs
  nbRecovs <- reForm$nbRecovs
  nbG0covs <- reForm$nbG0covs
  
  gamInd<-(length(m$mod$estimate)-(nbCovs+1)*nbStates*(nbStates-1)*mixtures+1):(length(m$mod$estimate))-(ncol(m$covsPi)*(mixtures-1))-ifelse(nbRecovs,nbRecovs+1+nbG0covs+1,0)-ncol(m$covsDelta)*(nbStates-1)*(!m$conditions$stationary)*mixtures
  
  if(is.null(m$mod$Sigma)) stop("Covariance matrix is required -- has 'hessian' been set to FALSE for this model?")
  Sigma <- m$mod$Sigma
  
  if(!inherits(Sigma,"error")){
    
    tmpSplineInputs <- getSplineFormula(reForm$newformula,m$data,cbind(m$data,reForm$newdata))
    desMat <- stats::model.matrix(tmpSplineInputs$formula,data=tmpSplineInputs$covs)
    
    VC <- Sigma[gamInd[unique(c(m$conditions$betaCons))],gamInd[unique(c(m$conditions$betaCons))]]
    

  } else {
    VC <- matrix(NA,length(gamInd[unique(c(m$conditions$betaCons))]),length(gamInd[unique(c(m$conditions$betaCons))]))
  }
  return(VC)
}

# variance components functions
opt <- function(sigma2,W,betahat,X,df){
  D <- sigma2 * diag(length(betahat)) + W
  Dinv <- MASS::ginv(D)
  muhat <- MASS::ginv(t(X)%*%Dinv%*%X)%*%t(X)%*%Dinv%*%betahat
  r <- log((df - (t(betahat-X%*%muhat) %*% Dinv %*% (betahat-X%*%muhat)))^2)
  if(isTRUE(all.equal(Im(r),0))) r <- Re(r)
  return(r)
}

opt12 <- function(sigma2,W,betahat,X,df,opt="opt"){
  r <- rep(NA, 2)
  sig2 <- sigma2[1]
  if(sigma2[1]<0){
    sig2 <- 0
    H <- matrix(0,nrow(W),ncol(W)) #sigma$est * MASS::ginv(expm::sqrtm(D))
  } else {
    H <- MASS::ginv(expm::sqrtm(diag(length(betahat))+1/sigma2[1]*W))
  }
  if(any(is.complex(c(H)))){
    if(isTRUE(all.equal(Im(H),matrix(0,nrow(H),ncol(H))))) H <- Re(H)
  }
  D <- sig2 * diag(length(betahat)) + W
  Dinv <- MASS::ginv(D)
  muhat <- MASS::ginv(t(X)%*%Dinv%*%X)%*%t(X)%*%Dinv%*%betahat
  ztilde <- H %*% (betahat-X%*%muhat) + X %*% muhat # + (diag(length(betahat))-H) %*% X %*% muhat
  
  Dinvorig <- MASS::ginv(sigma2[1] * diag(length(betahat)) + W)
  muhatorig <- MASS::ginv(t(X)%*%Dinvorig%*%X)%*%t(X)%*%Dinvorig%*%betahat
  r[1] <- df - t(betahat-X%*%muhatorig) %*% Dinvorig %*% (betahat-X%*%muhatorig)
  r[2] <- (t((ztilde-X%*%muhat))%*%(ztilde-X%*%muhat))/df - sig2
  if(isTRUE(all.equal(Im(r),rep(0,length(r))))) r <- Re(r)
  if(opt=="solve") return(r)
  else return(log(sum(r^2)))
}
