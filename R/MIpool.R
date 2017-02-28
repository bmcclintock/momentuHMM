#'
#' Calculate pooled estimates across multiple imputations
#' 
#' @param HMMfits List comprised of \code{\link{momentuHMM}} objects
#' 
#' @return A \code{\link{miSum}} object, i.e., a list comprised of model and pooled parameter summaries, including \code{data} (averaged across imputations), \code{conditions}, \code{Par}, and \code{MIcombine} 
#' (as returned by \code{\link[mitools]{MIcombine}} for working parameters).
#' 
#' \code{miSum$Par} is a list comprised of:
#' \item{beta}{Pooled estimates for the working parameters}
#' \item{real}{Pooled estimates for the natural parameters based on covariate means across imputations (if applicable)}
#' \item{timeInStates}{The proportion of time steps assigned to each state}
#' \item{states}{The most freqent state assignment for each time step based on the \code{\link{viterbi}} algorithm for each model fit}
#' \item{stateProbs}{Pooled state probability estimates for each time step}
#' 
#' @details
#' Pooled estimates, standard errors, and confidence intervals are calculated using standard multiple imputation formulas. Working scale parameters are pooled
#' using \code{\link[mitools]{MIcombine}} and t-distributed confidence intervals. Natural scale parameters are pooled using similar formulae but are based on 
#' covariate means across all imputations (if applicable) and normally-distributed confidence intervals.
#' 
#' Note that pooled estimates for \code{timeInStates} and \code{stateProbs} do not include within-model uncertainty and are based entirely on across-model variability.
#' 
#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats median var qt
#' @importFrom boot logit inv.logit
#' @importFrom CircStats circ.mean
#' @importFrom car dataEllipse
#' @importFrom mitools MIcombine
#' @importFrom MASS ginv
MIpool<-function(HMMfits,alpha=0.95,ncores){
  
  im <- HMMfits
  simind <- which((unlist(lapply(im,is.momentuHMM))))
  nsims <- length(simind)
  if(nsims<1) stop("'HMMfits' must be a list comprised of momentuHMM objects")
  
  checkmove <- which(!(unlist(lapply(im,is.momentuHMM))))
  if(length(checkmove)) {
    im[checkmove]<-NULL
    warning("The following imputations are not momentuHMM objects and will be ignored: ",paste(checkmove,collapse=", "))
  }
  checksims <- lapply(im,function(x) x[match("conditions",names(x))])
  ident <- !unlist(lapply(checksims,function(x) isTRUE(all.equal(x,checksims[[1]]))))
  if(any(ident)){
    # check that only differences are in the design matrix covariate values
    checksims2 <- lapply(checksims, function(x) x$conditions[-match("fullDM",names(x$conditions))])
    ident2 <- !unlist(lapply(checksims2,function(x) isTRUE(all.equal(x,checksims2[[1]]))))
    if(any(ident2)) stop("Model conditions for each imputation must be identical. Imputations that do not match the first: ",paste(which(ident),collapse=", "))
  }
  
  tmpDet <- which(unlist(lapply(im,function(x) det(x$mod$hessian)))==0)
  if(length(tmpDet)){
    warning("Hessian is singular for HMM fit(s): ",paste0(tmpDet,collapse=", "))
  }
  
  m <- im[[1]]
  data <- m$data
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  estAngleMean <- m$conditions$estAngleMean
  zeroInflation <- m$conditions$zeroInflation
  DM <- m$conditions$DM
  DMind <- m$conditions$DMind
  
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,DM,m$conditions$bounds)
  
  if(nbStates>1) {
    cat("Decoding state sequences and probabilities for each imputation... ")
    registerDoParallel(cores=ncores)
    im_states <- foreach(i = 1:nsims, .combine = rbind) %dopar% {momentuHMM::viterbi(im[[i]])}
    stopImplicitCluster()
    states <- apply(im_states,2,function(x) which.max(hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts))
    registerDoParallel(cores=ncores)
    im_stateProbs <- foreach(i = 1:nsims) %dopar% {momentuHMM::stateProbs(im[[i]])}
    stopImplicitCluster()
    cat("DONE\n")
  } else states <- rep(1,nrow(data))
  
  # pool estimates on working scale
  parms <- names(m$CI_beta)
  nparms <- length(parms)
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parmcols <- lapply(m$conditions$fullDM,ncol)#
  parmcols$beta <- ncol(m$mle$beta)
  parmcols <- unlist(parmcols[parms])
  
  parindex <- c(0,cumsum(c(unlist(lapply(m$conditions$fullDM,ncol)),length(m$mle$beta),nbStates-1)))
  names(parindex)[1:length(distnames)] <- distnames
  if(nbStates>1) names(parindex)[length(distnames)+1] <- "beta"
  names(parindex)[length(parindex)-1] <- "delta"
  
  miBeta <- mitools::MIcombine(results=lapply(im,function(x) x$mod$estimate),variances=lapply(im,function(x) ginv(x$mod$hessian)))
  
  for(parm in 1:nparms){
    
    parnames <- rownames(m$CI_beta[[parms[parm]]]$est)
    coeffs <- matrix(miBeta$coefficients[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    vars <- matrix(diag(miBeta$variance)[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    dfs <- matrix(miBeta$df[(parindex[parm]+1):parindex[parm+1]],nrow=length(parnames),dimnames=list(parnames))
    
    xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
    rownames(xbar[[parms[parm]]]) <- parnames
    MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    
    for(j in parnames){
      
      xbar[[parms[parm]]][j,] <- coeffs[j,]
      MI_se[[parms[parm]]][j,] <- sqrt(vars[j,])
      
      quantSup<-qt(1-(1-alpha)/2,df=dfs[j,])
      lower[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]-quantSup*MI_se[[parms[parm]]][j,]
      upper[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]+quantSup*MI_se[[parms[parm]]][j,]   
      
    }
  }
  
  Par <- list()
  Par$beta <- list()
  for(i in parms){
    Par$beta[[i]] <- mi_parm_list(xbar[[i]],MI_se[[i]],lower[[i]],upper[[i]],m$CI_beta[[i]]$est)
  }
  
  #average all numeric variables in imputed data
  mhdata<-m$data
  for(i in distnames){
    if(dist[[i]] %in% angledists) {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,CircStats::circ.mean)
    } else if(dist[[i]] %in% "pois"){
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,median)     
    } else {
      mhdata[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(m$data[[i]]),byrow=TRUE),2,mean)
    }
  }
  for(j in names(m$data)[which(unlist(lapply(m$data,class))!="factor" & !(names(m$data) %in% distnames))]){
    mhdata[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(m$data[[j]]),byrow=TRUE),2,mean)
  }
  mhrawCovs<-m$rawCovs
  for(j in names(m$rawCovs)[which(unlist(lapply(m$rawCovs,class))!="factor")]){
    mhrawCovs[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$rawCovs[[j]])),ncol=length(m$rawCovs[[j]]),byrow=TRUE),2,mean)
  }
  
  tempCovs <- mhdata[1,]
  for(j in names(m$data)[which(unlist(lapply(m$data,class))!="factor")]){
    tempCovs[[j]]<-mean(mhdata[[j]],na.rm=TRUE)
  }
  
  tmPar <- m$mle[distnames]
  parindex <- c(0,cumsum(unlist(lapply(m$conditions$fullDM,ncol)))[-length(m$conditions$fullDM)])
  names(parindex) <- distnames
  for(i in distnames){
    if(!is.null(m$conditions$DM[[i]])){# & m$conditions$DMind[[i]]){
      tmPar[[i]] <- m$mod$estimate[parindex[[i]]+1:ncol(m$conditions$fullDM[[i]])]
      names(tmPar[[i]])<-colnames(m$conditions$fullDM[[i]])
    } else if((m$conditions$dist[[i]] %in% angledists) & (!m$conditions$estAngleMean[[i]])){
      tmPar[[i]] <- tmPar[[i]][-1,]
    }
  }
  
  inputs <- checkInputs(nbStates,m$conditions$dist,tmPar,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$zeroInflation,m$conditions$DM,m$conditions$userBounds,m$conditions$cons,m$conditions$workcons,m$stateNames)
  p<-inputs$p
  DMinputs<-getDM(tempCovs,inputs$DM,m$conditions$dist,nbStates,p$parNames,p$bounds,tmPar,m$conditions$cons,m$conditions$workcons,m$conditions$zeroInflation,m$conditions$circularAngleMean)
  fullDM<-DMinputs$fullDM
  nbCovs <- ncol(model.matrix(m$conditions$formula,m$data))-1 # substract intercept column
  
  Par$real<-list()
  for(i in distnames){
    
    DMind[[i]] <- FALSE
    par <- c(w2n(miBeta$coefficients,p$bounds,p$parSize,nbStates,nbCovs,m$conditions$estAngleMean,m$conditions$circularAngleMean,m$conditions$stationary,m$conditions$cons,fullDM,DMind,m$conditions$workcons,1,dist[i],m$conditions$Bndind)[[i]])

    if(!(dist[[i]] %in% angledists) | (dist[[i]] %in% angledists & m$conditions$estAngleMean[[i]] & !m$conditions$Bndind[[i]])) {
      Par$real[[i]] <- get_CI(miBeta$coefficients,par,m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],nbStates,alpha,p$parNames[[i]],m$stateNames)
    } else {
      if(!m$conditions$estAngleMean[[i]]){
        Par$real[[i]] <- get_CI(miBeta$coefficients,par[-(1:nbStates)],m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],nbStates,alpha,p$parNames[[i]],m$stateNames)
        Par$real[[i]]$est <- matrix(c(rep(0,nbStates),Par$real[[i]]$est),ncol=nbStates,byrow=T)
        Par$real[[i]]$se <- matrix(c(rep(NA,nbStates),Par$real[[i]]$se),ncol=nbStates,byrow=T)
        Par$real[[i]]$lower <- matrix(c(rep(NA,nbStates),Par$real[[i]]$lower),ncol=nbStates,byrow=T)
        Par$real[[i]]$upper <- matrix(c(rep(NA,nbStates),Par$real[[i]]$upper),ncol=nbStates,byrow=T)  
        dimnames(Par$real[[i]]$est) <- dimnames(Par$real[[i]]$se) <- dimnames(Par$real[[i]]$lower) <- dimnames(Par$real[[i]]$upper) <- list(c("mean",p$parNames[[i]]),m$stateNames)
      } else {
        if(m$conditions$Bndind[[i]]){
          Par$real[[i]] <- CI_angle(miBeta$coefficients,par,m,parindex[[i]]+1:ncol(fullDM[[i]]),fullDM[[i]],DMind[[i]],p$bounds[[i]],m$conditions$cons[[i]],m$conditions$workcons[[i]],miBeta$variance,m$conditions$circularAngleMean[[i]],nbStates,alpha,p$parNames[[i]],m$stateNames)
        }
      }
    }
  }
  
  quantSup<-qnorm(1-(1-alpha)/2)
  
  # pooled delta estimates
  deltInd<-(length(miBeta$coefficients)-nbStates+2):length(miBeta$coefficients)
  est <- get_delta(miBeta$coefficients[deltInd],1:nbStates)
  lower<-upper<-se<-numeric(length(est))
  for(k in 1:length(est)){
    dN<-numDeriv::grad(get_delta,miBeta$coefficients[deltInd],i=k)
    se[k]<-suppressWarnings(sqrt(dN%*%miBeta$variance[deltInd,deltInd]%*%dN))
    lower[k] <- probCI(est[k],se[k],quantSup,bound="lower")
    upper[k] <- probCI(est[k],se[k],quantSup,bound="upper")
  }
  Par$real$delta <- list(est=est,se=se,lower=lower,upper=upper)
  names(Par$real$delta$est) <- names(Par$real$delta$se) <- names(Par$real$delta$lower) <- names(Par$real$delta$upper) <- m$stateNames
    
  # pooled gamma estimates
  if(nbStates>1){
    gamInd<-(length(miBeta$coefficients)-(nbCovs+1)*nbStates*(nbStates-1)+1):(length(miBeta$coefficients))-(nbStates-1)
    est <- get_gamma(matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),model.matrix(m$conditions$formula,mhdata),nbStates,1:nbStates,1:nbStates)
    lower<-upper<-se<-matrix(0,nrow(est),ncol(est))
    for(i in 1:nrow(est)){
      for(j in 1:ncol(est)){
        dN<-numDeriv::grad(get_gamma,matrix(miBeta$coefficients[gamInd],nrow=nbCovs+1),covs=model.matrix(m$conditions$formula,mhdata),nbStates=nbStates,i=i,j=j)
        se[i,j]<-suppressWarnings(sqrt(dN%*%miBeta$variance[gamInd,gamInd]%*%dN))
        lower[i,j]<-1/(1+exp(-(log(est[i,j]/(1-est[i,j]))-quantSup*(1/(est[i,j]-est[i,j]^2))*se[i,j])))#est[i,j]-quantSup*se[i,j]
        upper[i,j]<-1/(1+exp(-(log(est[i,j]/(1-est[i,j]))+quantSup*(1/(est[i,j]-est[i,j]^2))*se[i,j])))#est[i,j]+quantSup*se[i,j]
      }
    }
    Par$real$gamma <- list(est=est,se=se,lower=lower,upper=upper)
    dimnames(Par$real$gamma$est) <- dimnames(Par$real$gamma$se) <- dimnames(Par$real$gamma$lower) <- dimnames(Par$real$gamma$upper) <- list(m$stateNames,m$stateNames)
  }
  
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  
  if(nbStates>1){
    xmat[["stateProbs"]] <- array(unlist(im_stateProbs),c(nrow(data),nbStates,nsims))
    xvar[["stateProbs"]] <- array(0,c(nrow(data),nbStates,nsims)) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["stateProbs"]])+is.na(xvar[["stateProbs"]])),1:2,sum)
    
    if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates for stateProbs")
    
    xbar[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["stateProbs"]] <- apply( xvar[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["stateProbs"]] <- sqrt(W_m[["stateProbs"]] + (n+1)/n * B_m[["stateProbs"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["stateProbs"]]/B_m[["stateProbs"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    lower[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"lower"))
    upper[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"upper"))
    
    xmat[["timeInStates"]] <- t(apply(im_states,1,function(x) {counts<-hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts;counts/sum(counts)}))
    xvar[["timeInStates"]] <- matrix(0 , ncol=nbStates, nrow=nsims, byrow=TRUE) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["timeInStates"]])+is.na(xvar[["timeInStates"]])),2,sum)
    
    if(any(n<2)) warning("need at least 2 simulations with valid point and variance estimates for timeInStates")
    
    xbar[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,mean,na.rm=TRUE)
    B_m[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,var,na.rm=TRUE)
    
    W_m[["timeInStates"]] <- apply(xvar[["timeInStates"]],2,mean,na.rm=TRUE)
    MI_se[["timeInStates"]] <- sqrt(W_m[["timeInStates"]] + (n+1)/n * B_m[["timeInStates"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["timeInStates"]]/B_m[["timeInStates"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    lower[["timeInStates"]] <- probCI(xbar[["timeInStates"]],MI_se[["timeInStates"]],quantSup,"lower")
    upper[["timeInStates"]] <- probCI(xbar[["timeInStates"]],MI_se[["timeInStates"]],quantSup,"upper")
  }

  if(nbStates>1) {
    Par$timeInStates <- list(est=xbar$timeInStates,se=MI_se$timeInStates,lower=lower$timeInStates,upper=upper$timeInStates)
    names(Par$timeInStates$est) <- m$stateNames
    names(Par$timeInStates$se) <- m$stateNames
    names(Par$timeInStates$lower) <- m$stateNames
    names(Par$timeInStates$upper) <- m$stateNames
    
    Par$states <- states
    
    Par$stateProbs <- list(est=xbar$stateProbs,se=MI_se$stateProbs,lower=lower$stateProbs,upper=upper$stateProbs)
    rownames(Par$stateProbs$est) <- data$ID
    rownames(Par$stateProbs$se) <- data$ID
    rownames(Par$stateProbs$lower) <- data$ID
    rownames(Par$stateProbs$upper) <- data$ID
    colnames(Par$stateProbs$est) <- m$stateNames
    colnames(Par$stateProbs$se) <- m$stateNames
    colnames(Par$stateProbs$lower) <- m$stateNames
    colnames(Par$stateProbs$upper) <- m$stateNames
  }
  
  mh <- im[[1]]
  attr(mh,"class") <- NULL
  mh$mle <- NULL
  mh$mod <- NULL
  mh$CI_real <- NULL
  mh$CI_beta <- NULL
  if(any(ident)) mh$conditions$fullDM <- fullDM
  
  mh$data<-mhdata
  mh$rawCovs<-mhrawCovs

  errorEllipse<-NULL
  if(all(c("x","y") %in% names(mh$data))){
    checkerrs <- lapply(im,function(x) x$data[match(c("x","y"),names(x$data))])
    ident <- !unlist(lapply(checkerrs,function(x) isTRUE(all.equal(x,checkerrs[[1]]))))
    if(any(ident)){
      # calculate location alpha% error ellipses
      cat("Calculating location",paste0(alpha*100,"%"),"error ellipses... ")
      registerDoParallel(cores=ncores)
      errorEllipse<-foreach(i = 1:nrow(mh$data)) %dopar% {
        car::dataEllipse(cbind(unlist(lapply(im,function(x) x$data$x[i])),unlist(lapply(im,function(x) x$data$y[i]))),levels=alpha,draw=FALSE,segments=100)
      }
      stopImplicitCluster()
      cat("DONE\n")
    }
  }
  mh$errorEllipse <- errorEllipse
  mh$Par <- Par
  mh$MIcombine <- miBeta
  
  return(miSum(mh))
}

mi_parm_list<-function(est,se,lower,upper,m){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  dimnames(Par$est) <- dimnames(Par$se) <- dimnames(Par$lower) <- dimnames(Par$upper) <- list(rownames(m),colnames(m))
  Par
}

probCI<-function(x,se,z,bound="lower"){
  if(bound=="lower")
    ci<-boot::inv.logit(boot::logit(x)-z*(1/(x-x^2))*se) 
  else if(bound=="upper")
    ci<-boot::inv.logit(boot::logit(x)+z*(1/(x-x^2))*se) 
  ci
}