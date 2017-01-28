#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats var qt
#' @importFrom boot logit inv.logit
#' @importFrom CircStats circ.mean
MI_summary<-function(im,alpha=0.95,ncores=4){
  
  simind <- which((unlist(lapply(im,is.momentuHMM))))
  nsims <- length(simind)
  if(nsims<1) stop("'im' must be a list comprised of momentuHMM objects")
  
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
  
  if(is.null(m$CI_real) | is.null(m$CI_beta)){
    #get variances for each parameter on natural scale
    if(is.null(m$CI_real)){
      registerDoParallel(cores=ncores)
      if(nsims>5) cat("Calculating standard errors on natural scale; this might take a while for large models or simulations...")
      CI <- foreach(i = 1:nsims) %dopar% {
        ci<-momentuHMM::CI_real(im[[i]])
      }
      if(nsims>5) cat("DONE\n")
      stopImplicitCluster()
    }
    #get variances for each parameter on working scale
    if(is.null(m$CI_beta)){
      registerDoParallel(cores=ncores)
      if(nsims>5) cat("Calculating standard errors on working scale; this might take a while for large models or simulations...")
      CIbeta <- foreach(i = 1:nsims) %dopar% {
        ci<-momentuHMM::CI_beta(im[[i]])
      }
      if(nsims>5) cat("DONE\n")
      stopImplicitCluster()
    }
  } else {
    CI <- lapply(im,function(x) x$CI_real)
    CIbeta <- lapply(im,function(x) x$CI_beta)
  }
  
  #DMind$beta<-FALSE
  
  # pool estimates on real scale
  parms <- distnames
  nparms <- length(parms)
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parmcols <- rep(nbStates,nparms)
  
  for(parm in 1:nparms){
    
    parnames <- rownames(CI[[1]][[parms[parm]]]$est)
    #parnames <- rownames(m$mle[[parms[parm]]])
    
    #if(!is.null(dist[[parms[parm]]]))
    #  if(dist[[parms[parm]]] %in% angledists & !estAngleMean[[parms[parm]]] & DMind[[parms[parm]]])
    #    parnames <- rownames(m$mle[[parms[parm]]])[-1]
      
    xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
    rownames(xbar[[parms[parm]]]) <- parnames
    MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    
    for(j in parnames){
      
      xmat[[parms[parm]]][[j]] <- matrix(unlist(lapply(CI,function(x) x[[parms[parm]]]$est[j,])) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
      xvar[[parms[parm]]][[j]] <- matrix(unlist(lapply(CI,function(x) x[[parms[parm]]]$se[j,]^2)) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)

      n <- apply(!(is.na(xmat[[parms[parm]]][[j]])+is.na(xvar[[parms[parm]]][[j]])),2,sum)
      
      if(any(n<2)) stop("need at least 2 simulations with valid estimates for ",parms[parm]," ",j)
      
      xbar[[parms[parm]]][j,] <- apply(xmat[[parms[parm]]][[j]],2,mean,na.rm=TRUE)
      B_m[[parms[parm]]][[j]] <- apply(xmat[[parms[parm]]][[j]],2,var,na.rm=TRUE)
      
      W_m[[parms[parm]]][[j]] <- apply(xvar[[parms[parm]]][[j]],2,mean,na.rm=TRUE)
      MI_se[[parms[parm]]][j,] <- sqrt(W_m[[parms[parm]]][[j]] + (n+1)/n * B_m[[parms[parm]]][[j]])
      
      dfs<-(n-1)*(1+1/(n+1)*W_m[[parms[parm]]][[j]]/B_m[[parms[parm]]][[j]])^2
      quantSup<-qt(1-(1-alpha)/2,df=dfs)
      lower[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]-quantSup*MI_se[[parms[parm]]][j,]
      upper[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]+quantSup*MI_se[[parms[parm]]][j,]   
      
    }
  }
  
  xbar[["delta"]] <- matrix(NA,nrow=1,ncol=nbStates)
  MI_se[["delta"]] <- lower[["delta"]] <- upper[["delta"]] <- xbar[["delta"]]

  xmat[["delta"]] <- matrix(unlist(lapply(im,function(x) x$mle[["delta"]])) , ncol=nbStates, nrow=nsims, byrow=TRUE)
  xvar[["delta"]] <- matrix(unlist(lapply(CI,function(x) x[["delta"]]$se^2)) , ncol=nbStates, nrow=nsims, byrow=TRUE)
  n <- apply(!(is.na(xmat[["delta"]])+is.na(xvar[["delta"]])),2,sum)
  
  if(any(n<2)) stop("need at least 2 simulations with valid estimates for delta")
  
  xbar[["delta"]] <- apply(xmat[["delta"]],2,mean,na.rm=TRUE)
  B_m[["delta"]] <- apply(xmat[["delta"]],2,var,na.rm=TRUE)
  
  W_m[["delta"]] <- apply(xvar[["delta"]],2,mean,na.rm=TRUE)
  MI_se[["delta"]] <- sqrt(W_m[["delta"]] + (n+1)/n * B_m[["delta"]])
  
  dfs<-(n-1)*(1+1/(n+1)*W_m[["delta"]]/B_m[["delta"]])^2
  quantSup<-qt(1-(1-alpha)/2,df=dfs)
  #lower[["delta"]] <- xbar[["delta"]]-quantSup*MI_se[["delta"]]
  #upper[["delta"]] <- xbar[["delta"]]+quantSup*MI_se[["delta"]] 
  lower[["delta"]] <- probCI(xbar[["delta"]],MI_se[["delta"]],quantSup,"lower")
  upper[["delta"]] <- probCI(xbar[["delta"]],MI_se[["delta"]],quantSup,"upper")
    
  if(!is.null(m$mle$gamma)){
    xmat[["gamma"]] <- array(unlist(lapply(im,function(x) x$mle[["gamma"]])),c(nbStates,nbStates,nsims))
    xvar[["gamma"]] <- array(unlist(lapply(CI,function(x) x[["gamma"]]$se^2)),c(nbStates,nbStates,nsims))
    
    n <- apply(!(is.na(xmat[["gamma"]])+is.na(xvar[["gamma"]])),1:2,sum)
    
    if(any(n<2)) stop("need at least 2 simulations with valid estimates for gamma")
    
    xbar[["gamma"]] <-   apply( xmat[["gamma"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["gamma"]] <-   apply( xmat[["gamma"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["gamma"]] <- apply( xvar[["gamma"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["gamma"]] <- sqrt(W_m[["gamma"]] + (n+1)/n * B_m[["gamma"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["gamma"]]/B_m[["gamma"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
 
    #lower[["gamma"]] <- xbar[["gamma"]]-quantSup*MI_se[["gamma"]]
    #upper[["gamma"]] <- xbar[["gamma"]]+quantSup*MI_se[["gamma"]]    
    lower[["gamma"]] <- probCI(xbar[["gamma"]],MI_se[["gamma"]],quantSup,"lower")
    upper[["gamma"]] <- probCI(xbar[["gamma"]],MI_se[["gamma"]],quantSup,"upper")
  }
  Par <- list()
  Par$real <- list()
  for(i in distnames){
    pnames <- CI[[1]][[i]]$est
    if(dist[[i]] %in% angledists) {
      if(!estAngleMean[[i]]){
        xbar[[i]]<-matrix(c(rep(0,nbStates),xbar[[i]]),ncol=nbStates,byrow=T)
        lower[[i]]<-matrix(c(rep(NA,nbStates),lower[[i]]),ncol=nbStates,byrow=T)
        upper[[i]]<-matrix(c(rep(NA,nbStates),upper[[i]]),ncol=nbStates,byrow=T)  
        MI_se[[i]]<-matrix(c(rep(NA,nbStates),MI_se[[i]]),ncol=nbStates,byrow=T)
        tmp<- rbind(matrix(rep(0,nbStates),1,nbStates),pnames)
        rownames(tmp)[1]<-"mean"
        pnames <- tmp
        
      }
    }
    Par$real[[i]] <- mi_parm_list(xbar[[i]],MI_se[[i]],lower[[i]],upper[[i]],pnames)
  }
  
  if(!is.null(m$mle$gamma)){
    Par$real$gamma <- list(est=xbar$gamma,se=MI_se$gamma,lower=lower$gamma,upper=upper$gamma)
    rownames(Par$real$gamma$est) <- m$stateNames
    rownames(Par$real$gamma$se) <- m$stateNames
    rownames(Par$real$gamma$lower) <- m$stateNames
    rownames(Par$real$gamma$upper) <- m$stateNames
    colnames(Par$real$gamma$est) <- m$stateNames
    colnames(Par$real$gamma$se) <- m$stateNames
    colnames(Par$real$gamma$lower) <- m$stateNames
    colnames(Par$real$gamma$upper) <- m$stateNames
  }
  
  Par$real$delta <- list(est=xbar$delta,se=MI_se$delta,lower=lower$delta,upper=upper$delta)
  names(Par$real$delta$est) <- names(Par$real$delta$se) <- names(Par$real$delta$lower) <- names(Par$real$delta$upper) <- m$stateNames
  
  # pool estimates on working scale
  parms <- names(CIbeta[[1]])
  nparms <- length(parms)
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parmcols <- lapply(m$conditions$fullDM,ncol)#
  parmcols$beta <- ncol(m$mle$beta)#unlist(lapply(m$mle,function(x) ncol(x)))
  parmcols <- unlist(parmcols[parms])
  
  for(parm in 1:nparms){
    
    parnames <- rownames(CIbeta[[1]][[parms[parm]]]$est)
    #parnames <- rownames(m$mle[[parms[parm]]])
    
    #if(!is.null(dist[[parms[parm]]]))
    #  if(dist[[parms[parm]]] %in% angledists & !estAngleMean[[parms[parm]]] & DMind[[parms[parm]]])
    #    parnames <- rownames(m$mle[[parms[parm]]])[-1]
    
    xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
    rownames(xbar[[parms[parm]]]) <- parnames
    MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    
    for(j in parnames){
      
      xmat[[parms[parm]]][[j]] <- matrix(unlist(lapply(CIbeta,function(x) x[[parms[parm]]]$est[j,])) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
      
      #if(DMind[[parms[parm]]])
      #  xvar[[parms[parm]]][[j]] <- matrix(unlist(lapply(CI,function(x) x[[parms[parm]]]$se[j,]^2)) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
      #else
      xvar[[parms[parm]]][[j]] <- matrix(unlist(lapply(CIbeta,function(x) x[[parms[parm]]]$se[j,]^2)) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
      
      n <- apply(!(is.na(xmat[[parms[parm]]][[j]])+is.na(xvar[[parms[parm]]][[j]])),2,sum)
      
      if(any(n<2)) stop("need at least 2 simulations with valid estimates for ",parms[parm]," ",j)
      
      xbar[[parms[parm]]][j,] <- apply(xmat[[parms[parm]]][[j]],2,mean,na.rm=TRUE)
      B_m[[parms[parm]]][[j]] <- apply(xmat[[parms[parm]]][[j]],2,var,na.rm=TRUE)
      
      W_m[[parms[parm]]][[j]] <- apply(xvar[[parms[parm]]][[j]],2,mean,na.rm=TRUE)
      MI_se[[parms[parm]]][j,] <- sqrt(W_m[[parms[parm]]][[j]] + (n+1)/n * B_m[[parms[parm]]][[j]])
      
      dfs<-(n-1)*(1+1/(n+1)*W_m[[parms[parm]]][[j]]/B_m[[parms[parm]]][[j]])^2
      quantSup<-qt(1-(1-alpha)/2,df=dfs)
      lower[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]-quantSup*MI_se[[parms[parm]]][j,]
      upper[[parms[parm]]][j,] <- xbar[[parms[parm]]][j,]+quantSup*MI_se[[parms[parm]]][j,]   
      
    }
  }
  
  Par$beta <- list()
  for(i in parms){
    Par$beta[[i]] <- mi_parm_list(xbar[[i]],MI_se[[i]],lower[[i]],upper[[i]],CIbeta[[1]][[i]]$est)
  }
  
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  
  if(nbStates>1){
    xmat[["stateProbs"]] <- array(unlist(im_stateProbs),c(nrow(data),nbStates,nsims))
    xvar[["stateProbs"]] <- array(0,c(nrow(data),nbStates,nsims)) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["stateProbs"]])+is.na(xvar[["stateProbs"]])),1:2,sum)
    
    if(any(n<2)) stop("need at least 2 simulations with valid estimates for stateProbs")
    
    xbar[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["stateProbs"]] <- apply( xvar[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["stateProbs"]] <- sqrt(W_m[["stateProbs"]] + (n+1)/n * B_m[["stateProbs"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["stateProbs"]]/B_m[["stateProbs"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    #lower[["stateProbs"]] <- xbar[["stateProbs"]]-quantSup*MI_se[["stateProbs"]]
    #upper[["stateProbs"]] <- xbar[["stateProbs"]]+quantSup*MI_se[["stateProbs"]]  
    lower[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"lower"))
    upper[["stateProbs"]] <- suppressWarnings(probCI(xbar[["stateProbs"]],MI_se[["stateProbs"]],quantSup,"upper"))
    
    xmat[["timeInStates"]] <- t(apply(im_states,1,function(x) {counts<-hist(x,breaks=seq(0.5,nbStates+0.5),plot=FALSE)$counts;counts/sum(counts)}))
    xvar[["timeInStates"]] <- matrix(0 , ncol=nbStates, nrow=nsims, byrow=TRUE) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["timeInStates"]])+is.na(xvar[["timeInStates"]])),2,sum)
    
    if(any(n<2)) stop("need at least 2 simulations with valid estimates for timeInStates")
    
    xbar[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,mean,na.rm=TRUE)
    B_m[["timeInStates"]] <- apply(xmat[["timeInStates"]],2,var,na.rm=TRUE)
    
    W_m[["timeInStates"]] <- apply(xvar[["timeInStates"]],2,mean,na.rm=TRUE)
    MI_se[["timeInStates"]] <- sqrt(W_m[["timeInStates"]] + (n+1)/n * B_m[["timeInStates"]])
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["timeInStates"]]/B_m[["timeInStates"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    #lower[["timeInStates"]] <- xbar[["timeInStates"]]-quantSup*MI_se[["timeInStates"]]
    #upper[["timeInStates"]] <- xbar[["timeInStates"]]+quantSup*MI_se[["timeInStates"]]  
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
  if(any(ident)) mh$conditions$fullDM <- NULL
  
  #average all numeric variables in imputed data
  for(i in distnames){
    if(dist[[i]] %in% angledists) {
      mh$data[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(mh$data[[i]]),byrow=TRUE),2,CircStats::circ.mean)
    } else if(dist[[i]] %in% "pois"){
      mh$data[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(mh$data[[i]]),byrow=TRUE),2,median)     
    } else {
      mh$data[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(mh$data[[i]]),byrow=TRUE),2,mean)
    }
  }
  for(j in names(mh$data)[which(unlist(lapply(mh$data,class))!="factor" & !(names(mh$data) %in% distnames))]){
    mh$data[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(mh$data[[j]]),byrow=TRUE),2,mean)
    mh$data[[j]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[j]])),ncol=length(mh$data[[j]]),byrow=TRUE),2,mean)
  }
  mh$Par <- Par
  return(momentuHMMMI(mh))
  
}

mi_parm_list<-function(est,se,lower,upper,m){
  Par <- list(est=est,se=se,lower=lower,upper=upper)
  rownames(Par$est) <- rownames(m)
  rownames(Par$se) <- rownames(m)
  rownames(Par$lower) <- rownames(m)
  rownames(Par$upper) <- rownames(m)
  colnames(Par$est) <- colnames(m)
  colnames(Par$se) <- colnames(m)
  colnames(Par$lower) <- colnames(m)
  colnames(Par$upper) <- colnames(m)
  Par
}

probCI<-function(x,se,z,bound="lower"){
  if(bound=="lower")
    ci<-boot::inv.logit(boot::logit(x)-z*(1/(x-x^2))*se) 
  else if(bound=="upper")
    ci<-boot::inv.logit(boot::logit(x)+z*(1/(x-x^2))*se) 
  ci
}