#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom stats var qt
MI_summary<-function(im,alpha=0.95,ncores=4){
  
  simind <- which((unlist(lapply(im,is.momentuHMM))))
  nsims <- length(simind)
  if(nsims<1) stop("'im' must be a list comprised of momentuHMM objects")
  
  checkmove <- which(!(unlist(lapply(im,is.momentuHMM))))
  if(length(checkmove)) {
    im[checkmove]<-NULL
    warning("The following imputations are not momentuHMM objects and will be ignored: ",paste(checkmove,collapse=", "))
  }
  checksims <- lapply(im,function(x) x[match(c("conditions","bounds"),names(x))])
  ident <- !unlist(lapply(checksims,function(x) isTRUE(all.equal(x,checksims[[1]]))))
  if(any(ident)) stop("Model conditions and bounds for each imputation must be identical. Imputations that do not match the first: ",paste(which(ident),collapse=", "))
  
  m <- im[[1]]
  data <- m$data
  nbStates <- length(m$stateNames)
  dist <- m$conditions$dist
  distnames <- names(dist)
  estAngleMean <- m$conditions$estAngleMean
  zeroInflation <- m$conditions$zeroInflation
  DM <- m$conditions$DM
  
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,m$conditions$bounds,DM)
  
  parms <- names(m$mle)[which(!unlist(lapply(m$mle,is.null)))]
  nparms <- length(parms)-(!is.null(m$mle$gamma))
  xmat <- xbar <- xvar <- W_m <- B_m <- MI_se <- lower <- upper <- list()
  parmcols <- unlist(lapply(m$mle,function(x) ncol(x)))
  

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
  #get variances for each parameter on natural scale
  if(is.null(m$CI_real) | is.null(m$CI_beta)){
    registerDoParallel(cores=ncores)
    if(nsims>5) cat("Calculating standard errors on real scale; this might take a while for large models or simulations...")
    CI <- foreach(i = 1:nsims) %dopar% {
      ci<-momentuHMM::CI_real(im[[i]])
      ci$beta$se<-momentuHMM::CI_beta(im[[i]])$beta$se
      ci
    }
    if(nsims>5) cat("DONE\n")
    stopImplicitCluster()
  } else {
    CI <- lapply(im,function(x) x$CI_real)
    for(i in 1:nsims){
      CI[[i]]$beta$se <- im[[i]]$CI_beta$beta$se
    }
  }
  
  for(parm in 1:nparms){
    
    parnames <- rownames(m$mle[[parms[parm]]])
    if(!is.null(dist[[parms[parm]]]))
      if(dist[[parms[parm]]] %in% c("wrpcauchy","vm") & !estAngleMean[[parms[parm]]])
        parnames <- rownames(m$mle[[parms[parm]]])[-1]

    
    if(!is.null(parnames)){
      xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
      rownames(xbar[[parms[parm]]]) <- parnames
      MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    }
    
    for(j in parnames){
      if(is.character(parnames)){
        xmat[[parms[parm]]][[j]] <- matrix(unlist(lapply(im,function(x) x$mle[[parms[parm]]][j,])) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
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
  lower[["delta"]] <- xbar[["delta"]]-quantSup*MI_se[["delta"]]
  upper[["delta"]] <- xbar[["delta"]]+quantSup*MI_se[["delta"]]   
    
  if(!is.null(m$mle$gamma)){
    xmat[["gamma"]] <- array(unlist(lapply(im,function(x) x$mle[["gamma"]])),c(nbStates,nbStates,nsims))
    xvar[["gamma"]] <- array(unlist(lapply(CI,function(x) x[["gamma"]]$se^2)),c(nbStates,nbStates,nsims))
    
    n <- apply(!(is.na(xmat[["gamma"]])+is.na(xvar[["gamma"]])),1:2,sum)
    
    if(any(n<2)) stop("need at least 2 simulations with valid estimates for gamma")
    
    xbar[["gamma"]] <-   apply( xmat[["gamma"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["gamma"]] <-   apply( xmat[["gamma"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["gamma"]] <- apply( xvar[["gamma"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["gamma"]] <- sqrt(W_m[["gamma"]] + (n+1)/n * B_m[["gamma"]])
    
    #MI_se_beta <- (1/xbar[["gamma"]]+1/(1-xbar[["gamma"]]))*MI_se[["gamma"]]  #get se on logit scale
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["gamma"]]/B_m[["gamma"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    #lower[["gamma"]] <- inv.logit(logit(xbar[["gamma"]])-quantSup*MI_se_beta)
    #upper[["gamma"]] <- inv.logit(logit(xbar[["gamma"]])+quantSup*MI_se_beta) 
    
    lower[["gamma"]] <- xbar[["gamma"]]-quantSup*MI_se[["gamma"]]
    upper[["gamma"]] <- xbar[["gamma"]]+quantSup*MI_se[["gamma"]]   
  }
  if(nbStates>1){
    xmat[["stateProbs"]] <- array(unlist(im_stateProbs),c(nrow(data),nbStates,nsims))
    xvar[["stateProbs"]] <- array(0,c(nrow(data),nbStates,nsims)) # don't have se's; might be a way to get these but probably quite complicated
    n <- apply(!(is.na(xmat[["stateProbs"]])+is.na(xvar[["stateProbs"]])),1:2,sum)
    
    if(any(n<2)) stop("need at least 2 simulations with valid estimates for stateProbs")
    
    xbar[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    B_m[["stateProbs"]] <-   apply( xmat[["stateProbs"]] , 1:2 , var,na.rm=TRUE)
    
    W_m[["stateProbs"]] <- apply( xvar[["stateProbs"]] , 1:2 , mean,na.rm=TRUE)
    MI_se[["stateProbs"]] <- sqrt(W_m[["stateProbs"]] + (n+1)/n * B_m[["stateProbs"]])
    
    #MI_se_beta <- (1/xbar[["stateProbs"]]+1/(1-xbar[["stateProbs"]]))*MI_se[["stateProbs"]]  #get se on logit scale
    
    dfs<-(n-1)*(1+1/(n+1)*W_m[["stateProbs"]]/B_m[["stateProbs"]])^2
    quantSup<-qt(1-(1-alpha)/2,df=dfs)
    
    #lower[["stateProbs"]] <- inv.logit(logit(xbar[["stateProbs"]])-quantSup*MI_se_beta)
    #upper[["stateProbs"]] <- inv.logit(logit(xbar[["stateProbs"]])+quantSup*MI_se_beta) 
    
    lower[["stateProbs"]] <- xbar[["stateProbs"]]-quantSup*MI_se[["stateProbs"]]
    upper[["stateProbs"]] <- xbar[["stateProbs"]]+quantSup*MI_se[["stateProbs"]]   
    
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
    
    lower[["timeInStates"]] <- xbar[["timeInStates"]]-quantSup*MI_se[["timeInStates"]]
    upper[["timeInStates"]] <- xbar[["timeInStates"]]+quantSup*MI_se[["timeInStates"]]   
  }
  Par <- list()
  for(i in distnames){
    if(dist[[i]] %in% c("wrpcauchy","vm")) {
      if(!estAngleMean[[i]]){
        xbar[[i]]<-matrix(c(rep(0,nbStates),xbar[[i]]),ncol=nbStates,byrow=T)
        lower[[i]]<-matrix(c(rep(NA,nbStates),lower[[i]]),ncol=nbStates,byrow=T)
        upper[[i]]<-matrix(c(rep(NA,nbStates),upper[[i]]),ncol=nbStates,byrow=T)  
        MI_se[[i]]<-matrix(c(rep(NA,nbStates),MI_se[[i]]),ncol=nbStates,byrow=T)
      }
    }
    Par[[i]] <- mi_parm_list(xbar[[i]],MI_se[[i]],lower[[i]],upper[[i]],m$mle[[i]])
  }
  
  Par$delta <- list(est=xbar$delta,se=MI_se$delta,lower=lower$delta,upper=upper$delta)

  if(nbStates>1) {
    Par$beta <- list(est=xbar$beta,se=MI_se$beta,lower=lower$beta,upper=upper$beta)
    rownames(Par$beta$est) <- rownames(m$mle$beta)
    rownames(Par$beta$se) <- rownames(m$mle$beta)
    rownames(Par$beta$lower) <- rownames(m$mle$beta)
    rownames(Par$beta$upper) <- rownames(m$mle$beta)
    colnames(Par$beta$est) <- colnames(m$mle$beta)
    colnames(Par$beta$se) <- colnames(m$mle$beta)
    colnames(Par$beta$lower) <- colnames(m$mle$beta)
    colnames(Par$beta$upper) <- colnames(m$mle$beta)
    
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
  for(i in distnames){
    if(dist[[i]] %in% c("wrpcauchy","vm")) {
      mh$data[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(mh$data[[i]]),byrow=TRUE),2,CircStats::circ.mean)
    } else {
      mh$data[[i]]<-apply(matrix(unlist(lapply(im,function(x) x$data[[i]])),ncol=length(mh$data[[i]]),byrow=TRUE),2,mean)
    }
  }
  if(all(c("x","y") %in% names(mh$data))){
    mh$data$x<-apply(matrix(unlist(lapply(im,function(x) x$data$x)),ncol=length(mh$data$x),byrow=TRUE),2,mean)
    mh$data$y<-apply(matrix(unlist(lapply(im,function(x) x$data$y)),ncol=length(mh$data$y),byrow=TRUE),2,mean)
  }
  
  if(!is.null(m$mle$gamma)){
    Par$gamma <- list(est=xbar$gamma,se=MI_se$gamma,lower=lower$gamma,upper=upper$gamma)
    rownames(Par$gamma$est) <- m$stateNames
    rownames(Par$gamma$se) <- m$stateNames
    rownames(Par$gamma$lower) <- m$stateNames
    rownames(Par$gamma$upper) <- m$stateNames
    colnames(Par$gamma$est) <- m$stateNames
    colnames(Par$gamma$se) <- m$stateNames
    colnames(Par$gamma$lower) <- m$stateNames
    colnames(Par$gamma$upper) <- m$stateNames
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