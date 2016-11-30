#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
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
  nbStates <- ncol(m$mle$stepPar)
  p <- parDef(m$conditions$stepDist,m$conditions$angleDist,m$conditions$omegaDist,m$conditions$dryDist,m$conditions$diveDist,m$conditions$iceDist,m$conditions$landDist,nbStates,m$conditions$estAngleMean,
              m$conditions$zeroInflation,m$conditions$bounds,m$conditions$stepDM,m$conditions$angleDM,m$conditions$omegaDM,m$conditions$dryDM,m$conditions$diveDM,m$conditions$iceDM,m$conditions$landDM)
  
  parms <- names(m$mle)[which(!unlist(lapply(m$mle,is.null)))]
  nparms <- length(parms)-(!is.null(m$mle$gamma))
  nbStates <- ncol(m$mle$stepPar)
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
  } else states <- rep(1,nrow(m$data))
  #get variances for each parameter on natural scale
  if(is.null(m$CI_real) | is.null(m$CI_beta)){
    registerDoParallel(cores=ncores)
    if(nsims>5) cat("Calculating standard errors on real scale; this might take a while for large models or simulations...")
    CI <- foreach(i = 1:nsims) %dopar% {
      ci<-momentuHMM::CI_real(im[[i]])
      ci$beta$se<-momentuHMM::CI_beta(im[[i]])$betaPar$se
      ci
    }
    if(nsims>5) cat("DONE\n")
    stopImplicitCluster()
  } else {
    CI <- lapply(im,function(x) x$CI_real)
    for(i in 1:nsims){
      CI[[i]]$beta$se <- im[[i]]$CI_beta$betaPar$se
    }
  }
  
  for(parm in 1:nparms){
    
    if(parms[parm]=="anglePar"){
      if(!m$conditions$estAngleMean) parnames <- rownames(m$mle[[parms[parm]]])[-1]
    } else parnames <- rownames(m$mle[[parms[parm]]])
    
    if(!is.null(parnames)){
      xbar[[parms[parm]]] <- matrix(NA,nrow=length(parnames),ncol=parmcols[parm])
      rownames(xbar[[parms[parm]]]) <- parnames
      MI_se[[parms[parm]]] <- lower[[parms[parm]]] <- upper[[parms[parm]]] <- xbar[[parms[parm]]]
    }
    
    for(i in parnames){
      if(is.character(parnames)){
        xmat[[parms[parm]]][[i]] <- matrix(unlist(lapply(im,function(x) x$mle[[parms[parm]]][i,])) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
        xvar[[parms[parm]]][[i]] <- matrix(unlist(lapply(CI,function(x) x[[parms[parm]]]$se[i,]^2)) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
        n <- apply(!(is.na(xmat[[parms[parm]]][[i]])+is.na(xvar[[parms[parm]]][[i]])),2,sum)
        
        if(any(n<2)) stop("need at least 2 simulations with valid estimates for ",parms[parm]," ",i)
        
        xbar[[parms[parm]]][i,] <- apply(xmat[[parms[parm]]][[i]],2,mean,na.rm=TRUE)
        B_m[[parms[parm]]][[i]] <- apply(xmat[[parms[parm]]][[i]],2,var,na.rm=TRUE)
        
        W_m[[parms[parm]]][[i]] <- apply(xvar[[parms[parm]]][[i]],2,mean,na.rm=TRUE)
        MI_se[[parms[parm]]][i,] <- sqrt(W_m[[parms[parm]]][[i]] + (n+1)/n * B_m[[parms[parm]]][[i]])
        
        dfs<-(n-1)*(1+1/(n+1)*W_m[[parms[parm]]][[i]]/B_m[[parms[parm]]][[i]])^2
        quantSup<-qt(1-(1-alpha)/2,df=dfs)
        lower[[parms[parm]]][i,] <- xbar[[parms[parm]]][i,]-quantSup*MI_se[[parms[parm]]][i,]
        upper[[parms[parm]]][i,] <- xbar[[parms[parm]]][i,]+quantSup*MI_se[[parms[parm]]][i,]   
        
      } #else { #delta parameter not currently handled by CI_real
      #xmat[[parms[parm]]] <- matrix(unlist(lapply(im,function(x) x$mle[[parms[parm]]][i])) , ncol=nbStates, nrow=nsims, byrow=TRUE) 
      #xvar[[parms[parm]]] <- matrix(unlist(lapply(CI,function(x) x[[parms[parm]]]$se[i]^2)) , ncol=parmcols[parm], nrow=nsims, byrow=TRUE)
      #n <- apply(!(is.na(xmat[[parms[parm]]])+is.na(xvar[[parms[parm]]])),2,sum)
      
      #if(any(n<2)) stop("need at least 2 simulations with valid estimates for ",parms[parm])
      
      #xbar[[parms[parm]]] <- apply(xmat[[parms[parm]]],2,mean,na.rm=TRUE)
      #B_m[[parms[parm]]] <- apply(xmat[[parms[parm]]],2,var,na.rm=TRUE)
      
      #W_m[[parms[parm]]] <- apply(xvar[[parms[parm]]],2,mean,na.rm=TRUE)
      #MI_se[[parms[parm]]] <- sqrt(W_m[[parms[parm]]] + (n+1)/n * B_m[[parms[parm]]])
      
      #dfs<-(n-1)*(1+1/(n+1)*W_m[[parms[parm]]]/B_m[[parms[parm]]])^2
      #quantSup<-qt(1-(1-alpha)/2,df=dfs)
      #lower[[parms[parm]]]<- xbar[[parms[parm]]]-quantSup*MI_se[[parms[parm]]]
      #upper[[parms[parm]]] <- xbar[[parms[parm]]]+quantSup*MI_se[[parms[parm]]]
      #}
    }
  }
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
    xmat[["stateProbs"]] <- array(unlist(im_stateProbs),c(nrow(m$data),nbStates,nsims))
    xvar[["stateProbs"]] <- array(0,c(nrow(m$data),nbStates,nsims)) # don't have se's; might be a way to get these but probably quite complicated
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
  stepPar <- mi_parm_list(xbar$stepPar,MI_se$stepPar,lower$stepPar,upper$stepPar,m$mle$stepPar)
  anglePar <- omegaPar <- dryPar <- divePar <- icePar <- landPar <- beta <- stateProbs <- timeInStates <- NULL
  if(m$conditions$angleDist!="none") {
    if(!m$conditions$estAngleMean){
      xbar$anglePar<-matrix(c(rep(0,nbStates),xbar$anglePar),ncol=nbStates,byrow=T)
      lower$anglePar<-matrix(c(rep(NA,nbStates),lower$anglePar),ncol=nbStates,byrow=T)
      upper$anglePar<-matrix(c(rep(NA,nbStates),upper$anglePar),ncol=nbStates,byrow=T)  
      MI_se$anglePar<-matrix(c(rep(NA,nbStates),MI_se$anglePar),ncol=nbStates,byrow=T)
    }
    anglePar <- mi_parm_list(xbar$anglePar,MI_se$anglePar,lower$anglePar,upper$anglePar,m$mle$anglePar)
  }
  if(m$conditions$omegaDist!="none") omegaPar <- mi_parm_list(xbar$omegaPar,MI_se$omegaPar,lower$omegaPar,upper$omegaPar,m$mle$omegaPar)
  if(m$conditions$dryDist!="none") dryPar <- mi_parm_list(xbar$dryPar,MI_se$dryPar,lower$dryPar,upper$dryPar,m$mle$dryPar)  
  if(m$conditions$diveDist!="none") divePar <- mi_parm_list(xbar$divePar,MI_se$divePar,lower$divePar,upper$divePar,m$mle$divePar)  
  if(m$conditions$iceDist!="none") icePar <- mi_parm_list(xbar$icePar,MI_se$icePar,lower$icePar,upper$icePar,m$mle$icePar)
  if(m$conditions$landDist!="none") landPar <- mi_parm_list(xbar$landPar,MI_se$landPar,lower$landPar,upper$landPar,m$mle$landPar)
  
  if(nbStates>1) {
    beta <- list(est=xbar$beta,se=MI_se$beta,lower=lower$beta,upper=upper$beta)
    rownames(beta$est) <- rownames(m$mle$beta)
    rownames(beta$se) <- rownames(m$mle$beta)
    rownames(beta$lower) <- rownames(m$mle$beta)
    rownames(beta$upper) <- rownames(m$mle$beta)
    colnames(beta$est) <- colnames(m$mle$beta)
    colnames(beta$se) <- colnames(m$mle$beta)
    colnames(beta$lower) <- colnames(m$mle$beta)
    colnames(beta$upper) <- colnames(m$mle$beta)
    
    timeInStates <- list(est=xbar$timeInStates,se=MI_se$timeInStates,lower=lower$timeInStates,upper=upper$timeInStates)
    names(timeInStates$est) <- m$stateNames
    names(timeInStates$se) <- m$stateNames
    names(timeInStates$lower) <- m$stateNames
    names(timeInStates$upper) <- m$stateNames
    
    stateProbs <- list(est=xbar$stateProbs,se=MI_se$stateProbs,lower=lower$stateProbs,upper=upper$stateProbs)
    rownames(stateProbs$est) <- m$data$ID
    rownames(stateProbs$se) <- m$data$ID
    rownames(stateProbs$lower) <- m$data$ID
    rownames(stateProbs$upper) <- m$data$ID
    colnames(stateProbs$est) <- m$stateNames
    colnames(stateProbs$se) <- m$stateNames
    colnames(stateProbs$lower) <- m$stateNames
    colnames(stateProbs$upper) <- m$stateNames
  }
  
  mh <- im[[1]]
  attr(mh,"class") <- NULL
  mh$mle <- NULL
  mh$mod <- NULL
  mh$CI_real <- NULL
  mh$CI_beta <- NULL
  mh$data$step<-apply(matrix(unlist(lapply(im,function(x) x$data$step)),ncol=length(mh$data$step),byrow=TRUE),2,mean)
  mh$data$angle<-apply(matrix(unlist(lapply(im,function(x) x$data$angle)),ncol=length(mh$data$angle),byrow=TRUE),2,CircStats::circ.mean)
  mh$data$x<-apply(matrix(unlist(lapply(im,function(x) x$data$x)),ncol=length(mh$data$x),byrow=TRUE),2,mean)
  mh$data$y<-apply(matrix(unlist(lapply(im,function(x) x$data$y)),ncol=length(mh$data$y),byrow=TRUE),2,mean)
  mh$data$ice<-apply(matrix(unlist(lapply(im,function(x) x$data$ice)),ncol=length(mh$data$ice),byrow=TRUE),2,mean)
  
  if(!is.null(m$mle$gamma)){
    gammaPar <- list(est=xbar$gamma,se=MI_se$gamma,lower=lower$gamma,upper=upper$gamma)
    rownames(gammaPar$est) <- m$stateNames
    rownames(gammaPar$se) <- m$stateNames
    rownames(gammaPar$lower) <- m$stateNames
    rownames(gammaPar$upper) <- m$stateNames
    colnames(gammaPar$est) <- m$stateNames
    colnames(gammaPar$se) <- m$stateNames
    colnames(gammaPar$lower) <- m$stateNames
    colnames(gammaPar$upper) <- m$stateNames
    mh$Par<-list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar,beta=beta,gamma=gammaPar,timeInStates=timeInStates,states=states,stateProbs=stateProbs)
  } else mh$Par<-list(stepPar=stepPar,anglePar=anglePar,omegaPar=omegaPar,dryPar=dryPar,divePar=divePar,icePar=icePar,landPar=landPar,beta=beta,timeInStates=timeInStates,states=states,stateProbs=stateProbs)
  
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