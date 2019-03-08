hierGamma <- function(m){
  
  if(is.miSum(m) && inherits(m,"hierarchical")){
    m$mle <- lapply(m$Par$real,function(x) x$est)
    m$mle$beta <- m$Par$beta$beta$est
    m$mle$pi <- m$Par$real$pi$est
    m$mle$delta <- m$Par$real$delta$est
    m$mod <- list()
    if(!is.null(m$conditions$recharge)){
      nbRecovs <- ncol(m$g0covs) + ncol(m$reCovs)
      m$mle$g0 <- c(m$Par$beta$g0$est)
      names(m$mle$g0) <- colnames(m$Par$beta$g0$est)
      m$mle$theta <- c(m$Par$beta$theta$est)
      names(m$mle$theta) <- colnames(m$Par$beta$theta$est)
    } else nbRecovs <- 0
    m$mod$estimate <- expandPar(m$MIcombine$coefficients,m$conditions$optInd,unlist(m$conditions$fixPar),m$conditions$wparIndex,m$conditions$betaCons,m$conditions$deltaCons,length(m$stateNames),ncol(m$covsDelta)-1,m$conditions$stationary,nrow(m$Par$beta$beta$est)/m$conditions$mixtures-1,nbRecovs,m$conditions$mixtures,ncol(m$covsPi)-1)
    m$mod$hessian <- NA
    m$mod$Sigma <- matrix(0,length(m$mod$estimate),length(m$mod$estimate))
    m$mod$Sigma[-m$conditions$optInd,-m$conditions$optInd] <- m$MIcombine$variance
    m$CIreal <- m$Par$real
    m <- momentuHMM(m)
  } else if(!inherits(m,"momentuHierHMM")) stop("m must be a 'momentuHierHMM' object")
  
  CIgamma <- data.tree::Node$new("hierGamma")
  CIdelta <- data.tree::Node$new("hierDelta")
  
  hierStates <- m$conditions$hierStates
  
  mixtures <- m$conditions$mixtures
  
  ref <- hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2)
  mixref <- rep(seq(0,mixtures*length(m$stateNames)-1,length(m$stateNames)),each=length(ref))+ref
  if(mixtures>1) nameref <- paste0(rep(names(ref),mixtures),"_mix",rep(1:mixtures,each=length(ref)))
  else nameref <- names(ref)
  
  CIgamma$AddChild("level1",gamma=lapply(m$CIreal$gamma,function(x) matrix(x[mixref,ref],length(mixref),length(ref),dimnames=list(nameref,names(ref)))))
  CIdelta$AddChild("level1",delta=lapply(m$CIreal$delta,function(x) matrix(x[,ref],nrow(x),dimnames=list(rownames(x),names(ref)))))
  
  for(j in 2:(hierStates$height-1)){
    
    t <- data.tree::Traverse(hierStates,filterFun=function(x) x$level==j)
    names(t) <- hierStates$Get("name",filterFun=function(x) x$level==j)
    
    tmpGamma <- CIreal(m,covs=data.frame(level=j),parms=c("gamma"))$gamma
    tmpDelta <- CIreal(m,covs=data.frame(level=paste0(j,"i")),parms=c("gamma"))$gamma
    CIgamma$AddChild(paste0("level",j),gamma=list())
    CIdelta$AddChild(paste0("level",j),delta=list())
    
    ref <- hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j)
    
    for(k in names(t)){
      levelStates <- t[[k]]$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==j+1)#t[[k]]$Get("state",filterFun = data.tree::isLeaf)
      if(!is.null(levelStates)){
        mixref <- rep(seq(0,mixtures*length(m$stateNames)-1,length(m$stateNames)),each=length(levelStates))+levelStates
        if(mixtures>1) nameref <- paste0(rep(names(levelStates),mixtures),"_mix",rep(1:mixtures,each=length(levelStates)))
        else nameref <- names(levelStates)
        
        mixrefk <- rep(seq(0,mixtures*length(m$stateNames)-1,length(m$stateNames)),each=length(ref[[k]]))+ref[[k]]
        if(mixtures>1) namerefk <- paste0(rep(names(ref[k]),mixtures),"_mix",1:mixtures)
        else namerefk <- k
        
        CIgamma[[paste0("level",j)]]$gamma[[k]] <- lapply(tmpGamma,function(x) matrix(x[mixref,levelStates],length(mixref),length(levelStates),dimnames=list(nameref,names(levelStates))))
        CIdelta[[paste0("level",j)]]$delta[[k]] <- lapply(tmpDelta,function(x) matrix(x[mixrefk,levelStates],mixtures,dimnames=list(namerefk,names(levelStates))))
      }
    }  
  }
  list(hierDelta=CIdelta,hierGamma=CIgamma)
}
  