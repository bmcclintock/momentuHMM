hierGamma <- function(m){
  
  if(!inherits(m,"momentuHierHMM")) stop("m must be a 'momentuHierHMM' object")
  
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
        else nameref <- names(ref)
        
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
  