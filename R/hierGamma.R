hierGamma <- function(m){
  
  if(!inherits(m,"momentuHierHMM")) stop("m must be a 'momentuHierHMM' object")
  
  CIgamma <- data.tree::Node$new("hierGamma")
  CIdelta <- data.tree::Node$new("hierDelta")
  
  hierStates <- m$conditions$hierStates
  
  ref1 <- hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2)
  
  CIgamma$AddChild("level1",gamma=lapply(m$CIreal$gamma,function(x) matrix(x[ref1,ref1],length(ref1),length(ref1),dimnames=list(names(ref1),names(ref1)))))
  CIdelta$AddChild("level1",delta=lapply(m$CIreal$delta,function(x) matrix(x[,ref1],nrow(x),dimnames=list(rownames(x),names(ref1)))))
  
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
        CIgamma[[paste0("level",j)]]$gamma[[k]] <- lapply(tmpGamma,function(x) matrix(x[levelStates,levelStates],length(levelStates),length(levelStates),dimnames=list(names(levelStates),names(levelStates))))
        CIdelta[[paste0("level",j)]]$delta[[k]] <- lapply(tmpDelta,function(x) matrix(x[ref[[k]],levelStates],1,dimnames=list(k,names(levelStates))))
      }
    }  
  }
  list(hierDelta=CIdelta,hierGamma=CIgamma)
}
  