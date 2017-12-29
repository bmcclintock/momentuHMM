checkInputs<-function(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,cons,workcons,stateNames,checkInflation=FALSE)
{
  distnames<-names(dist)
  
  for(i in distnames){
    dist[[i]]<-match.arg(dist[[i]],momentuHMMdists)
  }

  if(nbStates<0)
    stop("nbStates should be at least 1.")
  
  if(is.null(estAngleMean)){
    estAngleMean <- vector('list',length(distnames))
    names(estAngleMean) <- distnames
  } else {
    if(!is.list(estAngleMean) | is.null(names(estAngleMean))) stop("'estAngleMean' must be a named list")
  }
  for(i in distnames[which(!(dist %in% angledists))]){
    estAngleMean[[i]] <- FALSE
  }
  for(i in distnames){
    if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
    if(!is.logical(estAngleMean[[i]])) stop("estAngleMean$",i," must be logical")
  }
  estAngleMean<-estAngleMean[distnames]
  
  if(is.null(circularAngleMean)){
    circularAngleMean <- vector('list',length(distnames))
    names(circularAngleMean) <- distnames
  } else {
    if(!is.list(circularAngleMean) | is.null(names(circularAngleMean))) stop("'circularAngleMean' must be a named list")
  }

  consensus <- vector('list',length(distnames))
  names(consensus) <- distnames
  
  ndist <- dist
  for(i in distnames){
    if(dist[[i]]=="vmConsensus"){
      consensus[[i]] <- TRUE
      estAngleMean[[i]] <- circularAngleMean[[i]] <- TRUE
      if(is.null(DM[[i]])) stop("DM$",i," must be specified when dist$",i,"=vmConsensus")
      ndist[[i]] <- gsub("Consensus","",dist[[i]])
    } else consensus[[i]] <- FALSE
    if(is.null(circularAngleMean[[i]]) | !estAngleMean[[i]]) circularAngleMean[[i]] <- FALSE
    if(!is.logical(circularAngleMean[[i]])) stop("circularAngleMean$",i," must be logical")
    if(circularAngleMean[[i]] & is.null(DM[[i]])) stop("DM$",i," must be specified when circularAngleMean$",i,"=TRUE")
  }
  circularAngleMean<-circularAngleMean[distnames]
  consensus<-consensus[distnames]
  
  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)

  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,oneInflation,DM,userBounds)
  parSize <- p$parSize
  
  bounds <- p$bounds

  if(is.null(DM)){
    DM <- cons <- workcons <- vector('list',length(distnames))
    names(DM) <- names(cons) <- names(workcons) <- distnames
  } else {
    if(!is.list(DM) | is.null(names(DM))) stop("'DM' must be a named list")
    if(!any(names(DM) %in% distnames)) stop("DM names must include at least one of: ",paste0(distnames,collapse=", "))
  }
  
  for(i in distnames){
    if(is.null(DM[[i]]) & length(Par[[i]])!=(parSize[[i]]*nbStates)){
      error<-paste0("Wrong number of initial parameters -- there should be ",parSize[[i]]*nbStates," initial ",i," parameters")
      if(zeroInflation[[i]]) error<-paste0(error," -- zero-mass parameters should be included")
      if(oneInflation[[i]]) error<-paste0(error," -- one-mass parameters should be included")
      stop(error)
    }
    if(checkInflation & zeroInflation[[i]] & oneInflation[[i]]){
      ztmp <- p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates*oneInflation[[i]]-nbStates:1+1,,drop=FALSE]
      if(any(ztmp[,1]!=0 | ztmp[,2]!=1)) warning("when both zero- and one-inflation are included, then userBounds for zeromass parameters are ignored")
      otmp <- p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates:1+1,,drop=FALSE]
      if(any(otmp[,1]!=0 | otmp[,2]!=1)) warning("when both zero- and one-inflation are included, then userBounds for onemass parameters are ignored")
      p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates*oneInflation[[i]]-nbStates:1+1,1] <- p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates:1+1,1] <- 0
      p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates*oneInflation[[i]]-nbStates:1+1,2] <- p$bounds[[i]][(parSize[[i]] * nbStates)-nbStates:1+1,2] <- 1
    }
    #if(!is.null(DM[[i]]) !is.null(userBounds[[i]])) stop("either userBounds$",i," or DM$",i," must be NULL")
  }

  return(list(p=p,dist=ndist,estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,consensus=consensus,DM=DM,cons=cons,workcons=workcons))
}