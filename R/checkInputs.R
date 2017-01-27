checkInputs<-function(nbStates,dist,Par,estAngleMean,zeroInflation,DM,userBounds,cons,workcons,stateNames)
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
  for(i in distnames){
    if(is.null(estAngleMean[[i]])) estAngleMean[[i]] <- FALSE
  }
  for(i in distnames[which(!(dist %in% angledists))]){
    estAngleMean[[i]] <- FALSE
  }
  estAngleMean<-estAngleMean[distnames]
  
  
  if(!is.null(stateNames) & length(stateNames)!=nbStates)
    stop("stateNames must have length ",nbStates)

  par0 <- unlist(Par)
  p <- parDef(dist,nbStates,estAngleMean,zeroInflation,DM,userBounds)
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
      stop(error)
    }
  }

  return(list(p=p,estAngleMean=estAngleMean,DM=DM,cons=cons,workcons=workcons))
}