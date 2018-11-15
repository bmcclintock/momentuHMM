get_inflation <- function(data,dist,DM,Par0,nbStates){
  distnames <- names(dist)
  zeroInflation <- oneInflation <- vector('list',length(distnames))
  names(zeroInflation) <- names(oneInflation) <- distnames
  for(i in distnames){
    if(dist[[i]] %in% zeroInflationdists){
      if(length(which(data[[i]]==0))>0) {
        zeroInflation[[i]]<-TRUE
      }
      else 
        zeroInflation[[i]]<-FALSE
    }
    else zeroInflation[[i]]<-FALSE
    if(dist[[i]] %in% oneInflationdists){
      if(length(which(data[[i]]==1))>0) {
        oneInflation[[i]]<-TRUE
      }
      else 
        oneInflation[[i]]<-FALSE
    }
    else oneInflation[[i]]<-FALSE
    if(is.null(DM[[i]])){
      if(zeroInflation[[i]]){
        # check that zero-mass is in the open interval (0,1)
        zm0 <- Par0[[i]][(length(Par0[[i]])-nbStates*oneInflation[[i]]-nbStates+1):(length(Par0[[i]])-nbStates*oneInflation[[i]])]
        zm0[which(zm0==0)] <- 1e-8
        zm0[which(zm0==1)] <- 1-1e-8
        Par0[[i]][(length(Par0[[i]])-nbStates*oneInflation[[i]]-nbStates+1):(length(Par0[[i]])-nbStates*oneInflation[[i]])] <- zm0
      }
      if(oneInflation[[i]]){
        # check that one-mass is in the open interval (0,1)
        om0 <- Par0[[i]][(length(Par0[[i]])-nbStates+1):length(Par0[[i]])]
        om0[which(om0==0)] <- 1e-8
        om0[which(om0==1)] <- 1-1e-8
        Par0[[i]][(length(Par0[[i]])-nbStates+1):length(Par0[[i]])] <- om0
      }
    }
  }
  list(Par0=Par0,zeroInflation=zeroInflation,oneInflation=oneInflation)
}