getDM<-function(data,DM,dist,nbStates,parNames,bounds,Par,cons,workcons,zeroInflation){
  
  distnames<-names(dist)
  fullDM <- vector('list',length(dist))
  names(fullDM) <- distnames
  fullbounds <- vector('list',length(dist))
  names(fullbounds) <- distnames
  nbObs<-nrow(data)
  parSize<-lapply(parNames,length)
  
  for(i in distnames){
    if(is.null(DM[[i]])){
      tmpDM <- diag(length(Par[[i]]))
      tmpDM <- array(tmpDM,dim=c(nrow(tmpDM),ncol(tmpDM),nbObs))
      DMnames <- paste0(rep(parNames[[i]],each=nbStates),"_",1:nbStates,":(Intercept)")
    } else if(is.list(DM[[i]])){
      if(!all(parNames[[i]] %in% names(DM[[i]])) | !all(unlist(lapply(DM[[i]],is.formula)))) stop('DM for ',i,' must include formula for ',paste(parNames[[i]],collapse=" and "))
      DM[[i]]<-DM[[i]][parNames[[i]]]
      if(any(unlist(lapply(DM[[i]],function(x) attr(terms(x),"response")!=0))))
        stop("The response variable should not be specified in the DM formula for ",i)
      parSizeDM<-unlist(lapply(DM[[i]],function(x) length(attr(terms.formula(x),"term.labels"))))+1
      tmpDM<-array(0,dim=c(parSize[[i]]*nbStates,sum(parSizeDM)*nbStates,nbObs))
      DMnames<-character(sum(parSizeDM)*nbStates)
      parInd<-0
      for(j in 1:length(parNames[[i]])){
        tmpCov<-model.matrix(DM[[i]][[parNames[[i]][j]]],data)
        if(nrow(tmpCov)!=nbObs) stop("covariates cannot contain missing values")
        for(state in 1:nbStates){
          tmpDM[(j-1)*nbStates+state,(state-1)*parSizeDM[j]+parInd*nbStates+1:parSizeDM[j],]<-t(tmpCov)
          DMnames[(state-1)*parSizeDM[j]+parInd*nbStates+1:parSizeDM[j]]<-paste0(parNames[[i]][j],"_",state,":",colnames(tmpCov))
        }
        parInd<-sum(parSizeDM[1:j])
      }
      if(ncol(tmpDM)!=length(Par[[i]])) stop("Based on DM$",i,", Par$",i," must be of length ",ncol(tmpDM))
    } else {
      if(is.null(dim(DM[[i]]))) stop("DM for ",i," is not specified correctly")
      tmpDM<-array(DM[[i]],dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs))
      DMnames<-colnames(DM[[i]])
      DMterms<-unique(DM[[i]][suppressWarnings(which(is.na(as.numeric(DM[[i]]))))])
      for(cov in DMterms){
        covs<-model.matrix(formula(paste("~",cov)),data)[,2]
        if(length(covs)!=nbObs) stop("covariates cannot contain missing values")
        for(k in 1:nbObs){
          ind<-which(tmpDM[,,k]==cov)
          tmpDM[,,k][ind]<-covs[k]
        }
      }
      tmpDM<-array(as.numeric(tmpDM),dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs))
    }
    colnames(tmpDM)<-DMnames
    fullDM[[i]]<-tmpDM
  }
  tmp<-simpDM<-lapply(fullDM,function(x) apply(x,1:2,unique))
  for(i in distnames){
    k <- which(matrix(mapply(length,tmp[[i]])>1 & mapply(length,tmp[[i]])<nbObs,nrow(tmp[[i]]),ncol(tmp[[i]])),arr.ind=TRUE)
    if(length(k)){
      for(j in 1:nrow(k)){
        simpDM[[i]][[k[j,1],k[j,2]]]<-fullDM[[i]][k[j,1],k[j,2],]
      }
    }
  }
  simpDM<-simpDM[distnames]
  
  for(i in distnames){
    if(nrow(simpDM[[i]])!=(parSize[[i]]*nbStates)){
      error<- paste0("DM for ",i," should have ",(parSize[[i]]*nbStates)," rows")
      if(zeroInflation[[i]])
        stop(paste0(error,". Should zero inflation parameters be included?"))
      else stop(error)
    }
  }
  
  if(any(unlist(lapply(Par,length))!=unlist(lapply(simpDM,ncol))))
    stop("Dimension mismatch between Par and DM for: ",paste(names(which(unlist(lapply(Par,length))!=unlist(lapply(simpDM,ncol)))),collapse=", "))
  
  if(sum((unlist(parSize)>0)*unlist(lapply(simpDM,ncol)))!=length(unlist(Par))) {
    error <- "Wrong number of initial parameters"
    stop(error)
  }
  
  if(is.null(cons)){
    cons <- vector('list',length(distnames))
    names(cons) <- distnames
  } else {
    if(!is.list(cons) | is.null(names(cons))) stop("'cons' must be a named list")
  }
  for(i in distnames){
    if(is.null(cons[[i]])) cons[[i]] <- rep(1,ncol(simpDM[[i]]))
  }
  cons<-cons[distnames]
  if(any(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and cons for: ",paste(names(which(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  
  if(is.null(workcons)){
    workcons <- vector('list',length(distnames))
    names(workcons) <- distnames
  } else {
    if(!is.list(workcons) | is.null(names(workcons))) stop("'workcons' must be a named list")
  }
  for(i in distnames){
    if(is.null(workcons[[i]])) workcons[[i]] <- rep(0,ncol(simpDM[[i]]))
  }
  for(i in which(!(dist %in% "wrpcauchy"))){
    workcons[[distnames[i]]]<-rep(0,ncol(simpDM[[distnames[i]]]))
  }
  workcons<-workcons[distnames]
  if(any(unlist(lapply(workcons,length))!=unlist(lapply(Par,length)))) 
    stop("Length mismatch between Par and workcons for: ",paste(names(which(unlist(lapply(workcons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  
  DMind <- lapply(simpDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))
  
  return(list(fullDM=simpDM,DMind=DMind,cons=cons,workcons=workcons))
}