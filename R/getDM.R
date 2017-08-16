#' @importFrom stats formula terms.formula
getDM<-function(data,DM,dist,nbStates,parNames,bounds,Par,cons,workcons,zeroInflation,oneInflation,circularAngleMean,ParChecks=TRUE){
  
  distnames<-names(dist)
  fullDM <- vector('list',length(dist))
  names(fullDM) <- distnames
  fullbounds <- vector('list',length(dist))
  names(fullbounds) <- distnames
  nbObs<-nrow(data)
  parSize<-lapply(parNames,length)
  
  for(i in distnames){
    if(is.null(DM[[i]])){
      tmpDM <- diag(parSize[[i]]*nbStates)
      tmpDM <- array(tmpDM,dim=c(nrow(tmpDM),ncol(tmpDM),nbObs))
      DMnames <- paste0(rep(parNames[[i]],each=nbStates),"_",1:nbStates,":(Intercept)")
    } else if(is.list(DM[[i]])){
      if(!all(parNames[[i]] %in% names(DM[[i]])) | !all(unlist(lapply(DM[[i]],is.formula)))) stop('DM$',i,' must include formula for ',paste(parNames[[i]],collapse=" and "))
      if(!all(names(DM[[i]]) %in% parNames[[i]])){
        err <- paste0('DM$',i,' should only include formula for ',paste(parNames[[i]],collapse=" and ")," parameter(s)")
        if(dist[[i]] %in% angledists) stop(err,". Check DM$",i,", estAngleMean$",i,", and circularAngleMean$",i)
      }
      if(!zeroInflation[[i]] & "zeromass" %in% names(DM[[i]])) stop('zeromass should not be included in DM$',i)
      if(!oneInflation[[i]] & "onemass" %in% names(DM[[i]])) stop('onemass should not be included in DM$',i)
      DM[[i]]<-DM[[i]][parNames[[i]]]
      if(any(unlist(lapply(DM[[i]],function(x) attr(terms(x),"response")!=0))))
        stop("The response variable should not be specified in the DM formula for ",i)
      
      formulaStates <- vector('list',length(parNames[[i]]))
      names(formulaStates) <- parNames[[i]]
      for(j in parNames[[i]])
        formulaStates[[j]]<- stateFormulas(DM[[i]][[j]],nbStates)
      
      #if(circularAngleMean[[i]]){
      tmpCov <- vector('list',length(parNames[[i]]))
      names(tmpCov) <- parNames[[i]]
      for(j in names(DM[[i]])){
        tmpCov[[j]] <- vector('list',nbStates)
        for(state in 1:nbStates){
          tmpCov[[j]][[state]]<-model.matrix(formulaStates[[j]][[state]],data)
          if(circularAngleMean[[i]]){
            if(j=="mean" & attr(terms.formula(formulaStates[[j]][[state]]),"intercept")) tmpCov[[j]][[state]] <- tmpCov[[j]][[state]][,-1,drop=FALSE]
            #if(!length(tmpCov[[j]][[state]])) stop("invalid circular-circular regression formula for ",i," ",j)
          }
        }
      }
      #} else tmpCov<-lapply(DM[[i]],function(x) model.matrix(x,data))
      parSizeDM<-unlist(lapply(tmpCov,function(x) lapply(x,ncol)))
      tmpDM<-array(0,dim=c(parSize[[i]]*nbStates,sum(parSizeDM),nbObs))
      DMnames<-character(sum(parSizeDM))
      parInd<-0
      for(j in 1:length(parNames[[i]])){
        parmStates<-which(unlist(lapply(tmpCov[[j]],length))>0)
        for(state in parmStates){
          if(nrow(tmpCov[[j]][[state]])!=nbObs) stop("covariates cannot contain missing values")
          tmpDM[(j-1)*nbStates+state,parInd+1:parSizeDM[(j-1)*nbStates+state],]<-t(tmpCov[[j]][[state]])
          DMnames[parInd+1:parSizeDM[(j-1)*nbStates+state]]<-paste0(parNames[[i]][j],"_",state,":",colnames(tmpCov[[j]][[state]]))
          parInd<-sum(parSizeDM[1:((j-1)*nbStates+state)])
        }
      }
      if(ncol(tmpDM)!=length(Par[[i]]) & ParChecks) stop("Based on DM$",i,", Par$",i," must be of length ",ncol(tmpDM))
    } else {
      if(is.null(dim(DM[[i]]))) stop("DM for ",i," is not specified correctly")
      tmpDM<-suppressWarnings(array(as.numeric(DM[[i]]),dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs)))
      DMnames<-colnames(DM[[i]])
      if(is.null(DMnames)) DMnames<-paste0(i,"Beta",1:ncol(DM[[i]]))
      DMterms<-unique(DM[[i]][suppressWarnings(which(is.na(as.numeric(DM[[i]]))))])
      factorterms<-names(data)[unlist(lapply(data,is.factor))]
      factorcovs<-paste0(rep(factorterms,times=unlist(lapply(data[factorterms],nlevels))),unlist(lapply(data[factorterms],levels)))
      covs<-numeric()
      for(cov in DMterms){
        if(is.factor(data[[cov]])) stop('factor levels must be specified individually when using pseudo-design matrices')
        form<-formula(paste("~",cov))
        varform<-all.vars(form)
        if(any(varform %in% factorcovs)){
          factorvar<-factorcovs %in% varform
          tmpcov<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
          tmpcov<-gsub(factorcovs[factorvar],tmpcov,cov)
          tmpcovs<-model.matrix(formula(paste("~ 0 + ",tmpcov)),data)
          tmpcovs<-tmpcovs[,which(gsub(" ","",colnames(tmpcovs)) %in% gsub(" ","",cov))]
          covs<-cbind(covs,tmpcovs)
        } else {
          tmpcovs<-model.matrix(form,data)[,2]
          covs<-cbind(covs,tmpcovs)
        }
        if(length(tmpcovs)!=nbObs) stop("covariates cannot contain missing values")
      }
      if(length(DMterms)) tmpDM<-getDM_rcpp(tmpDM,covs,c(DM[[i]]),nrow(tmpDM),ncol(tmpDM),DMterms,nbObs)
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
      if(oneInflation[[i]])
        stop(paste0(error,". Should one inflation parameters be included?"))
      else stop(error)
    }
  }
  
  if(ParChecks){
    if(any(unlist(lapply(Par,length))!=unlist(lapply(simpDM,ncol))))
      stop("Dimension mismatch between Par and DM for: ",paste(names(which(unlist(lapply(Par,length))!=unlist(lapply(simpDM,ncol)))),collapse=", "))
    
    if(sum((unlist(parSize)>0)*unlist(lapply(simpDM,ncol)))!=length(unlist(Par))) {
      error <- "Wrong number of initial parameters"
      stop(error)
    }
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
  if(ParChecks){
    if(any(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))) 
      stop("Length mismatch between Par and cons for: ",paste(names(which(unlist(lapply(cons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  } else {
    if(any(unlist(lapply(cons,length))!=unlist(lapply(simpDM,ncol)))) 
      stop("Length mismatch between DM and cons for: ",paste(names(which(unlist(lapply(cons,length))!=unlist(lapply(simpDM,ncol)))),collapse=", "))    
  }
  if(is.null(workcons)){
    workcons <- vector('list',length(distnames))
    names(workcons) <- distnames
  } else {
    if(!is.list(workcons) | is.null(names(workcons))) stop("'workcons' must be a named list")
  }
  for(i in distnames){
    if(is.null(workcons[[i]])) workcons[[i]] <- rep(0,ncol(simpDM[[i]]))
  }
  #for(i in which(!(dist %in% "wrpcauchy"))){
  #  workcons[[distnames[i]]]<-rep(0,ncol(simpDM[[distnames[i]]]))
  #}
  workcons<-workcons[distnames]
  if(ParChecks){
    if(any(unlist(lapply(workcons,length))!=unlist(lapply(Par,length)))) 
      stop("Length mismatch between Par and workcons for: ",paste(names(which(unlist(lapply(workcons,length))!=unlist(lapply(Par,length)))),collapse=", "))
  } else {
    if(any(unlist(lapply(workcons,length))!=unlist(lapply(simpDM,ncol)))) 
      stop("Length mismatch between DM and workcons for: ",paste(names(which(unlist(lapply(workcons,length))!=unlist(lapply(simpDM,ncol)))),collapse=", "))      
  }
  DMind <- lapply(simpDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))
  
  for(i in distnames){
    #if(DMind[[i]]){
    #  getbndInd <- getboundInd(simpDM[[i]])
    #} else {
    #  getbndInd <- getboundInd(fullDM[[i]][,,1])
    #  #getbndInd <- apply(fullDM[[i]],3,getboundInd)
    #}
    getbndInd <- getboundInd(fullDM[[i]][,,1])
    bndInd <- which(!duplicated(getbndInd))
    if(any(bounds[[i]]!=bounds[[i]][bndInd,,drop=FALSE][getbndInd,,drop=FALSE])) stop('userBounds not consistent with DM for ',i)
    rownames(simpDM[[i]]) <- rownames(bounds[[i]])
  }
  
  return(list(fullDM=simpDM,DMind=DMind,cons=cons,workcons=workcons))
}