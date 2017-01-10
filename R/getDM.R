getDM<-function(data,DM,dist,nbStates,parNames,bounds){
  
  distnames<-names(dist)
  fullDM <- vector('list',length(dist))
  names(fullDM) <- distnames
  fullbounds <- vector('list',length(dist))
  names(fullbounds) <- distnames
  nbObs<-nrow(data)
  parSize<-lapply(parNames,length)
  
  is.formula <- function(x)
    tryCatch(inherits(x,"formula"),error= function(e) {FALSE})
  
  for(i in distnames){
    if(is.list(DM[[i]])){
      if(!all(names(DM[[i]]) %in% parNames[[i]]) | !all(unlist(lapply(DM[[i]],is.formula)))) stop('DM for ',i,' must include formula for ',paste(parNames[[i]],collapse=" and "))
      parSizeDM<-unlist(lapply(DM[[i]],function(x) length(attr(terms.formula(x),"term.labels"))))+1
      tmpDM<-array(0,dim=c(parSize[[i]]*nbStates,sum(parSizeDM)*nbStates,nbObs))
      DMnames<-character(sum(parSizeDM)*nbStates)
      #tmpbounds<-matrix(0,2,length(parNames[[i]])*nbStates*nbObs)
      parInd<-0
      for(j in 1:length(parNames[[i]])){
        tmpCov<-model.matrix(DM[[i]][[parNames[[i]][j]]],data)
        for(state in 1:nbStates){
          tmpDM[(j-1)*nbStates+state,(state-1)*parSizeDM[j]+parInd*nbStates+1:parSizeDM[j],]<-t(tmpCov)
          DMnames[(state-1)*parSizeDM[j]+parInd*nbStates+1:parSizeDM[j]]<-paste0(parNames[[i]][j],state,":",colnames(tmpCov))
        }
        parInd<-parSizeDM[j]
      }
      #tmpbounds[,(j-1)*nbObs+(state-1)*length(parNames[[i]])*nbObs+1:nbObs]<-bounds[[i]][(j-1)*nbStates+state,]
      #tmpbounds <- t(tmpbounds)
    } else {
      if(is.null(dim(DM[[i]]))) stop("DM for ",i," is not specified correctly")
      tmpDM<-array(DM[[i]],dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs))
      DMnames<-colnames(DM[[i]])
      if(is.null(DMnames)) warning("No names for the regression coeffs were provided in DM for ",i)
      #tmpbounds<-array(bounds[[i]],dim=c(parSize[[i]]*nbStates,2,nbObs))
      DMterms<-unique(DM[[i]][suppressWarnings(which(is.na(as.numeric(DM[[i]]))))])
      for(cov in DMterms){
        covs<-model.matrix(formula(paste("~",cov)),data)[,2]
        for(k in 1:nbObs){
          ind<-which(tmpDM[,,k]==cov)
          tmpDM[,,k][ind]<-covs[k]
        }
      }
      tmpDM<-array(as.numeric(tmpDM),dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs))
    }
    colnames(tmpDM)<-DMnames
    fullDM[[i]]<-tmpDM
    #fullbounds[[i]]<-tmpbounds
  }
  tmp<-simpDM<-lapply(fullDM,function(x) apply(x,1:2,unique))
  for(i in distnames){
    for(j in which(mapply(length,tmp[[i]])>1 & mapply(length,tmp[[i]])<nbObs)){
      simpDM[[i]][[j]]<-fullDM[[i]][ceiling(j/ncol(tmp[[i]])),ceiling(j/nrow(tmp[[i]])),]
    }
  }
  simpDM[distnames]
}