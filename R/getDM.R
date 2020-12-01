#' @importFrom stats formula terms.formula
getDM<-function(data,DM,dist,nbStates,parNames,bounds,Par,zeroInflation,oneInflation,circularAngleMean,ParChecks=TRUE,wlag=FALSE){
  
  distnames<-names(dist)
  fullDM <- vector('list',length(dist))
  names(fullDM) <- distnames
  fullbounds <- vector('list',length(dist))
  names(fullbounds) <- distnames
  nbObs<-nrow(data)
  parSize<-lapply(parNames,length)
  #parCount <- vector('list',length(dist))
  lag <- 0
  lanInd <- vector('list',length(dist))
  names(lanInd) <- distnames
  
  for(i in distnames){
    if(is.null(DM[[i]])){
      if(dist[[i]] %in% rwdists) stop("DM$",i," must be specified")
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
      if(any(unlist(lapply(DM[[i]],function(x) attr(stats::terms(x),"response")!=0))))
        stop("The response variable should not be specified in the DM formula for ",i)
      
      formulaStates <- vector('list',length(parNames[[i]]))
      names(formulaStates) <- parNames[[i]]
      for(j in parNames[[i]]){
        formulaStates[[j]]<- stateFormulas(DM[[i]][[j]],nbStates,angleMean=(j=="mean" & !isFALSE(circularAngleMean[[i]])),data=data)
        if(dist[[i]] %in% rwdists){
          #if(dist[[i]] %in% mvndists){
          for(k in 1:nbStates){
            formterms<-attr(stats::terms.formula(formulaStates[[j]][[k]]),"term.labels")
            if(grepl("mean",j)) formulaStates[[j]][[k]] <- stats::as.formula(paste("~ + 0 + ",paste(c(paste0(i,gsub("mean","",j),"_tm1"),formterms),collapse=" + ")))
          }
          #}
        }
      }
      #if(circularAngleMean[[i]]){
      tmpCov <- vector('list',length(parNames[[i]]))
      names(tmpCov) <- parNames[[i]]
      for(j in names(DM[[i]])){
        tmpCov[[j]] <- vector('list',nbStates)
        for(state in 1:nbStates){
          if(wlag) lag <- get_crwlag(formulaStates[[j]][[state]],lag)
          tmpCov[[j]][[state]]<-stats::model.matrix(formulaStates[[j]][[state]],data)
          if(!isFALSE(circularAngleMean[[i]])){
            if(j=="mean"){
              if(attr(stats::terms.formula(formulaStates[[j]][[state]]),"intercept")) tmpCov[[j]][[state]] <- tmpCov[[j]][[state]][,-1,drop=FALSE]
              cnames <- colnames(tmpCov[[j]][[state]])
              # make sure columns are arranged in sin/cos pairs
              tmpCov[[j]][[state]] <- tmpCov[[j]][[state]][,cnames[order(match(gsub("cos","",gsub("sin","",cnames)),unique(gsub("cos","",gsub("sin","",cnames)))))],drop=FALSE]
            }
            #if(!length(tmpCov[[j]][[state]])) stop("invalid circular-circular regression formula for ",i," ",j)
          }
        }
      }
      #} else tmpCov<-lapply(DM[[i]],function(x) stats::model.matrix(x,data))
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
      for(j in DMnames){
        if(grepl("langevin\\(",j)) {
          if(!(dist[[i]] %in% rwdists)) stop("langevin() special formula can only be used for normal random walk data stream distributions")
          if(!grepl("mean",j))
            stop("langevin() special formula can only be used for the means of normal random walk data stream distributions")
          lanInd[[i]] <- c(lanInd[[i]],j)
        }
      }
      #parCount[[i]]<-ncol(tmpDM)
      #if(circularAngleMean[[i]]) parCount[[i]] <- parCount[[i]] - sum(parSizeDM[grepl("mean",names(parSizeDM))])/2
      #if(parCount[[i]]!=length(Par[[i]]) & ParChecks) stop("Based on DM$",i,", Par$",i," must be of length ",ncol(tmpDM))
    } else {
      if(is.null(dim(DM[[i]]))) stop("DM for ",i," is not specified correctly")
      if(nrow(DM[[i]])!=parSize[[i]]*nbStates) stop("DM$",i," should consist of ",parSize[[i]]*nbStates," rows")
      DMnames<-colnames(DM[[i]])
      if(is.null(DMnames)) DMnames<-paste0(i,"Beta",1:ncol(DM[[i]]))
      DMterms<-unique(DM[[i]][suppressWarnings(which(is.na(as.numeric(DM[[i]]))))])
      meanind <- NULL
      if(isFALSE(circularAngleMean[[i]])){
        tmpDM<-suppressWarnings(array(as.numeric(DM[[i]]),dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs)))
        newDM <- DM[[i]]
      } else {
        meanind<-which(!apply(DM[[i]][1:nbStates,,drop=FALSE],2,function(x) all(x==0)))
        if(length(meanind)){
          sdind<-which(!apply(DM[[i]][1:nbStates+nbStates,,drop=FALSE],2,function(x) all(x==0)))
          newDM <- matrix(0,nrow(DM[[i]]),ncol=length(meanind)*2+length(sdind))
          if(any(grepl("angleFormula",DM[[i]][,sdind]))) stop("angleFormula function only applies to angle mean")
          newDM[,length(meanind)*2+1:length(sdind)]<-DM[[i]][,sdind]
          colnames(newDM)[length(meanind)*2+1:length(sdind)] <- DMnames[sdind]
          colnames(newDM)[seq(1,length(meanind)*2-1,2)] <- DMnames[meanind]
          for(j in meanind){
            tmpcolname <- colnames(newDM)[seq(1,length(meanind)*2-1,2)[j]]
            for(jj in 1:nbStates){
              if(DM[[i]][jj,j]!=0){
                terms <- sort(attr(stats::terms(stateFormulas(formula(paste0("~",DM[[i]][jj,j])),nbStates,angleMean=TRUE,data=data)[[jj]]),"term.labels"),decreasing=TRUE)
                newDM[jj,seq(1,length(meanind)*2-1,2)[j]]<-terms[1]
                colnames(newDM)[seq(1,length(meanind)*2-1,2)[j]]<-paste0(tmpcolname,"sin")
                newDM[jj,seq(1,length(meanind)*2-1,2)[j]+1]<-terms[2]
                colnames(newDM)[seq(1,length(meanind)*2-1,2)[j]+1]<-paste0(tmpcolname,"cos")
              }
            }
          }
          tmpDM<-suppressWarnings(array(as.numeric(newDM),dim=c(nrow(newDM),ncol(newDM),nbObs)))
          DMnames<-colnames(newDM)
          DMterms<-unique(newDM[suppressWarnings(which(is.na(as.numeric(newDM))))])
        } else {
          tmpDM<-suppressWarnings(array(as.numeric(DM[[i]]),dim=c(nrow(DM[[i]]),ncol(DM[[i]]),nbObs)))
          newDM <- DM[[i]]
        }
      }
      if(!length(meanind) & any(grepl("angleFormula",DMterms))) stop("angleFormula function only applies to circular-circular regression model for angle mean")
      factorterms<-names(data)[unlist(lapply(data,is.factor))]
      factorcovs<-paste0(rep(factorterms,times=unlist(lapply(data[factorterms],nlevels))),unlist(lapply(data[factorterms],levels)))
      covs<-numeric()
      for(cov in DMterms){
        if(grepl("langevin\\(",cov)) {
          if(!(dist[[i]] %in% rwdists)) stop("langevin() special formula can only be used for normal random walk data stream distributions")
          if(!all(grepl("mean",rownames(bounds[[i]])[which(apply(DM[[i]],1,function(x) any(grepl("langevin\\(",x))))])))
            stop("langevin() special formula can only be used for the means of normal random walk data stream distributions")
          lanInd[[i]] <- c(lanInd[[i]],cov)
        }
        form<-stats::formula(paste("~",cov))
        varform<-all.vars(form)
        if(!all((varform %in% names(data)) | (varform %in% factorcovs))){
          varform <- varform[(varform %in% names(data)) | (varform %in% factorcovs)]
        }
        if(any(unlist(lapply(data[varform[!(varform %in% factorcovs)]],function(x) inherits(x,"factor"))))) stop('factor levels must be specified individually when using pseudo-design matrices')
        if(any(varform %in% factorcovs)){
          factorvar<-factorcovs %in% varform
          tmpcov<-rep(factorterms,times=unlist(lapply(data[factorterms],nlevels)))[which(factorvar)]
          tmpcovj <- cov
          for(j in 1:length(tmpcov)){
            tmpcovj <- gsub(factorcovs[factorvar][j],tmpcov[j],tmpcovj)
          }
          tform <- stats::formula(paste("~ 0 + ",tmpcovj))
          if(wlag) lag <- get_crwlag(tform,lag)
          tmpcovs<-stats::model.matrix(tform,data)
          tmpcovs<-tmpcovs[,which(gsub(" ","",colnames(tmpcovs)) %in% gsub(" ","",cov))]
          covs<-cbind(covs,tmpcovs)
        } else {
          if(wlag) lag <- get_crwlag(form,lag)
          tmpcovs<-stats::model.matrix(form,data)[,2]
          covs<-cbind(covs,tmpcovs)
        }
        if(length(tmpcovs)!=nbObs) stop("covariates cannot contain missing values")
      }
      if(length(DMterms)) tmpDM<-getDM_rcpp(tmpDM,covs,c(newDM),nrow(tmpDM),ncol(tmpDM),DMterms,nbObs)
      #parCount[[i]] <- ncol(tmpDM)
    }
    colnames(tmpDM)<-DMnames
    if(!wlag) fullDM[[i]]<-tmpDM
    else fullDM[[i]]<-tmpDM[,,dim(tmpDM)[3],drop=FALSE]
  }
  tmp<-simpDM<-lapply(fullDM,function(x) apply(x,1:2,unique))
  for(i in distnames){
    k <- which(matrix(mapply(length,tmp[[i]])>1 & mapply(length,tmp[[i]])<nbObs,nrow(tmp[[i]]),ncol(tmp[[i]])),arr.ind=TRUE)
    if(length(k)){
      for(j in 1:nrow(k)){
        simpDM[[i]][[k[j,1],k[j,2]]]<-fullDM[[i]][k[j,1],k[j,2],]
      }
    }
    if(!is.null(lanInd[[i]])){
      if(is.list(DM[[i]])){
        comment(simpDM[[i]]) <- paste0("langevin",which(grepl("langevin\\(",colnames(simpDM[[i]]))))
      } else {
        comment(simpDM[[i]]) <- paste0("langevin",which(apply(DM[[i]],2,function(x) any(grepl("langevin\\(",x)))))
      }
    }
  }
  simpDM<-simpDM[distnames]
  
  parCount<- lapply(simpDM,ncol)
  for(i in distnames[!unlist(lapply(circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(simpDM[[i]])))))
  }

  for(i in distnames){
    if(nrow(simpDM[[i]])!=(parSize[[i]]*nbStates)){
      error<- paste0("DM for ",i," should have ",(parSize[[i]]*nbStates)," rows")
      if(zeroInflation[[i]])
        stop(paste0(error,". Should zero inflation parameters be included?"))
      if(oneInflation[[i]])
        stop(paste0(error,". Should one inflation parameters be included?"))
      else stop(error)
    }
    #if(isFALSE(circularAngleMean[[i]])) parCount[[i]] <- ncol(simpDM[[i]])
  }
  parCount <- parCount[distnames]
  
  if(ParChecks){
    if(any(unlist(lapply(Par,length))!=unlist(parCount)))
      stop("Dimension mismatch between Par and DM for: ",paste(names(which(unlist(lapply(Par,length))!=unlist(parCount))),collapse=", "))
    
    if(sum((unlist(parSize)>0)*unlist(parCount))!=length(unlist(Par))) {
      error <- "Wrong number of initial parameters"
      stop(error)
    }
  }
  
  DMind <- lapply(simpDM,function(x) all(unlist(apply(x,1,function(y) lapply(y,length)))==1))
  
  for(i in distnames){
    #if(DMind[[i]]){
    #  getbndInd <- getboundInd(simpDM[[i]])
    #} else {
    #  getbndInd <- getboundInd(fullDM[[i]][,,1])
    #  #getbndInd <- apply(fullDM[[i]],3,getboundInd)
    #}
    getbndInd <- getboundInd(matrix(fullDM[[i]][,,1],dim(fullDM[[i]])[1],dim(fullDM[[i]])[2]))
    bndInd <- which(!duplicated(getbndInd))
    if(any(c(t(bounds[[i]]))!=c(t(bounds[[i]][bndInd,,drop=FALSE][getbndInd,,drop=FALSE])))) stop('userBounds not consistent with DM for ',i)
    rownames(simpDM[[i]]) <- rownames(bounds[[i]])
  }
  
  return(list(fullDM=simpDM,DMind=DMind,lag=lag))
}

get_crwlag <- function(form,clag){
  lag <- 0
  if(length(attr(stats::terms(form, specials="crw"),"specials")$crw)){
    Terms <- stats::terms(form,specials="crw")
    lag <- as.numeric(unlist(attr(prodlim::strip.terms(Terms[attr(Terms,"specials")$crw],specials="crw",arguments=list(crw=list("lag"=1))),"stripped.arguments")))
  }
  max(c(clag,lag))
}