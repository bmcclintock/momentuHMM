getDM<-function(data,DM,dist,nbStates,estAngleMean){
  
  distnames<-names(dist)
  fullDM <- vector('list',length(dist))
  names(fullDM) <- distnames
  nbObs<-nrow(data)
  
  for(i in distnames){
    switch(dist[[i]],
           "beta"={
             parNames<-c("shape1","shape2")
           },
           "pois"={
             parNames<-"lambda"
           },
           "weibull"={
             parNames<-c("shape","scale")
           },
           "gamma"={
             parNames<-c("mean","sd")
           },
           "lnorm"={
             parNames <- c("location","scale")
           },
           "exp"={
             parNames <- c("rate")
           },
           "vm"={
             if(estAngleMean[[i]]) { # if the angle mean is estimated
               parNames <- c("mean","sd") 
             }
             else {
               parNames <- c("sd")
             }
           },
           "wrpcauchy"={
             if(estAngleMean[[i]]) {
               parNames <- c("mean","concentration")
             }
             else {
               parNames <- c("concentration")
             }
           }
    )
    if(is.list(DM[[i]])){
      if(!all(names(DM[[i]]) %in% parNames)) stop('DM for ',i,' must include formula for ',paste(parNames,collapse=" and "))
      parSize<-unlist(lapply(DM[[i]],function(x) length(attr(terms.formula(x),"term.labels"))))+1
      tmpDM<-matrix(0,length(parNames)*nbStates*nbObs,sum(parSize)*nbStates)
      parInd<-0
      for(j in 1:length(parNames)){
        for(state in 1:nbStates){
          tmpDM[(state-1)*nbObs+(j-1)*nbStates*nbObs+1:nbObs,parInd*nbStates+state]<-1
          if(parSize[j]>1){
            for(jj in 2:parSize[j]){
              tmpDM[(state-1)*nbObs+(j-1)*nbStates*nbObs+1:nbObs,parInd*nbStates+(jj-1)*nbStates+state]<-model.matrix(DM[[i]][[parNames[j]]],data)[,jj]
              print(c(j,state,jj,range((state-1)*nbObs+(j-1)*nbStates*nbObs+1:nbObs),parInd*nbStates+(jj-1)*nbStates+state))
            }
          }
        }
        parInd<-parSize[j]
      }
    } else {
      tmpDM<-matrix(0,ncol(DM[[i]]),length(parNames)*nbStates*nbObs)
      DMterms<-unique(DM[[i]][suppressWarnings(which(is.na(as.numeric(DM[[i]]))))])
      for(j in 1:length(parNames)){
        for(state in 1:nbStates){
          tmpDM[,(state-1)*nbObs+(j-1)*nbStates*nbObs+1:nbObs]<-DM[[i]][(j-1)*nbStates+state,]
        }
      }
      tmpDM<-t(tmpDM)
      for(cov in DMterms){
        covs<-model.matrix(formula(paste("~",cov)),data)[,2]
        for(j in 1:length(parNames)){
          for(state in 1:nbStates){
            for(k in 1:nbObs){
              ind<-which(tmpDM[(state-1)*nbObs+(j-1)*nbStates*nbObs+k,]==cov)
              tmpDM[(state-1)*nbObs+(j-1)*nbStates*nbObs+k,][ind]<-covs[k]
            }
          }
        }
      }
      tmpDM<-matrix(as.numeric(tmpDM),nrow(tmpDM),ncol(tmpDM))
    }
    fullDM[[i]]<-tmpDM
  }
  fullDM[distnames]
}